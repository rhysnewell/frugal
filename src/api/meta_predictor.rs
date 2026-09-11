//! Reusable metagenomic gene predictor.
//!
//! Caches the 50 metagenomic models and uses a rayon thread pool
//! with a large stack for prediction.

use std::alloc::{alloc_zeroed, handle_alloc_error, Layout};
use std::os::raw::c_int;
use std::sync::Arc;

use rayon::prelude::*;
use rayon::ThreadPool;

use super::convert::gene_to_predicted;
use super::encode::SequenceBuffer;
use super::types::{PredictedGene, ProdigalConfig, ProdigalError};
use crate::types::{Gene, Training, MAX_SEQ, NUM_META};

use super::predict::validate_config;

/// Recommended stack size for Rayon pools used by `MetaPredictor`.
///
/// The low-level translated prediction routines use deep call stacks, so the
/// default Rust thread stack can be too small for worker threads.
pub const META_PREDICTOR_STACK_SIZE: usize = 32 * 1024 * 1024; // 32 MB

use crate::dprog::{dprog, eliminate_bad_genes};
use crate::gene::{add_genes, record_gene_data, tweak_final_starts};
use crate::node::{
    add_nodes, record_overlapping_starts, record_rbs_masks, reset_node_scores,
    score_nodes_with_rbs,
};

/// Reusable metagenomic gene predictor.
///
/// Pre-loads 50 metagenomic models and evaluates qualifying models with
/// Prodigal-compatible state preservation across adjacent model bins.
pub struct MetaPredictor {
    pool: Arc<ThreadPool>,
    models: Arc<Vec<Box<Training>>>,
    config: ProdigalConfig,
}

impl MetaPredictor {
    /// Create a new predictor with default config.
    pub fn new() -> Result<Self, ProdigalError> {
        Self::with_config(ProdigalConfig::default())
    }

    /// Create a new predictor with custom config.
    pub fn with_config(config: ProdigalConfig) -> Result<Self, ProdigalError> {
        validate_config(&config)?;

        let pool = rayon::ThreadPoolBuilder::new()
            .stack_size(META_PREDICTOR_STACK_SIZE)
            .build()
            .map_err(|e| {
                ProdigalError::Io(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    e.to_string(),
                ))
            })?;

        Self::with_config_and_thread_pool(config, Arc::new(pool))
    }

    /// Create a new predictor with default config and a shared Rayon thread pool.
    ///
    /// The supplied pool is reused for parallel metagenomic model evaluation.
    /// For large inputs, build the pool with at least `META_PREDICTOR_STACK_SIZE`
    /// stack bytes per worker thread.
    pub fn with_thread_pool(pool: Arc<ThreadPool>) -> Result<Self, ProdigalError> {
        Self::with_config_and_thread_pool(ProdigalConfig::default(), pool)
    }

    /// Create a new predictor with custom config and a shared Rayon thread pool.
    ///
    /// The supplied pool is reused for parallel metagenomic model evaluation.
    /// For large inputs, build the pool with at least `META_PREDICTOR_STACK_SIZE`
    /// stack bytes per worker thread.
    pub fn with_config_and_thread_pool(
        config: ProdigalConfig,
        pool: Arc<ThreadPool>,
    ) -> Result<Self, ProdigalError> {
        validate_config(&config)?;
        let models = Arc::new(load_meta_models());

        Ok(MetaPredictor {
            pool,
            models,
            config,
        })
    }

    /// Predict genes in the given sequence.
    pub fn predict(&self, seq: &[u8]) -> Result<Vec<PredictedGene>, ProdigalError> {
        validate_sequence(seq)?;

        self.pool.install(|| {
            let mut worker = Worker::new(&self.models);
            predict_parallel(seq, &mut worker.models, &self.config, &mut worker.buf)
        })
    }

    /// Predict genes in a batch of sequences, preserving input order.
    pub fn predict_batch<S>(&self, seqs: &[S]) -> Result<Vec<Vec<PredictedGene>>, ProdigalError>
    where
        S: AsRef<[u8]> + Sync,
    {
        for seq in seqs {
            validate_sequence(seq.as_ref())?;
        }

        let models = Arc::clone(&self.models);
        self.pool.install(move || {
            seqs.par_iter()
                .map_init(
                    || Worker::new(&models),
                    |worker, seq| {
                        predict_parallel(
                            seq.as_ref(),
                            &mut worker.models,
                            &self.config,
                            &mut worker.buf,
                        )
                    },
                )
                .collect()
        })
    }
}

/// Per-thread state. Each worker owns its own copy of the 50 models so the scoring path can
/// take them by mutable reference without copying one per candidate.
struct Worker {
    buf: SequenceBuffer,
    models: Vec<Box<Training>>,
}

impl Worker {
    fn new(models: &[Box<Training>]) -> Self {
        Worker {
            buf: SequenceBuffer::reusable(),
            models: models.to_vec(),
        }
    }
}

/// Reject empty or oversized input sequences before encoding.
fn validate_sequence(seq: &[u8]) -> Result<(), ProdigalError> {
    if seq.is_empty() {
        return Err(ProdigalError::EmptySequence);
    }
    if seq.len() > MAX_SEQ {
        return Err(ProdigalError::SequenceTooLong {
            length: seq.len(),
            max: MAX_SEQ,
        });
    }
    Ok(())
}

/// Allocate and initialize all `NUM_META` pre-trained metagenomic models on the heap.
fn load_meta_models() -> Vec<Box<Training>> {
    let mut models: Vec<Box<Training>> = Vec::with_capacity(NUM_META);
    for i in 0..NUM_META {
        unsafe {
            let layout = Layout::new::<Training>();
            let ptr = alloc_zeroed(layout) as *mut Training;
            if ptr.is_null() {
                handle_alloc_error(layout);
            }
            crate::training_data::load_metagenome(i, ptr);
            models.push(Box::from_raw(ptr));
        }
    }
    models
}

/// Run the metagenomic gene-prediction pipeline on a single sequence.
///
/// Encodes the input, rebuilds node arrays when the translation table changes,
/// scores all GC-compatible models, keeps the highest-scoring solution, and
/// converts its genes to `PredictedGene` values.
fn predict_parallel(
    seq: &[u8],
    models: &mut [Box<Training>],
    config: &ProdigalConfig,
    buf: &mut SequenceBuffer,
) -> Result<Vec<PredictedGene>, ProdigalError> {
    let closed = if config.closed_ends { 1 } else { 0 };

    let (slen, gc) = unsafe { buf.encode(seq, config.mask_n_runs) };
    if slen == 0 {
        return Err(ProdigalError::EmptySequence);
    }
    buf.ensure_node_capacity(slen);

    // GC window for model selection
    let mut low = 0.88495 * gc - 0.0102337;
    if low > 0.65 {
        low = 0.65;
    }
    let mut high = 0.86596 * gc + 0.1131991;
    if high < 0.35 {
        high = 0.35;
    }

    let mut best_score = f64::NEG_INFINITY;
    let mut best_nodes = Vec::new();
    let mut best_genes: Vec<Gene> = Vec::new();
    let mut best_tinf: Option<usize> = None;
    let mut nn: c_int = 0;
    let mut masks_fresh = false;
    let mut built_table: Option<c_int> = None;

    unsafe {
        for i in 0..NUM_META {
            if models[i].gc < low || models[i].gc > high {
                continue;
            }
            let tinf: &mut Training = &mut models[i];

            if built_table != Some(tinf.trans_table) {
                built_table = Some(tinf.trans_table);
                buf.clear_nodes(nn);
                nn = add_nodes(
                    buf.seq.as_mut_ptr(),
                    buf.rseq.as_mut_ptr(),
                    slen,
                    buf.nodes.as_mut_ptr(),
                    closed,
                    buf.masks.as_mut_ptr(),
                    buf.nmask,
                    tinf,
                );
                buf.nodes[..nn as usize]
                    .sort_unstable_by(|a, b| a.ndx.cmp(&b.ndx).then(b.strand.cmp(&a.strand)));
                masks_fresh = false;
            }

            if tinf.uses_sd == 1 && !masks_fresh {
                buf.ensure_rbs_capacity(nn);
                record_rbs_masks(
                    buf.seq.as_mut_ptr(),
                    buf.rseq.as_mut_ptr(),
                    slen,
                    buf.nodes.as_mut_ptr(),
                    nn,
                    buf.rbs_masks.as_mut_ptr(),
                );
                masks_fresh = true;
            }
            let rbs_masks = if masks_fresh {
                buf.rbs_masks.as_ptr()
            } else {
                std::ptr::null()
            };
            reset_node_scores(buf.nodes.as_mut_ptr(), nn);
            score_nodes_with_rbs(
                buf.seq.as_mut_ptr(),
                buf.rseq.as_mut_ptr(),
                slen,
                buf.nodes.as_mut_ptr(),
                nn,
                tinf,
                closed,
                1,
                rbs_masks,
            );
            record_overlapping_starts(buf.nodes.as_mut_ptr(), nn, tinf, 1);
            let ipath = dprog(buf.nodes.as_mut_ptr(), nn, tinf, 1);
            if ipath < 0 || ipath >= nn {
                continue;
            }

            let score = buf.nodes[ipath as usize].score;
            if score > best_score {
                best_score = score;
                eliminate_bad_genes(buf.nodes.as_mut_ptr(), ipath, tinf);

                let ng = add_genes(buf.genes.as_mut_ptr(), buf.nodes.as_mut_ptr(), ipath);
                tweak_final_starts(
                    buf.genes.as_mut_ptr(),
                    ng,
                    buf.nodes.as_mut_ptr(),
                    nn,
                    tinf,
                );
                let kept = (ng.max(0) as usize).min(buf.genes.len());
                best_genes.clear();
                best_genes.extend_from_slice(&buf.genes[..kept]);
                let dirty = (kept + 1).min(buf.genes.len());
                for gene in &mut buf.genes[..dirty] {
                    *gene = std::mem::zeroed();
                }

                best_nodes.clear();
                best_nodes.extend_from_slice(&buf.nodes[..nn as usize]);
                best_tinf = Some(i);
            }
        }
    }

    buf.mark_nodes(nn);

    let Some(best) = best_tinf else {
        return Ok(Vec::new());
    };
    let tinf: &mut Training = &mut models[best];

    unsafe {
        record_gene_data(
            best_genes.as_mut_ptr(),
            best_genes.len() as c_int,
            best_nodes.as_mut_ptr(),
            tinf,
            1,
        );

        let mut result = Vec::with_capacity(best_genes.len());
        for gene in &best_genes {
            result.push(gene_to_predicted(
                gene,
                best_nodes.as_ptr(),
                tinf,
                slen as usize,
            ));
        }
        Ok(result)
    }
}

//! Reusable metagenomic gene predictor.
//!
//! Caches the 50 metagenomic models and uses a rayon thread pool
//! with a large stack for prediction.

use std::alloc::{Layout, alloc_zeroed, handle_alloc_error};
use std::os::raw::c_int;
use std::sync::{Arc, Mutex};

use rayon::ThreadPool;
use rayon::prelude::*;

use super::convert::gene_to_predicted;
use super::encode::SequenceBuffer;
use super::types::{PredictedGene, ProdigalConfig, ProdigalError};
use crate::types::{Gene, MAX_SEQ, NUM_META, Training};

use super::predict::validate_config;

/// Recommended stack size for Rayon pools used by `MetaPredictor`.
///
/// The low-level translated prediction routines use deep call stacks, so the
/// default Rust thread stack can be too small for worker threads.
pub const META_PREDICTOR_STACK_SIZE: usize = 32 * 1024 * 1024; // 32 MB

use super::meta_score::{Scored, Shared, adopt, better, score_against, score_prepared};
use crate::gene::record_gene_data;
use crate::node::{add_nodes, calc_orf_gc, mer_table, record_rbs_masks, reset_node_scores};

/// Reusable metagenomic gene predictor.
///
/// Pre-loads 50 metagenomic models and evaluates qualifying models with
/// Prodigal-compatible state preservation across adjacent model bins.
pub struct MetaPredictor {
    pool: Arc<ThreadPool>,
    models: Arc<Vec<Box<Training>>>,
    config: ProdigalConfig,
    // A buffer holds about 48 MB that the allocator zeroes on every fresh one, and rayon would
    // build one per job, so buffers go back here between contigs and batches.
    spare: Mutex<Vec<SequenceBuffer>>,
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
            spare: Mutex::new(Vec::new()),
        })
    }

    fn predict_one(&self, seq: &[u8]) -> Result<Vec<PredictedGene>, ProdigalError> {
        let held = self.spare.lock().ok().and_then(|mut spare| spare.pop());
        let mut buf = held.unwrap_or_else(SequenceBuffer::reusable);
        let predicted = predict_parallel(seq, &self.models, &self.config, &mut buf);
        if let Ok(mut spare) = self.spare.lock() {
            spare.push(buf);
        }
        predicted
    }

    /// Predict genes in the given sequence.
    pub fn predict(&self, seq: &[u8]) -> Result<Vec<PredictedGene>, ProdigalError> {
        validate_sequence(seq)?;

        self.pool.install(|| self.predict_one(seq))
    }

    /// Predict genes in a batch of sequences, preserving input order.
    pub fn predict_batch<S>(&self, seqs: &[S]) -> Result<Vec<Vec<PredictedGene>>, ProdigalError>
    where
        S: AsRef<[u8]> + Sync,
    {
        for seq in seqs {
            validate_sequence(seq.as_ref())?;
        }

        self.pool.install(|| {
            seqs.par_iter()
                .map(|seq| self.predict_one(seq.as_ref()))
                .collect()
        })
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

struct NodeState {
    nn: c_int,
    masks_fresh: bool,
    built_table: Option<c_int>,
}

unsafe fn prepare_nodes(
    buf: &mut SequenceBuffer,
    tinf: &Training,
    slen: c_int,
    closed: c_int,
    state: &mut NodeState,
) -> *const u32 {
    if state.built_table != Some(tinf.trans_table) {
        state.built_table = Some(tinf.trans_table);
        buf.clear_nodes(state.nn);
        state.nn = add_nodes(
            buf.seq.as_mut_ptr(),
            buf.rseq.as_mut_ptr(),
            slen,
            buf.nodes.as_mut_ptr(),
            closed,
            buf.masks.as_mut_ptr(),
            buf.nmask,
            tinf,
        );
        buf.nodes[..state.nn as usize]
            .sort_unstable_by(|a, b| a.ndx.cmp(&b.ndx).then(b.strand.cmp(&a.strand)));
        calc_orf_gc(
            buf.seq.as_mut_ptr(),
            buf.rseq.as_mut_ptr(),
            slen,
            buf.nodes.as_mut_ptr(),
            state.nn,
            tinf,
        );
        state.masks_fresh = false;
    }

    if tinf.uses_sd == 1 && !state.masks_fresh {
        buf.ensure_rbs_capacity(state.nn);
        record_rbs_masks(
            buf.seq.as_mut_ptr(),
            buf.rseq.as_mut_ptr(),
            slen,
            buf.nodes.as_mut_ptr(),
            state.nn,
            buf.rbs_masks.as_mut_ptr(),
        );
        state.masks_fresh = true;
    }
    match state.masks_fresh {
        true => buf.rbs_masks.as_ptr(),
        false => std::ptr::null(),
    }
}

/// A long contig is one indivisible unit of work to `predict_batch`, so with one left the box
/// runs a single core. Splitting its models is the only parallelism available at that point.
const MODEL_SPLIT_BASES: c_int = 100_000;

fn predict_parallel(
    seq: &[u8],
    models: &[Box<Training>],
    config: &ProdigalConfig,
    buf: &mut SequenceBuffer,
) -> Result<Vec<PredictedGene>, ProdigalError> {
    let closed = if config.closed_ends { 1 } else { 0 };

    let (slen, gc) = unsafe { buf.encode(seq, config.mask_n_runs) };
    if slen == 0 {
        return Err(ProdigalError::EmptySequence);
    }
    buf.ensure_node_capacity(slen);
    unsafe {
        mer_table(buf.seq.as_mut_ptr(), slen, &mut buf.forward_mers);
        mer_table(buf.rseq.as_mut_ptr(), slen, &mut buf.reverse_mers);
    }

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
    let mut state = NodeState {
        nn: 0,
        masks_fresh: false,
        built_table: None,
    };

    let mut chosen: Vec<usize> = (0..NUM_META)
        .filter(|i| models[*i].gc >= low && models[*i].gc <= high)
        .collect();

    unsafe {
        if config.model_depth > 0 && chosen.len() > config.model_depth {
            let mut ranked: Vec<(usize, f64)> = Vec::with_capacity(chosen.len());
            for i in &chosen {
                let tinf: &Training = &models[*i];
                let rbs_masks = prepare_nodes(buf, tinf, slen, closed, &mut state);
                reset_node_scores(buf.nodes.as_mut_ptr(), state.nn);
                let shared = Shared::over(buf, rbs_masks);
                score_prepared(
                    &shared,
                    slen,
                    buf.nodes.as_mut_ptr(),
                    state.nn,
                    tinf,
                    closed,
                );
                let mut peak = f64::NEG_INFINITY;
                for node in &buf.nodes[..state.nn as usize] {
                    peak = peak.max(node.cscore + node.sscore);
                }
                ranked.push((*i, peak));
            }
            ranked.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
            ranked.truncate(config.model_depth);
            chosen = ranked.iter().map(|(i, _)| *i).collect();
            chosen.sort_unstable();
        }

        if slen < MODEL_SPLIT_BASES || chosen.len() < 2 {
            for i in &chosen {
                let tinf: &Training = &models[*i];
                let masks = prepare_nodes(buf, tinf, slen, closed, &mut state);
                let nn = state.nn;
                let shared = Shared::over(buf, masks);
                let mut nodes = std::mem::take(&mut buf.nodes);
                let (score, ipath) = score_against(&shared, slen, &mut nodes, nn, tinf, closed);
                buf.nodes = nodes;
                if ipath < 0 || ipath >= nn || score <= best_score {
                    continue;
                }
                best_score = score;
                best_tinf = Some(*i);
                let mut entry = Scored {
                    model: *i,
                    score,
                    ipath,
                    nodes: buf.nodes[..nn as usize].to_vec(),
                };
                adopt(
                    &mut entry,
                    tinf,
                    &mut buf.genes,
                    &mut best_genes,
                    &mut best_nodes,
                    nn,
                );
            }
        } else {
            let mut plan: Vec<(usize, usize, bool)> = Vec::with_capacity(chosen.len());
            let mut generation = 0usize;
            let mut table = state.built_table;
            let mut fresh = state.masks_fresh;
            for i in &chosen {
                if table != Some(models[*i].trans_table) {
                    table = Some(models[*i].trans_table);
                    fresh = false;
                    generation += 1;
                }
                fresh |= models[*i].uses_sd == 1;
                plan.push((*i, generation, fresh));
            }

            let mut at = 0;
            while at < plan.len() {
                let mut upto = at;
                while upto < plan.len() && plan[upto].1 == plan[at].1 {
                    upto += 1;
                }
                let group = &plan[at..upto];
                for (i, _, _) in group {
                    prepare_nodes(buf, &models[*i], slen, closed, &mut state);
                }
                let nn = state.nn;
                let masks = buf.rbs_masks.as_ptr();
                let with_masks = Shared::over(buf, masks);
                let bare = Shared::over(buf, std::ptr::null());
                let base = &buf.nodes[..nn as usize];
                let winner = group
                    .par_iter()
                    .fold(
                        || None::<Scored>,
                        |held, (i, _, masks)| {
                            let shared = match masks {
                                true => &with_masks,
                                false => &bare,
                            };
                            let mut nodes = base.to_vec();
                            let (score, ipath) =
                                score_against(shared, slen, &mut nodes, nn, &models[*i], closed);
                            match ipath >= 0 && ipath < nn {
                                false => held,
                                true => {
                                    let found = Scored {
                                        model: *i,
                                        score,
                                        ipath,
                                        nodes,
                                    };
                                    Some(match held {
                                        Some(held) => better(held, found),
                                        None => found,
                                    })
                                }
                            }
                        },
                    )
                    .reduce(
                        || None::<Scored>,
                        |a, b| match (a, b) {
                            (Some(a), Some(b)) => Some(better(a, b)),
                            (held, None) | (None, held) => held,
                        },
                    );
                if let Some(mut entry) = winner {
                    if entry.score > best_score {
                        best_score = entry.score;
                        best_tinf = Some(entry.model);
                        let tinf: &Training = &models[entry.model];
                        adopt(
                            &mut entry,
                            tinf,
                            &mut buf.genes,
                            &mut best_genes,
                            &mut best_nodes,
                            nn,
                        );
                    }
                }
                at = upto;
            }
        }
    }

    buf.mark_nodes(state.nn);

    let Some(best) = best_tinf else {
        return Ok(Vec::new());
    };
    let tinf: &Training = &models[best];

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

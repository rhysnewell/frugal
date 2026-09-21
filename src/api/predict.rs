use std::os::raw::c_int;

use super::convert::gene_to_predicted;
use super::encode::SequenceBuffer;
use super::training::TrainingData;
use super::types::{PredictedGene, ProdigalConfig, ProdigalError};
use crate::types::{Training, MAX_SEQ, NUM_META};

const MIN_SINGLE_GENOME: usize = 20_000;
const STACK_SIZE: usize = 32 * 1024 * 1024; // 32 MB

/// Run a closure on a scoped thread with a large stack to accommodate
/// the deep call stacks of the internal prediction functions.
fn with_large_stack<F, T>(f: F) -> T
where
    F: FnOnce() -> T + Send,
    T: Send,
{
    std::thread::scope(|scope| {
        std::thread::Builder::new()
            .stack_size(STACK_SIZE)
            .spawn_scoped(scope, f)
            .expect("failed to spawn worker thread")
            .join()
            .expect("worker thread panicked")
    })
}

use crate::dprog::{dprog, eliminate_bad_genes};
use crate::gene::{add_genes, record_gene_data, tweak_final_starts};
use crate::node::{
    add_nodes, calc_dicodon_gene, determine_sd_usage, raw_coding_score, rbs_score, record_gc_bias,
    record_overlapping_starts, reset_node_scores, score_nodes, train_starts_nonsd, train_starts_sd,
};
use crate::sequence::calc_most_gc_frame;

/// Validate that `config.translation_table` is one of the NCBI tables Prodigal supports.
pub(crate) fn validate_config(config: &ProdigalConfig) -> Result<(), ProdigalError> {
    let tt = config.translation_table;
    if tt < 1 || tt > 25 || tt == 7 || tt == 8 || (tt >= 17 && tt <= 20) {
        return Err(ProdigalError::InvalidTranslationTable(tt));
    }
    Ok(())
}

/// Predict genes in metagenomic mode with default settings.
pub fn predict_meta(seq: &[u8]) -> Result<Vec<PredictedGene>, ProdigalError> {
    predict_meta_with(seq, &ProdigalConfig::default())
}

/// Predict genes in metagenomic mode with custom settings.
pub fn predict_meta_with(
    seq: &[u8],
    config: &ProdigalConfig,
) -> Result<Vec<PredictedGene>, ProdigalError> {
    validate_config(config)?;
    if seq.is_empty() {
        return Err(ProdigalError::EmptySequence);
    }
    with_large_stack(|| predict_meta_inner(seq, config))
}

/// Single-threaded metagenomic prediction core.
///
/// Loads all 50 models, encodes the sequence once, then iterates models in
/// order, rebuilding nodes only when the translation table changes and
/// retaining the best-scoring DP solution to convert into `PredictedGene`s.
fn predict_meta_inner(
    seq: &[u8],
    config: &ProdigalConfig,
) -> Result<Vec<PredictedGene>, ProdigalError> {
    if seq.is_empty() {
        return Err(ProdigalError::EmptySequence);
    }
    if seq.len() > MAX_SEQ {
        return Err(ProdigalError::SequenceTooLong {
            length: seq.len(),
            max: MAX_SEQ,
        });
    }

    let mut buf = SequenceBuffer::with_base_capacity(seq.len());
    let closed = if config.closed_ends { 1 } else { 0 };

    // Initialize 50 metagenomic models
    let mut models: Vec<Box<Training>> = Vec::with_capacity(NUM_META);
    for i in 0..NUM_META {
        let mut tinf: Box<Training> = Box::new(unsafe { std::mem::zeroed() });
        unsafe {
            crate::training_data::load_metagenome(i, &mut *tinf as *mut Training);
        }
        models.push(tinf);
    }

    // Encode the input sequence
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

    let mut max_score: f64 = -100.0;
    let mut best_nodes = Vec::new();
    let mut best_genes: Vec<crate::types::Gene> = Vec::new();
    let mut best_tinf: Option<Training> = None;
    let mut nn: c_int = 0;
    unsafe {
        for i in 0..NUM_META {
            let need_rebuild = i == 0 || models[i].trans_table != models[i - 1].trans_table;
            let tinf = &mut *models[i];

            if need_rebuild {
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
            }

            if tinf.gc < low || tinf.gc > high {
                continue;
            }

            reset_node_scores(buf.nodes.as_mut_ptr(), nn);
            score_nodes(
                buf.seq.as_mut_ptr(),
                buf.rseq.as_mut_ptr(),
                slen,
                buf.nodes.as_mut_ptr(),
                nn,
                tinf,
                closed,
                1,
            );
            record_overlapping_starts(buf.nodes.as_mut_ptr(), nn, tinf, 1);
            let ipath = dprog(buf.nodes.as_mut_ptr(), nn, tinf, 1);
            if ipath < 0 || ipath >= nn {
                continue;
            }

            if buf.nodes[ipath as usize].score > max_score {
                max_score = buf.nodes[ipath as usize].score;
                eliminate_bad_genes(buf.nodes.as_mut_ptr(), ipath, tinf);
                let ng = add_genes(buf.genes.as_mut_ptr(), buf.nodes.as_mut_ptr(), ipath);
                tweak_final_starts(buf.genes.as_mut_ptr(), ng, buf.nodes.as_mut_ptr(), nn, tinf);

                best_nodes = buf.nodes[..nn as usize].to_vec();
                best_genes = buf.genes[..ng as usize].to_vec();
                best_tinf = Some((*tinf).clone());
            }
        }

        let Some(mut tinf) = best_tinf else {
            return Ok(Vec::new());
        };
        record_gene_data(
            best_genes.as_mut_ptr(),
            best_genes.len() as c_int,
            best_nodes.as_mut_ptr(),
            &mut tinf,
            1,
        );

        let mut result = Vec::with_capacity(best_genes.len());
        for gene in &best_genes {
            result.push(gene_to_predicted(
                gene,
                best_nodes.as_ptr(),
                &tinf,
                slen as usize,
            ));
        }
        Ok(result)
    }
}

/// Train on a genome sequence (single-genome mode).
///
/// The sequence should be >= 20,000 bp. Multiple contigs can be provided
/// as separate entries; they will be concatenated with stop-codon spacers.
pub fn train(seq: &[u8]) -> Result<TrainingData, ProdigalError> {
    train_with(seq, &ProdigalConfig::default())
}

/// Train with custom settings.
pub fn train_with(seq: &[u8], config: &ProdigalConfig) -> Result<TrainingData, ProdigalError> {
    validate_config(config)?;
    if seq.is_empty() {
        return Err(ProdigalError::EmptySequence);
    }
    with_large_stack(|| train_inner(seq, config))
}

/// Core single-genome training routine.
///
/// Encodes the input, finds start/stop nodes, computes GC frame bias, runs an
/// initial dynamic-programming pass, trains dicodon coding statistics, then
/// trains RBS/start-codon weights (SD or non-SD) and returns the learned
/// `TrainingData`.
fn train_inner(seq: &[u8], config: &ProdigalConfig) -> Result<TrainingData, ProdigalError> {
    let mut buf = SequenceBuffer::with_base_capacity(seq.len() + 24);
    let mut tinf = Box::new(unsafe { std::mem::zeroed::<Training>() });
    tinf.st_wt = 4.35;
    tinf.trans_table = config.translation_table as c_int;
    let closed = if config.closed_ends { 1 } else { 0 };

    // Encode sequence(s) for training
    let (slen, gc) = unsafe { buf.encode_training(&[seq], config.mask_n_runs) };
    tinf.gc = gc;

    if (slen as usize) < MIN_SINGLE_GENOME {
        return Err(ProdigalError::SequenceTooShort {
            length: slen as usize,
            min: MIN_SINGLE_GENOME,
        });
    }

    buf.ensure_node_capacity(slen);

    unsafe {
        // Find all potential starts and stops
        let nn = add_nodes(
            buf.seq.as_mut_ptr(),
            buf.rseq.as_mut_ptr(),
            slen,
            buf.nodes.as_mut_ptr(),
            closed,
            buf.masks.as_mut_ptr(),
            buf.nmask,
            &mut *tinf,
        );
        buf.nodes[..nn as usize]
            .sort_unstable_by(|a, b| a.ndx.cmp(&b.ndx).then(b.strand.cmp(&a.strand)));

        // GC frame bias
        let gc_frame = calc_most_gc_frame(buf.seq.as_mut_ptr(), slen);
        if gc_frame.is_null() {
            return Err(ProdigalError::EmptySequence);
        }
        record_gc_bias(gc_frame, buf.nodes.as_mut_ptr(), nn, &mut *tinf);
        drop(Vec::from_raw_parts(gc_frame, slen as usize, slen as usize));

        // Initial DP with GC frame bias
        record_overlapping_starts(buf.nodes.as_mut_ptr(), nn, &mut *tinf, 0);
        let ipath = dprog(buf.nodes.as_mut_ptr(), nn, &mut *tinf, 0);

        // Dicodon statistics from initial gene set
        calc_dicodon_gene(
            &mut *tinf,
            buf.seq.as_mut_ptr(),
            buf.rseq.as_mut_ptr(),
            slen,
            buf.nodes.as_mut_ptr(),
            ipath,
        );
        raw_coding_score(
            buf.seq.as_mut_ptr(),
            buf.rseq.as_mut_ptr(),
            slen,
            buf.nodes.as_mut_ptr(),
            nn,
            &mut *tinf,
        );

        // RBS and start training
        rbs_score(
            buf.seq.as_mut_ptr(),
            buf.rseq.as_mut_ptr(),
            slen,
            buf.nodes.as_mut_ptr(),
            nn,
            &mut *tinf,
        );
        train_starts_sd(
            buf.seq.as_mut_ptr(),
            buf.rseq.as_mut_ptr(),
            slen,
            buf.nodes.as_mut_ptr(),
            nn,
            &mut *tinf,
        );
        determine_sd_usage(&mut *tinf);

        if config.force_non_sd {
            tinf.uses_sd = 0;
        }
        if tinf.uses_sd == 0 {
            train_starts_nonsd(
                buf.seq.as_mut_ptr(),
                buf.rseq.as_mut_ptr(),
                slen,
                buf.nodes.as_mut_ptr(),
                nn,
                &mut *tinf,
            );
        }
    }

    Ok(TrainingData { inner: tinf })
}

/// Predict genes using pre-trained model (single-genome mode).
pub fn predict(seq: &[u8], training: &TrainingData) -> Result<Vec<PredictedGene>, ProdigalError> {
    predict_with(seq, training, &ProdigalConfig::default())
}

/// Predict with custom settings.
pub fn predict_with(
    seq: &[u8],
    training: &TrainingData,
    config: &ProdigalConfig,
) -> Result<Vec<PredictedGene>, ProdigalError> {
    validate_config(config)?;
    if seq.is_empty() {
        return Err(ProdigalError::EmptySequence);
    }
    if seq.len() > MAX_SEQ {
        return Err(ProdigalError::SequenceTooLong {
            length: seq.len(),
            max: MAX_SEQ,
        });
    }
    with_large_stack(|| predict_inner(seq, training, config))
}

/// Core single-genome prediction routine.
///
/// Encodes the sequence, builds and scores nodes against `training`, runs the
/// dynamic-programming gene-path search, finalizes gene starts, and converts
/// the resulting genes into `PredictedGene`s.
fn predict_inner(
    seq: &[u8],
    training: &TrainingData,
    config: &ProdigalConfig,
) -> Result<Vec<PredictedGene>, ProdigalError> {
    let mut buf = SequenceBuffer::with_base_capacity(seq.len());
    let closed = if config.closed_ends { 1 } else { 0 };

    // We need a mutable copy of training for the internal functions (keep on heap — 558KB)
    let mut tinf_box = training.inner.clone();
    let tinf: *mut Training = &mut *tinf_box;

    let (slen, _gc) = unsafe { buf.encode(seq, config.mask_n_runs) };
    if slen == 0 {
        return Err(ProdigalError::EmptySequence);
    }
    buf.ensure_node_capacity(slen);

    unsafe {
        let nn = add_nodes(
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

        score_nodes(
            buf.seq.as_mut_ptr(),
            buf.rseq.as_mut_ptr(),
            slen,
            buf.nodes.as_mut_ptr(),
            nn,
            tinf,
            closed,
            0,
        );
        record_overlapping_starts(buf.nodes.as_mut_ptr(), nn, tinf, 1);
        let ipath = dprog(buf.nodes.as_mut_ptr(), nn, tinf, 1);
        eliminate_bad_genes(buf.nodes.as_mut_ptr(), ipath, tinf);
        let ng = add_genes(buf.genes.as_mut_ptr(), buf.nodes.as_mut_ptr(), ipath);
        tweak_final_starts(buf.genes.as_mut_ptr(), ng, buf.nodes.as_mut_ptr(), nn, tinf);
        record_gene_data(buf.genes.as_mut_ptr(), ng, buf.nodes.as_mut_ptr(), tinf, 1);

        let mut result = Vec::with_capacity(ng as usize);
        for i in 0..ng {
            result.push(gene_to_predicted(
                &buf.genes[i as usize],
                buf.nodes.as_ptr(),
                &*tinf,
                slen as usize,
            ));
        }
        Ok(result)
    }
}

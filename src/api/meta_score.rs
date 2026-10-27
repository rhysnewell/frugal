use std::os::raw::c_int;

use super::encode::SequenceBuffer;
use crate::dprog::{dprog, eliminate_bad_genes};
use crate::gene::{add_genes, tweak_final_starts};
use crate::node::{
    mer_from, raw_coding_score_with, record_overlapping_starts, reset_node_scores, score_starts,
};
use crate::types::{Gene, Node, Training};

pub(super) struct Scored {
    pub(super) model: usize,
    pub(super) score: f64,
    pub(super) ipath: c_int,
    pub(super) nodes: Vec<Node>,
}

/// Read-only for the length of a scoring pass, and raw because it crosses into the C port.
pub(super) struct Shared {
    seq: *mut u8,
    rseq: *mut u8,
    masks: *const u32,
    forward: *const u16,
    reverse: *const u16,
}

impl Shared {
    pub(super) fn over(buf: &mut SequenceBuffer, masks: *const u32) -> Self {
        Shared {
            seq: buf.seq.as_mut_ptr(),
            rseq: buf.rseq.as_mut_ptr(),
            masks,
            forward: buf.forward_mers.as_ptr(),
            reverse: buf.reverse_mers.as_ptr(),
        }
    }
}

unsafe impl Send for Shared {}
unsafe impl Sync for Shared {}

/// The serial loop keeps the first model to reach a score, so a tie falls to the lower index.
pub(super) fn better(a: Scored, b: Scored) -> Scored {
    match b.score > a.score || (b.score == a.score && b.model < a.model) {
        true => b,
        false => a,
    }
}

pub(super) unsafe fn score_against(
    shared: &Shared,
    slen: c_int,
    nodes: &mut [Node],
    nn: c_int,
    tinf: &Training,
    closed: c_int,
) -> (f64, c_int) {
    reset_node_scores(nodes.as_mut_ptr(), nn);
    score_prepared(shared, slen, nodes.as_mut_ptr(), nn, tinf, closed);
    record_overlapping_starts(nodes.as_mut_ptr(), nn, tinf, 1);
    let ipath = dprog(nodes.as_mut_ptr(), nn, tinf, 1);
    let score = match ipath >= 0 && ipath < nn {
        true => nodes[ipath as usize].score,
        false => f64::NEG_INFINITY,
    };
    (score, ipath)
}

pub(super) unsafe fn score_prepared(
    shared: &Shared,
    slen: c_int,
    nod: *mut Node,
    nn: c_int,
    tinf: &Training,
    closed: c_int,
) {
    let forward = std::slice::from_raw_parts(shared.forward, slen as usize);
    let reverse = std::slice::from_raw_parts(shared.reverse, slen as usize);
    raw_coding_score_with(
        nod,
        nn,
        tinf,
        |at| usize::from(forward[at as usize]),
        |at| usize::from(reverse[(slen - at - 1) as usize]),
    );
    let mer = |strand_forward: bool, len: c_int, at: c_int| match strand_forward {
        true => mer_from(forward, len, at),
        false => mer_from(reverse, len, at),
    };
    score_starts(
        shared.seq,
        shared.rseq,
        slen,
        nod,
        nn,
        tinf,
        closed,
        1,
        shared.masks,
        mer,
    );
}

pub(super) unsafe fn adopt(
    entry: &mut Scored,
    tinf: &Training,
    genes: &mut [Gene],
    best_genes: &mut Vec<Gene>,
    best_nodes: &mut Vec<Node>,
    nn: c_int,
) {
    eliminate_bad_genes(entry.nodes.as_mut_ptr(), entry.ipath, tinf);
    let ng = add_genes(genes.as_mut_ptr(), entry.nodes.as_mut_ptr(), entry.ipath);
    tweak_final_starts(genes.as_mut_ptr(), ng, entry.nodes.as_mut_ptr(), nn, tinf);
    let kept = (ng.max(0) as usize).min(genes.len());
    best_genes.clear();
    best_genes.extend_from_slice(&genes[..kept]);
    let dirty = (kept + 1).min(genes.len());
    for gene in &mut genes[..dirty] {
        *gene = std::mem::zeroed();
    }
    best_nodes.clear();
    best_nodes.append(&mut entry.nodes);
}

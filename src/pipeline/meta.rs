use std::os::raw::c_int;

use super::sorted_nodes;
use crate::dprog::{dprog, eliminate_bad_genes};
use crate::gene::{add_genes, record_gene_data, tweak_final_starts};
use crate::node::{record_overlapping_starts, reset_node_scores, score_nodes};
use crate::training_data::load_metagenome;
use crate::types::{Gene, Mask, NUM_META, Node, Training};

// The winning bin keeps its genes, but a later bin with another translation table rebuilt the
// nodes, so they are built again for the winner before anything is written.
pub(super) unsafe fn best_bin(
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    gc: f64,
    nodes: *mut Node,
    nn: &mut c_int,
    genes: *mut Gene,
    ng: &mut c_int,
    closed: c_int,
    mlist: *mut Mask,
    nmask: c_int,
    num_seq: c_int,
    meta_gc: &[f64; NUM_META],
    meta_trans_table: &[c_int; NUM_META],
    meta_tinf: *mut Training,
    best_tinf: *mut Training,
    max_phase: &mut c_int,
) {
    let mut low = 0.88495 * gc - 0.0102337;
    if low > 0.65 {
        low = 0.65;
    }
    let mut high = 0.86596 * gc + 0.1131991;
    if high < 0.35 {
        high = 0.35;
    }

    let mut max_score = -100.0;
    *ng = 0;
    for mi in 0..NUM_META {
        if mi == 0 || meta_trans_table[mi] != meta_trans_table[mi - 1] {
            (*meta_tinf).trans_table = meta_trans_table[mi];
            std::ptr::write_bytes(nodes, 0, *nn as usize);
            *nn = sorted_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, meta_tinf);
        }
        if meta_gc[mi] < low || meta_gc[mi] > high {
            continue;
        }
        load_metagenome(mi, meta_tinf);
        reset_node_scores(nodes, *nn);
        score_nodes(seq, rseq, slen, nodes, *nn, meta_tinf, closed, 1);
        record_overlapping_starts(nodes, *nn, meta_tinf, 1);
        let ipath = dprog(nodes, *nn, meta_tinf, 1);
        if ipath >= 0 && (*nodes.offset(ipath as isize)).score > max_score {
            *max_phase = mi as c_int;
            max_score = (*nodes.offset(ipath as isize)).score;
            eliminate_bad_genes(nodes, ipath, meta_tinf);
            *ng = add_genes(genes, nodes, ipath);
            tweak_final_starts(genes, *ng, nodes, *nn, meta_tinf);
            record_gene_data(genes, *ng, nodes, meta_tinf, num_seq);
        }
    }

    load_metagenome(*max_phase as usize, best_tinf);
    std::ptr::write_bytes(nodes, 0, *nn as usize);
    *nn = sorted_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, best_tinf);
    score_nodes(seq, rseq, slen, nodes, *nn, best_tinf, closed, 1);
}

use std::mem::MaybeUninit;
use std::os::raw::{c_int, c_void};

use super::sorted_nodes;
use crate::dprog::dprog;
use crate::node::{
    calc_dicodon_gene, determine_sd_usage, node_room, raw_coding_score, rbs_score, record_gc_bias,
    record_overlapping_starts, train_starts_nonsd, train_starts_sd,
};
use crate::sequence::{calc_most_gc_frame, rcom_seq, read_seq_training};
use crate::types::{Mask, Node, Training};

const MIN_SINGLE_GENOME: c_int = 20000;
const IDEAL_SINGLE_GENOME: c_int = 100000;

pub(super) unsafe fn train_single(
    input: *mut c_void,
    seq: *mut u8,
    rseq: *mut u8,
    useq: *mut u8,
    nodes_vec: &mut Vec<MaybeUninit<Node>>,
    tinf: *mut Training,
    do_mask: c_int,
    mlist: *mut Mask,
    nmask: &mut c_int,
    closed: c_int,
    force_nonsd: c_int,
    quiet: c_int,
) -> Result<(c_int, *mut Node, c_int), i32> {
    if quiet == 0 {
        eprintln!("Request:  Single Genome, Phase:  Training");
        eprint!("Reading in the sequence(s) to train...");
    }
    let slen = read_seq_training(input, seq, useq, &mut (*tinf).gc, do_mask, mlist, nmask);
    if slen == 0 {
        eprint!("\n\nSequence read failed (file must be Fasta, ");
        eprintln!("Genbank, or EMBL format).\n");
        return Err(9);
    }
    if slen < MIN_SINGLE_GENOME {
        eprint!("\n\nError:  Sequence must be {}", MIN_SINGLE_GENOME);
        eprint!(" characters (only {} read).\n(Consider", slen);
        eprint!(" running with the -p meta option or finding");
        eprintln!(" more contigs from the same genome.)\n");
        return Err(10);
    }
    if slen < IDEAL_SINGLE_GENOME {
        eprint!("\n\nWarning:  ideally Prodigal should be given at");
        eprint!(" least {} bases for ", IDEAL_SINGLE_GENOME);
        eprint!("training.\nYou may get better results with the ");
        eprintln!("-p meta option.\n");
    }
    rcom_seq(seq, rseq, useq, slen);
    if quiet == 0 {
        eprintln!("{} bp seq created, {:.2} pct GC", slen, (*tinf).gc * 100.0);
    }

    if quiet == 0 {
        eprint!("Locating all potential starts and stops...");
    }
    let nodes = node_room(nodes_vec, slen);
    let nn = sorted_nodes(seq, rseq, slen, nodes, closed, mlist, *nmask, tinf);
    if quiet == 0 {
        eprintln!("{} nodes", nn);
    }

    if quiet == 0 {
        eprint!("Looking for GC bias in different frames...");
    }
    let gc_frame = calc_most_gc_frame(seq, slen);
    if gc_frame.is_null() {
        eprintln!("Malloc failed on gc frame plot\n");
        return Err(11);
    }
    record_gc_bias(gc_frame, nodes, nn, tinf);
    if quiet == 0 {
        eprintln!(
            "frame bias scores: {:.2} {:.2} {:.2}",
            (*tinf).bias[0],
            (*tinf).bias[1],
            (*tinf).bias[2]
        );
    }
    drop(Vec::from_raw_parts(gc_frame, slen as usize, slen as usize));

    if quiet == 0 {
        eprint!("Building initial set of genes to train from...");
    }
    record_overlapping_starts(nodes, nn, tinf, 0);
    let ipath = dprog(nodes, nn, tinf, 0);
    if quiet == 0 {
        eprintln!("done!");
    }

    if quiet == 0 {
        eprint!("Creating coding model and scoring nodes...");
    }
    calc_dicodon_gene(tinf, seq, rseq, slen, nodes, ipath);
    raw_coding_score(seq, rseq, slen, nodes, nn, tinf);
    if quiet == 0 {
        eprintln!("done!");
    }

    if quiet == 0 {
        eprint!("Examining upstream regions and training starts...");
    }
    rbs_score(seq, rseq, slen, nodes, nn, tinf);
    train_starts_sd(seq, rseq, slen, nodes, nn, tinf);
    determine_sd_usage(tinf);
    if force_nonsd == 1 {
        (*tinf).uses_sd = 0;
    }
    if (*tinf).uses_sd == 0 {
        train_starts_nonsd(seq, rseq, slen, nodes, nn, tinf);
    }
    if quiet == 0 {
        eprintln!("done!");
    }
    Ok((slen, nodes, nn))
}

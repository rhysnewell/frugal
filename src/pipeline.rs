/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    Rust port of main.c — the pipeline orchestration module.
*******************************************************************************/

use std::ffi::CStr;
use std::mem::MaybeUninit;
use std::os::raw::{c_char, c_int};

use crate::types::{Gene, Mask, MetagenomicBin, Node, Training};
use crate::types::{MAX_GENES, MAX_LINE, MAX_MASKS, MAX_SEQ, NUM_META, STT_NOD};

mod emit;
mod handles;
mod meta;
mod train;

const VERSION: &str = "2.6.3";
const DATE: &str = "February, 2016";
const VERSION_CSTR: &[u8] = b"2.6.3\0";

/// Parsed pipeline configuration — constructed by the CLI (clap) layer.
pub struct PipelineConfig {
    pub input_file: Option<String>,
    pub output_file: Option<String>,
    pub trans_file: Option<String>,
    pub nuc_file: Option<String>,
    pub start_file: Option<String>,
    pub train_file: Option<String>,
    pub output_format: i32,
    pub trans_table: i32,
    pub closed: bool,
    pub do_mask: bool,
    pub force_nonsd: bool,
    pub is_meta: bool,
    pub quiet: bool,
}

use crate::dprog::{dprog, eliminate_bad_genes};
use crate::gene::{add_genes, record_gene_data, tweak_final_starts};
use crate::metagenomic::initialize_metagenomic_metadata;
use crate::node::{add_nodes, node_room, record_overlapping_starts, score_nodes, write_start_file};
use crate::reader::seq_reader_seek;
use crate::sequence::{calc_short_header, next_seq_multi, rcom_seq};
use crate::training::{read_training_file, write_training_file};

use emit::write_genes;
use handles::Handles;

unsafe fn sorted_nodes(
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nodes: *mut Node,
    closed: c_int,
    mlist: *mut Mask,
    nmask: c_int,
    tinf: *const Training,
) -> c_int {
    let nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, tinf);
    std::slice::from_raw_parts_mut(nodes, nn as usize)
        .sort_unstable_by(|a, b| a.ndx.cmp(&b.ndx).then(b.strand.cmp(&a.strand)));
    nn
}

// ---------------------------------------------------------------------------
// Main pipeline: replaces C main()
// ---------------------------------------------------------------------------
/// Orchestrate the full Prodigal pipeline (port of C `main`).
///
/// Opens input/output handles, optionally reads or writes a training file, and
/// then for each sequence either runs the single-genome path (self-training
/// followed by gene finding) or the metagenomic path (scoring against the 50
/// pre-built bins and keeping the best-scoring one). Returns the process exit
/// code: 0 on success, non-zero on failure (matching the C exit codes).
#[allow(unused_assignments)]
pub unsafe fn run_pipeline(config: &PipelineConfig) -> i32 {
    let train_cstr = config
        .train_file
        .as_ref()
        .map(|s| std::ffi::CString::new(s.as_str()).unwrap());

    // Variable declarations
    let mut rv: c_int;
    let mut slen: c_int;
    let mut nn: c_int;
    let mut ng: c_int;
    let mut ipath: c_int;
    let mut do_training: c_int;
    let output: c_int = config.output_format;
    let mut max_phase: c_int;
    let closed: c_int = config.closed as c_int;
    let do_mask: c_int = config.do_mask as c_int;
    let mut nmask: c_int;
    let force_nonsd: c_int = config.force_nonsd as c_int;
    let user_tt: c_int = config.trans_table;
    let is_meta: c_int = config.is_meta as c_int;
    let mut num_seq: c_int;
    let quiet: c_int = config.quiet as c_int;
    let mut gc: f64 = 0.0;

    let train_file: *const c_char = train_cstr.as_ref().map_or(std::ptr::null(), |c| c.as_ptr());

    let mut cur_header: [c_char; MAX_LINE] = [0; MAX_LINE];
    let mut new_header: [c_char; MAX_LINE] = [0; MAX_LINE];
    let mut short_header: [c_char; MAX_LINE] = [0; MAX_LINE];

    // Reserve the same C-style work buffers without touching every page up
    // front. Nodes and genes are zeroed as individual records are emitted.
    let mut seq_vec: Vec<u8> = vec![0u8; MAX_SEQ / 4];
    let mut rseq_vec: Vec<u8> = vec![0u8; MAX_SEQ / 4];
    let mut useq_vec: Vec<u8> = vec![0u8; MAX_SEQ / 8];
    let mut nodes_vec: Vec<MaybeUninit<Node>> = Vec::with_capacity(STT_NOD);
    let mut genes_vec: Vec<MaybeUninit<Gene>> = Vec::with_capacity(MAX_GENES);

    let seq = seq_vec.as_mut_ptr();
    let rseq = rseq_vec.as_mut_ptr();
    let useq = useq_vec.as_mut_ptr();
    let mut nodes = nodes_vec.as_mut_ptr() as *mut Node;
    let genes = genes_vec.as_mut_ptr() as *mut Gene;

    let mut tinf: Training = std::mem::zeroed();

    let mut meta: [MetagenomicBin; NUM_META] = std::mem::zeroed();
    let mut mlist: [Mask; MAX_MASKS] = [Mask { begin: 0, end: 0 }; MAX_MASKS];

    let mut meta_gc: [f64; NUM_META] = [0.0; NUM_META];
    let mut meta_trans_table: [c_int; NUM_META] = [0; NUM_META];
    let mut meta_tinf: Training = std::mem::zeroed();
    let mut best_tinf: Training = std::mem::zeroed();

    // Initialize variables
    nn = 0;
    slen = 0;
    ipath = 0;
    ng = 0;
    nmask = 0;
    num_seq = 0;
    max_phase = 0;
    do_training = 0;

    // Set default training parameters
    tinf.st_wt = 4.35;
    tinf.trans_table = if user_tt > 0 { user_tt } else { 11 };

    // Print header
    if quiet == 0 {
        eprintln!("-------------------------------------");
        eprintln!("PRODIGAL v{} [{}]         ", VERSION, DATE);
        eprintln!("Univ of Tenn / Oak Ridge National Lab");
        eprintln!("Doug Hyatt, Loren Hauser, et al.     ");
        eprintln!("-------------------------------------");
    }

    // Read in the training file (if specified)
    if !train_file.is_null() {
        if is_meta == 1 {
            eprint!("\nError: cannot specify metagenomic sequence with a");
            eprintln!(" training file.");
            return 2;
        }
        rv = read_training_file(train_file, &mut tinf);
        if rv == 1 {
            do_training = 1;
        } else {
            if force_nonsd == 1 {
                eprint!("\nError: cannot force non-SD finder with a training");
                eprintln!(" file already created!");
                return 3;
            }
            if quiet == 0 {
                let tf_str = CStr::from_ptr(train_file).to_string_lossy();
                eprint!("Reading in training data from file {}...", tf_str);
            }
            if user_tt > 0 && user_tt != tinf.trans_table {
                eprint!("\n\nWarning: user-specified translation table does");
                eprintln!("not match the one in the specified training file! \n");
            }
            if rv == -1 {
                eprintln!("\n\nError: training file did not read correctly!");
                return 4;
            }
            if quiet == 0 {
                eprintln!("done!");
                eprintln!("-------------------------------------");
            }
        }
    }

    let handles = match Handles::open(config) {
        Ok(handles) => handles,
        Err(code) => return code,
    };

    // =========================================================================
    // Single Genome Training
    // =========================================================================
    if is_meta == 0 && (do_training == 1 || (do_training == 0 && train_file.is_null())) {
        match train::train_single(
            handles.input,
            seq,
            rseq,
            useq,
            &mut nodes_vec,
            &mut tinf,
            do_mask,
            mlist.as_mut_ptr(),
            &mut nmask,
            closed,
            force_nonsd,
            quiet,
        ) {
            Ok((trained, trained_nodes, trained_nn)) => {
                slen = trained;
                nodes = trained_nodes;
                nn = trained_nn;
            }
            Err(code) => return code,
        }

        // If training specified, write the training file and exit
        if do_training == 1 {
            if quiet == 0 {
                let tf_str = CStr::from_ptr(train_file).to_string_lossy();
                eprint!("Writing data to training file {}...", tf_str);
            }
            rv = write_training_file(train_file, &mut tinf);
            if rv != 0 {
                eprintln!("\nError: could not write training file!");
                return 12;
            } else {
                if quiet == 0 {
                    eprintln!("done!");
                }
                handles.close();
                return 0;
            }
        }

        // Rewind input file
        if quiet == 0 {
            eprintln!("-------------------------------------");
        }
        if seq_reader_seek(handles.input, 0, 0) == -1 {
            // SEEK_SET = 0
            eprintln!("\nError: could not rewind input file.");
            return 13;
        }

        // Reset sequence/DP variables
        std::ptr::write_bytes(seq, 0, (slen as usize / 4 + 1).min(seq_vec.len()));
        std::ptr::write_bytes(rseq, 0, (slen as usize / 4 + 1).min(rseq_vec.len()));
        std::ptr::write_bytes(useq, 0, (slen as usize / 8 + 1).min(useq_vec.len()));
        std::ptr::write_bytes(nodes, 0, nn as usize);
        nn = 0;
        slen = 0;
        ipath = 0;
        nmask = 0;
    }
    // Initialize metagenomic bins
    else if is_meta == 1 {
        if quiet == 0 {
            eprintln!("Request:  Metagenomic, Phase:  Training");
            eprint!("Initializing training files...");
        }
        initialize_metagenomic_metadata(
            meta.as_mut_ptr(),
            meta_gc.as_mut_ptr(),
            meta_trans_table.as_mut_ptr(),
        );
        if quiet == 0 {
            eprintln!("done!");
            eprintln!("-------------------------------------");
        }
    }

    // Print header for gene finding phase
    if quiet == 0 {
        if is_meta == 1 {
            eprintln!("Request:  Metagenomic, Phase:  Gene Finding");
        } else {
            eprintln!("Request:  Single Genome, Phase:  Gene Finding");
        }
    }

    // Read and process each sequence
    {
        let s = std::ffi::CString::new("Prodigal_Seq_1").unwrap();
        let bytes = s.as_bytes_with_nul();
        std::ptr::copy_nonoverlapping(
            bytes.as_ptr() as *const c_char,
            cur_header.as_mut_ptr(),
            bytes.len(),
        );
    }
    {
        let s = std::ffi::CString::new("Prodigal_Seq_2").unwrap();
        let bytes = s.as_bytes_with_nul();
        std::ptr::copy_nonoverlapping(
            bytes.as_ptr() as *const c_char,
            new_header.as_mut_ptr(),
            bytes.len(),
        );
    }

    loop {
        slen = next_seq_multi(
            handles.input,
            seq,
            useq,
            &mut num_seq,
            &mut gc,
            do_mask,
            mlist.as_mut_ptr(),
            &mut nmask,
            cur_header.as_mut_ptr(),
            new_header.as_mut_ptr(),
        );
        if slen == -1 {
            break;
        }

        rcom_seq(seq, rseq, useq, slen);
        if slen == 0 {
            eprint!("\nSequence read failed (file must be Fasta, ");
            eprintln!("Genbank, or EMBL format).\n");
            return 14;
        }

        if quiet == 0 {
            eprint!("Finding genes in sequence #{} ({} bp)...", num_seq, slen);
        }

        nodes = node_room(&mut nodes_vec, slen);

        // Calculate short header
        calc_short_header(cur_header.as_mut_ptr(), short_header.as_mut_ptr(), num_seq);

        if is_meta == 0 {
            // ---- Single Genome Version ----
            nn = sorted_nodes(
                seq,
                rseq,
                slen,
                nodes,
                closed,
                mlist.as_mut_ptr(),
                nmask,
                &tinf,
            );
            score_nodes(seq, rseq, slen, nodes, nn, &tinf, closed, is_meta);
            if handles.start != handles.stdout {
                write_start_file(
                    handles.start,
                    nodes,
                    nn,
                    &tinf,
                    num_seq,
                    slen,
                    0,
                    std::ptr::null_mut(),
                    VERSION_CSTR.as_ptr() as *mut c_char,
                    cur_header.as_mut_ptr(),
                );
            }
            record_overlapping_starts(nodes, nn, &tinf, 1);
            ipath = dprog(nodes, nn, &tinf, 1);
            eliminate_bad_genes(nodes, ipath, &tinf);
            ng = add_genes(genes, nodes, ipath);
            tweak_final_starts(genes, ng, nodes, nn, &tinf);
            record_gene_data(genes, ng, nodes, &tinf, num_seq);
            if quiet == 0 {
                eprintln!("done!");
            }
            write_genes(
                &handles,
                genes,
                ng,
                nodes,
                seq,
                rseq,
                useq,
                slen,
                output,
                num_seq,
                std::ptr::null_mut(),
                &tinf,
                cur_header.as_mut_ptr(),
                short_header.as_mut_ptr(),
            );
        } else {
            // ---- Metagenomic Version ----
            meta::best_bin(
                seq,
                rseq,
                slen,
                gc,
                nodes,
                &mut nn,
                genes,
                &mut ng,
                closed,
                mlist.as_mut_ptr(),
                nmask,
                num_seq,
                &meta_gc,
                &meta_trans_table,
                &mut meta_tinf,
                &mut best_tinf,
                &mut max_phase,
            );
            if handles.start != handles.stdout {
                write_start_file(
                    handles.start,
                    nodes,
                    nn,
                    &best_tinf,
                    num_seq,
                    slen,
                    1,
                    meta[max_phase as usize].desc.as_mut_ptr(),
                    VERSION_CSTR.as_ptr() as *mut c_char,
                    cur_header.as_mut_ptr(),
                );
            }
            if quiet == 0 {
                eprintln!("done!");
            }
            write_genes(
                &handles,
                genes,
                ng,
                nodes,
                seq,
                rseq,
                useq,
                slen,
                output,
                num_seq,
                meta[max_phase as usize].desc.as_mut_ptr(),
                &best_tinf,
                cur_header.as_mut_ptr(),
                short_header.as_mut_ptr(),
            );
        }

        // Reset sequence/DP variables
        std::ptr::write_bytes(seq, 0, (slen as usize / 4 + 1).min(seq_vec.len()));
        std::ptr::write_bytes(rseq, 0, (slen as usize / 4 + 1).min(rseq_vec.len()));
        std::ptr::write_bytes(useq, 0, (slen as usize / 8 + 1).min(useq_vec.len()));
        std::ptr::write_bytes(nodes, 0, nn as usize);
        nn = 0;
        slen = 0;
        ipath = 0;
        nmask = 0;
        std::ptr::copy_nonoverlapping(new_header.as_ptr(), cur_header.as_mut_ptr(), MAX_LINE);
        {
            let s = format!("Prodigal_Seq_{}\n", num_seq + 1);
            let c = std::ffi::CString::new(s).unwrap();
            let bytes = c.as_bytes_with_nul();
            std::ptr::copy_nonoverlapping(
                bytes.as_ptr() as *const c_char,
                new_header.as_mut_ptr(),
                bytes.len().min(MAX_LINE),
            );
        }
    }

    if num_seq == 0 {
        eprintln!("\nError:  no input sequences to analyze.\n");
        return 18;
    }

    handles.close();
    0
}

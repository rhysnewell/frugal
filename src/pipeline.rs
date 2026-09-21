/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    Rust port of main.c — the pipeline orchestration module.
*******************************************************************************/

use std::ffi::CStr;
use std::mem::MaybeUninit;
use std::os::raw::{c_char, c_int, c_void};

use crate::types::{Gene, Mask, MetagenomicBin, Node, Training};
use crate::types::{MAX_GENES, MAX_LINE, MAX_MASKS, MAX_SEQ, NUM_META, STT_NOD};

const VERSION: &str = "2.6.3";
const DATE: &str = "February, 2016";
const VERSION_CSTR: &[u8] = b"2.6.3\0";
const MIN_SINGLE_GENOME: c_int = 20000;
const IDEAL_SINGLE_GENOME: c_int = 100000;

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
use crate::gene::{
    add_genes, print_genes, record_gene_data, tweak_final_starts, write_nucleotide_seqs,
    write_translations,
};
use crate::metagenomic::initialize_metagenomic_metadata;
use crate::node::{
    add_nodes, calc_dicodon_gene, determine_sd_usage, raw_coding_score, rbs_score, record_gc_bias,
    record_overlapping_starts, reset_node_scores, score_nodes, train_starts_nonsd, train_starts_sd,
    write_start_file,
};
use crate::output::{close_handle, create_file, stdout_handle};
use crate::reader::{seq_reader_close, seq_reader_open, seq_reader_open_stdin, seq_reader_seek};
use crate::sequence::{
    calc_most_gc_frame, calc_short_header, next_seq_multi, rcom_seq, read_seq_training,
};
use crate::training::{read_training_file, write_training_file};

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
    // Convert config file paths to CStrings for FFI
    let input_cstr = config
        .input_file
        .as_ref()
        .map(|s| std::ffi::CString::new(s.as_str()).unwrap());
    let output_cstr = config
        .output_file
        .as_ref()
        .map(|s| std::ffi::CString::new(s.as_str()).unwrap());
    let trans_cstr = config
        .trans_file
        .as_ref()
        .map(|s| std::ffi::CString::new(s.as_str()).unwrap());
    let nuc_cstr = config
        .nuc_file
        .as_ref()
        .map(|s| std::ffi::CString::new(s.as_str()).unwrap());
    let start_cstr = config
        .start_file
        .as_ref()
        .map(|s| std::ffi::CString::new(s.as_str()).unwrap());
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
    let gc_frame: *mut c_int;
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
    let mut max_slen: c_int;
    let mut max_score: f64;
    let mut gc: f64 = 0.0;
    let mut low: f64;
    let mut high: f64;

    // Raw C pointers for the extern "C" file path args
    let train_file: *const c_char = train_cstr.as_ref().map_or(std::ptr::null(), |c| c.as_ptr());
    let start_file: *const c_char = start_cstr.as_ref().map_or(std::ptr::null(), |c| c.as_ptr());
    let trans_file: *const c_char = trans_cstr.as_ref().map_or(std::ptr::null(), |c| c.as_ptr());
    let nuc_file: *const c_char = nuc_cstr.as_ref().map_or(std::ptr::null(), |c| c.as_ptr());
    let input_file: *const c_char = input_cstr.as_ref().map_or(std::ptr::null(), |c| c.as_ptr());
    let output_file: *const c_char = output_cstr
        .as_ref()
        .map_or(std::ptr::null(), |c| c.as_ptr());

    let mut cur_header: [c_char; MAX_LINE] = [0; MAX_LINE];
    let mut new_header: [c_char; MAX_LINE] = [0; MAX_LINE];
    let mut short_header: [c_char; MAX_LINE] = [0; MAX_LINE];

    let mut output_ptr: c_int;
    let mut start_ptr: c_int;
    let mut trans_ptr: c_int;
    let mut nuc_ptr: c_int;
    let mut input_ptr: *mut c_void = std::ptr::null_mut();

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
    max_score = -100.0;
    do_training = 0;
    max_slen = 0;

    // Set up file descriptors (stdout = fd 1)
    let stdout_fd: c_int = stdout_handle();
    start_ptr = stdout_fd;
    trans_ptr = stdout_fd;
    nuc_ptr = stdout_fd;
    output_ptr = stdout_fd;

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

    // Check i/o files and prepare them for reading/writing
    if !input_file.is_null() {
        input_ptr = seq_reader_open(input_file);
        if input_ptr.is_null() {
            let f = CStr::from_ptr(input_file).to_string_lossy();
            eprintln!("\nError: can't open input file {}.\n", f);
            return 5;
        }
    }
    if input_ptr.is_null() {
        input_ptr = seq_reader_open_stdin();
        if input_ptr.is_null() {
            eprintln!("\nError: can't open stdin.\n");
            return 5;
        }
    }
    if !output_file.is_null() {
        let path = CStr::from_ptr(output_file).to_string_lossy();
        match create_file(path.as_ref()) {
            Ok(handle) => {
                output_ptr = handle;
            }
            Err(_) => {
                eprintln!("\nError: can't open output file {}.\n", path);
                return 6;
            }
        }
    }
    if !start_file.is_null() {
        let path = CStr::from_ptr(start_file).to_string_lossy();
        match create_file(path.as_ref()) {
            Ok(handle) => {
                start_ptr = handle;
            }
            Err(_) => {
                eprintln!("\nError: can't open start file {}.\n", path);
                return 7;
            }
        }
    }
    if !trans_file.is_null() {
        let path = CStr::from_ptr(trans_file).to_string_lossy();
        match create_file(path.as_ref()) {
            Ok(handle) => {
                trans_ptr = handle;
            }
            Err(_) => {
                eprintln!("\nError: can't open translation file {}.\n", path);
                return 8;
            }
        }
    }
    if !nuc_file.is_null() {
        let path = CStr::from_ptr(nuc_file).to_string_lossy();
        match create_file(path.as_ref()) {
            Ok(handle) => {
                nuc_ptr = handle;
            }
            Err(_) => {
                eprintln!("\nError: can't open gene nucleotide file {}.\n", path);
                return 16;
            }
        }
    }

    // =========================================================================
    // Single Genome Training
    // =========================================================================
    if is_meta == 0 && (do_training == 1 || (do_training == 0 && train_file.is_null())) {
        if quiet == 0 {
            eprintln!("Request:  Single Genome, Phase:  Training");
            eprint!("Reading in the sequence(s) to train...");
        }
        slen = read_seq_training(
            input_ptr,
            seq,
            useq,
            &mut tinf.gc,
            do_mask,
            mlist.as_mut_ptr(),
            &mut nmask,
        );
        if slen == 0 {
            eprint!("\n\nSequence read failed (file must be Fasta, ");
            eprintln!("Genbank, or EMBL format).\n");
            return 9;
        }
        if slen < MIN_SINGLE_GENOME {
            eprint!("\n\nError:  Sequence must be {}", MIN_SINGLE_GENOME);
            eprint!(" characters (only {} read).\n(Consider", slen);
            eprint!(" running with the -p meta option or finding");
            eprintln!(" more contigs from the same genome.)\n");
            return 10;
        }
        if slen < IDEAL_SINGLE_GENOME {
            eprint!("\n\nWarning:  ideally Prodigal should be given at");
            eprint!(" least {} bases for ", IDEAL_SINGLE_GENOME);
            eprint!("training.\nYou may get better results with the ");
            eprintln!("-p meta option.\n");
        }
        rcom_seq(seq, rseq, useq, slen);
        if quiet == 0 {
            eprintln!("{} bp seq created, {:.2} pct GC", slen, tinf.gc * 100.0);
        }

        // Find all potential starts and stops
        if quiet == 0 {
            eprint!("Locating all potential starts and stops...");
        }
        if slen > max_slen && slen > (STT_NOD as c_int) * 8 {
            let new_size = slen as usize / 8;
            if nodes_vec.capacity() < new_size {
                nodes_vec.reserve_exact(new_size - nodes_vec.capacity());
            }
            nodes = nodes_vec.as_mut_ptr() as *mut Node;
            max_slen = slen;
        }
        nn = add_nodes(
            seq,
            rseq,
            slen,
            nodes,
            closed,
            mlist.as_mut_ptr(),
            nmask,
            &mut tinf,
        );
        std::slice::from_raw_parts_mut(nodes, nn as usize)
            .sort_unstable_by(|a, b| a.ndx.cmp(&b.ndx).then(b.strand.cmp(&a.strand)));
        if quiet == 0 {
            eprintln!("{} nodes", nn);
        }

        // Scan ORFs for GC bias
        if quiet == 0 {
            eprint!("Looking for GC bias in different frames...");
        }
        gc_frame = calc_most_gc_frame(seq, slen);
        if gc_frame.is_null() {
            eprintln!("Malloc failed on gc frame plot\n");
            return 11;
        }
        record_gc_bias(gc_frame, nodes, nn, &mut tinf);
        if quiet == 0 {
            eprintln!(
                "frame bias scores: {:.2} {:.2} {:.2}",
                tinf.bias[0], tinf.bias[1], tinf.bias[2]
            );
        }
        drop(Vec::from_raw_parts(gc_frame, slen as usize, slen as usize));

        // Initial DP with GC frame bias
        if quiet == 0 {
            eprint!("Building initial set of genes to train from...");
        }
        record_overlapping_starts(nodes, nn, &mut tinf, 0);
        ipath = dprog(nodes, nn, &mut tinf, 0);
        if quiet == 0 {
            eprintln!("done!");
        }

        // Gather dicodon statistics
        if quiet == 0 {
            eprint!("Creating coding model and scoring nodes...");
        }
        calc_dicodon_gene(&mut tinf, seq, rseq, slen, nodes, ipath);
        raw_coding_score(seq, rseq, slen, nodes, nn, &mut tinf);
        if quiet == 0 {
            eprintln!("done!");
        }

        // Determine SD usage and train starts
        if quiet == 0 {
            eprint!("Examining upstream regions and training starts...");
        }
        rbs_score(seq, rseq, slen, nodes, nn, &mut tinf);
        train_starts_sd(seq, rseq, slen, nodes, nn, &mut tinf);
        determine_sd_usage(&mut tinf);
        if force_nonsd == 1 {
            tinf.uses_sd = 0;
        }
        if tinf.uses_sd == 0 {
            train_starts_nonsd(seq, rseq, slen, nodes, nn, &mut tinf);
        }
        if quiet == 0 {
            eprintln!("done!");
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
                // Clean up and return
                // Vec memory freed automatically when dropped
                seq_reader_close(input_ptr);
                if output_ptr != stdout_fd {
                    close_handle(output_ptr);
                }
                if start_ptr != stdout_fd {
                    close_handle(start_ptr);
                }
                if trans_ptr != stdout_fd {
                    close_handle(trans_ptr);
                }
                if nuc_ptr != stdout_fd {
                    close_handle(nuc_ptr);
                }
                return 0;
            }
        }

        // Rewind input file
        if quiet == 0 {
            eprintln!("-------------------------------------");
        }
        if seq_reader_seek(input_ptr, 0, 0) == -1 {
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
            input_ptr,
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

        // Reallocate if this is the biggest sequence we've seen
        if slen > max_slen && slen > (STT_NOD as c_int) * 8 {
            let new_size = slen as usize / 8;
            if nodes_vec.capacity() < new_size {
                nodes_vec.reserve_exact(new_size - nodes_vec.capacity());
            }
            nodes = nodes_vec.as_mut_ptr() as *mut Node;
            max_slen = slen;
        }

        // Calculate short header
        calc_short_header(cur_header.as_mut_ptr(), short_header.as_mut_ptr(), num_seq);

        if is_meta == 0 {
            // ---- Single Genome Version ----
            nn = add_nodes(
                seq,
                rseq,
                slen,
                nodes,
                closed,
                mlist.as_mut_ptr(),
                nmask,
                &mut tinf,
            );
            std::slice::from_raw_parts_mut(nodes, nn as usize)
                .sort_unstable_by(|a, b| a.ndx.cmp(&b.ndx).then(b.strand.cmp(&a.strand)));

            score_nodes(seq, rseq, slen, nodes, nn, &mut tinf, closed, is_meta);
            if start_ptr != stdout_fd {
                write_start_file(
                    start_ptr,
                    nodes,
                    nn,
                    &mut tinf,
                    num_seq,
                    slen,
                    0,
                    std::ptr::null_mut(),
                    VERSION_CSTR.as_ptr() as *mut c_char,
                    cur_header.as_mut_ptr(),
                );
            }
            record_overlapping_starts(nodes, nn, &mut tinf, 1);
            ipath = dprog(nodes, nn, &mut tinf, 1);
            eliminate_bad_genes(nodes, ipath, &mut tinf);
            ng = add_genes(genes, nodes, ipath);
            tweak_final_starts(genes, ng, nodes, nn, &mut tinf);
            record_gene_data(genes, ng, nodes, &mut tinf, num_seq);
            if quiet == 0 {
                eprintln!("done!");
            }

            // Output the genes
            print_genes(
                output_ptr,
                genes,
                ng,
                nodes,
                slen,
                output,
                num_seq,
                0,
                std::ptr::null_mut(),
                &mut tinf,
                cur_header.as_mut_ptr(),
                short_header.as_mut_ptr(),
                VERSION_CSTR.as_ptr() as *mut c_char,
            );
            // fd writes go directly to kernel, no fflush needed
            if trans_ptr != stdout_fd {
                write_translations(
                    trans_ptr,
                    genes,
                    ng,
                    nodes,
                    seq,
                    rseq,
                    useq,
                    slen,
                    &mut tinf,
                    num_seq,
                    short_header.as_mut_ptr(),
                );
            }
            if nuc_ptr != stdout_fd {
                write_nucleotide_seqs(
                    nuc_ptr,
                    genes,
                    ng,
                    nodes,
                    seq,
                    rseq,
                    useq,
                    slen,
                    &mut tinf,
                    num_seq,
                    short_header.as_mut_ptr(),
                );
            }
        } else {
            // ---- Metagenomic Version ----
            low = 0.88495 * gc - 0.0102337;
            if low > 0.65 {
                low = 0.65;
            }
            high = 0.86596 * gc + 0.1131991;
            if high < 0.35 {
                high = 0.35;
            }

            max_score = -100.0;
            ng = 0;
            for mi in 0..NUM_META as c_int {
                if mi == 0 || meta_trans_table[mi as usize] != meta_trans_table[(mi - 1) as usize] {
                    meta_tinf.trans_table = meta_trans_table[mi as usize];
                    std::ptr::write_bytes(nodes, 0, nn as usize);
                    nn = add_nodes(
                        seq,
                        rseq,
                        slen,
                        nodes,
                        closed,
                        mlist.as_mut_ptr(),
                        nmask,
                        &mut meta_tinf,
                    );
                    std::slice::from_raw_parts_mut(nodes, nn as usize)
                        .sort_unstable_by(|a, b| a.ndx.cmp(&b.ndx).then(b.strand.cmp(&a.strand)));
                }
                if meta_gc[mi as usize] < low || meta_gc[mi as usize] > high {
                    continue;
                }
                crate::training_data::load_metagenome(mi as usize, &mut meta_tinf);
                reset_node_scores(nodes, nn);
                score_nodes(seq, rseq, slen, nodes, nn, &mut meta_tinf, closed, is_meta);
                record_overlapping_starts(nodes, nn, &mut meta_tinf, 1);
                ipath = dprog(nodes, nn, &mut meta_tinf, 1);
                if ipath >= 0 && (*nodes.offset(ipath as isize)).score > max_score {
                    max_phase = mi;
                    max_score = (*nodes.offset(ipath as isize)).score;
                    eliminate_bad_genes(nodes, ipath, &mut meta_tinf);
                    ng = add_genes(genes, nodes, ipath);
                    tweak_final_starts(genes, ng, nodes, nn, &mut meta_tinf);
                    record_gene_data(genes, ng, nodes, &mut meta_tinf, num_seq);
                }
            }

            // Recover the nodes for the best run
            crate::training_data::load_metagenome(max_phase as usize, &mut best_tinf);
            std::ptr::write_bytes(nodes, 0, nn as usize);
            nn = add_nodes(
                seq,
                rseq,
                slen,
                nodes,
                closed,
                mlist.as_mut_ptr(),
                nmask,
                &mut best_tinf,
            );
            std::slice::from_raw_parts_mut(nodes, nn as usize)
                .sort_unstable_by(|a, b| a.ndx.cmp(&b.ndx).then(b.strand.cmp(&a.strand)));
            score_nodes(seq, rseq, slen, nodes, nn, &mut best_tinf, closed, is_meta);
            if start_ptr != stdout_fd {
                write_start_file(
                    start_ptr,
                    nodes,
                    nn,
                    &mut best_tinf,
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

            // Output the genes
            print_genes(
                output_ptr,
                genes,
                ng,
                nodes,
                slen,
                output,
                num_seq,
                1,
                meta[max_phase as usize].desc.as_mut_ptr(),
                &mut best_tinf,
                cur_header.as_mut_ptr(),
                short_header.as_mut_ptr(),
                VERSION_CSTR.as_ptr() as *mut c_char,
            );
            // fd writes go directly to kernel, no fflush needed
            if trans_ptr != stdout_fd {
                write_translations(
                    trans_ptr,
                    genes,
                    ng,
                    nodes,
                    seq,
                    rseq,
                    useq,
                    slen,
                    &mut best_tinf,
                    num_seq,
                    short_header.as_mut_ptr(),
                );
            }
            if nuc_ptr != stdout_fd {
                write_nucleotide_seqs(
                    nuc_ptr,
                    genes,
                    ng,
                    nodes,
                    seq,
                    rseq,
                    useq,
                    slen,
                    &mut best_tinf,
                    num_seq,
                    short_header.as_mut_ptr(),
                );
            }
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

    // Close all filehandles
    seq_reader_close(input_ptr);
    if output_ptr != stdout_fd {
        close_handle(output_ptr);
    }
    if start_ptr != stdout_fd {
        close_handle(start_ptr);
    }
    if trans_ptr != stdout_fd {
        close_handle(trans_ptr);
    }
    if nuc_ptr != stdout_fd {
        close_handle(nuc_ptr);
    }

    0
}

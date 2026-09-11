/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    Rust translation maintaining exact C ABI and output compatibility.
*******************************************************************************/

use std::ffi::CStr;
use std::os::raw::{c_char, c_int};

use crate::types::{Gene, Node, Training, MAX_GENES, MAX_SAM_OVLP, STOP};

use crate::node::intergenic_mod;
use crate::sequence::{amino, is_a, is_c, is_g, is_n, is_t, mer_text};

/// Write a Rust string to a file descriptor.
#[inline]
unsafe fn fprint(fd: c_int, s: &str) {
    let _ = crate::output::write_to_handle(fd, s);
}

/// Convert a *const c_char (C string) to a &str.  Returns "" on null.
#[inline]
unsafe fn cstr(p: *const c_char) -> &'static str {
    if p.is_null() {
        ""
    } else {
        CStr::from_ptr(p).to_str().unwrap_or("")
    }
}

/// Copy a Rust string into a c_char buffer, NUL-terminated.
#[inline]
unsafe fn copy_to_cbuf(buf: &mut [c_char], s: &str) {
    let bytes = s.as_bytes();
    let len = bytes.len().min(buf.len() - 1);
    for i in 0..len {
        buf[i] = bytes[i] as c_char;
    }
    buf[len] = 0;
}

/// Append a Rust string to an existing NUL-terminated c_char buffer.
#[inline]
unsafe fn append_to_cbuf(buf: &mut [c_char], s: &str) {
    // Find current NUL position
    let mut pos = 0usize;
    while pos < buf.len() && buf[pos] != 0 {
        pos += 1;
    }
    let bytes = s.as_bytes();
    let len = bytes.len().min(buf.len() - 1 - pos);
    for i in 0..len {
        buf[pos + i] = bytes[i] as c_char;
    }
    buf[pos + len] = 0;
}

#[inline]
unsafe fn clear_gene(glist: *mut Gene, ctr: c_int) {
    std::ptr::write(glist.offset(ctr as isize), std::mem::zeroed());
}

/// Copies genes from the dynamic programming to a final array.
pub unsafe fn add_genes(glist: *mut Gene, nod: *mut Node, dbeg: c_int) -> c_int {
    if dbeg == -1 {
        return 0;
    }

    let mut path = dbeg;
    let mut ctr: c_int = 0;

    while (*nod.offset(path as isize)).traceb != -1 {
        path = (*nod.offset(path as isize)).traceb;
    }

    while path != -1 {
        let np = &*nod.offset(path as isize);
        if np.elim == 1 {
            path = np.tracef;
            continue;
        }
        if np.strand == 1 && np.type_ != STOP {
            clear_gene(glist, ctr);
            (*glist.offset(ctr as isize)).begin = np.ndx + 1;
            (*glist.offset(ctr as isize)).start_ndx = path;
        }
        if np.strand == -1 && np.type_ == STOP {
            clear_gene(glist, ctr);
            (*glist.offset(ctr as isize)).begin = np.ndx - 1;
            (*glist.offset(ctr as isize)).stop_ndx = path;
        }
        if np.strand == 1 && np.type_ == STOP {
            (*glist.offset(ctr as isize)).end = np.ndx + 3;
            (*glist.offset(ctr as isize)).stop_ndx = path;
            ctr += 1;
        }
        if np.strand == -1 && np.type_ != STOP {
            (*glist.offset(ctr as isize)).end = np.ndx + 1;
            (*glist.offset(ctr as isize)).start_ndx = path;
            ctr += 1;
        }
        path = (*nod.offset(path as isize)).tracef;
        if ctr == MAX_GENES as c_int {
            eprint!("warning, max # of genes exceeded, truncating...\n");
            return ctr;
        }
    }
    ctr
}

/// Attempts to solve the problem of extremely close starts. If two potential
/// starts are 5 amino acids or less away from each other, this routine sets
/// their coding equal to each other and lets the RBS/operon/ATG-GTG-TTG start
/// features determine which start to use, under the assumption that 5 or less
/// words of coding is too weak a signal to use to select the proper start.
///
/// In addition, attempts to correct TTG (or whatever start codon is rare)
/// starts that have an RBS score, an upstream score, and a coding score all
/// superior to whatever start was initially chosen.
///
/// This routine was tested on numerous genomes and found to increase overall
/// performance.
pub unsafe fn tweak_final_starts(
    genes: *mut Gene,
    ng: c_int,
    nod: *mut Node,
    nn: c_int,
    tinf: *const Training,
) {
    let mut maxndx: [c_int; 2];
    let mut maxsc: [f64; 2];
    let mut maxigm: [f64; 2];

    for i in 0..ng {
        let ndx = (*genes.offset(i as isize)).start_ndx;
        let sc = (*nod.offset(ndx as isize)).sscore + (*nod.offset(ndx as isize)).cscore;
        let mut igm: f64 = 0.0;

        if i > 0
            && (*nod.offset(ndx as isize)).strand == 1
            && (*nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize)).strand == 1
        {
            igm = intergenic_mod(
                nod.offset((*genes.offset((i - 1) as isize)).stop_ndx as isize),
                nod.offset(ndx as isize),
                tinf,
            );
        }
        if i > 0
            && (*nod.offset(ndx as isize)).strand == 1
            && (*nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize)).strand == -1
        {
            igm = intergenic_mod(
                nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize),
                nod.offset(ndx as isize),
                tinf,
            );
        }
        if i < ng - 1
            && (*nod.offset(ndx as isize)).strand == -1
            && (*nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize)).strand == 1
        {
            igm = intergenic_mod(
                nod.offset(ndx as isize),
                nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize),
                tinf,
            );
        }
        if i < ng - 1
            && (*nod.offset(ndx as isize)).strand == -1
            && (*nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize)).strand == -1
        {
            igm = intergenic_mod(
                nod.offset(ndx as isize),
                nod.offset((*genes.offset((i + 1) as isize)).stop_ndx as isize),
                tinf,
            );
        }

        // Search upstream and downstream for the #2 and #3 scoring starts
        maxndx = [-1, -1];
        maxsc = [0.0, 0.0];
        maxigm = [0.0, 0.0];

        for j in (ndx - 100)..(ndx + 100) {
            if j < 0 || j >= nn || j == ndx {
                continue;
            }
            let nj = &*nod.offset(j as isize);
            if nj.type_ == STOP || nj.stop_val != (*nod.offset(ndx as isize)).stop_val {
                continue;
            }

            let mut tigm: f64 = 0.0;
            if i > 0
                && nj.strand == 1
                && (*nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize)).strand == 1
            {
                if (*nod.offset((*genes.offset((i - 1) as isize)).stop_ndx as isize)).ndx - nj.ndx
                    > MAX_SAM_OVLP
                {
                    continue;
                }
                tigm = intergenic_mod(
                    nod.offset((*genes.offset((i - 1) as isize)).stop_ndx as isize),
                    nod.offset(j as isize),
                    tinf,
                );
            }
            if i > 0
                && nj.strand == 1
                && (*nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize)).strand == -1
            {
                if (*nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize)).ndx - nj.ndx
                    >= 0
                {
                    continue;
                }
                tigm = intergenic_mod(
                    nod.offset((*genes.offset((i - 1) as isize)).start_ndx as isize),
                    nod.offset(j as isize),
                    tinf,
                );
            }
            if i < ng - 1
                && nj.strand == -1
                && (*nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize)).strand == 1
            {
                if nj.ndx - (*nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize)).ndx
                    >= 0
                {
                    continue;
                }
                tigm = intergenic_mod(
                    nod.offset(j as isize),
                    nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize),
                    tinf,
                );
            }
            if i < ng - 1
                && nj.strand == -1
                && (*nod.offset((*genes.offset((i + 1) as isize)).start_ndx as isize)).strand == -1
            {
                if nj.ndx - (*nod.offset((*genes.offset((i + 1) as isize)).stop_ndx as isize)).ndx
                    > MAX_SAM_OVLP
                {
                    continue;
                }
                tigm = intergenic_mod(
                    nod.offset(j as isize),
                    nod.offset((*genes.offset((i + 1) as isize)).stop_ndx as isize),
                    tinf,
                );
            }

            if maxndx[0] == -1 {
                maxndx[0] = j;
                maxsc[0] = nj.cscore + nj.sscore;
                maxigm[0] = tigm;
            } else if nj.cscore + nj.sscore + tigm > maxsc[0] {
                maxndx[1] = maxndx[0];
                maxsc[1] = maxsc[0];
                maxigm[1] = maxigm[0];
                maxndx[0] = j;
                maxsc[0] = nj.cscore + nj.sscore;
                maxigm[0] = tigm;
            } else if maxndx[1] == -1 || nj.cscore + nj.sscore + tigm > maxsc[1] {
                maxndx[1] = j;
                maxsc[1] = nj.cscore + nj.sscore;
                maxigm[1] = tigm;
            }
        }

        // Change the start if it's a TTG with better coding/RBS/upstream score
        // Also change the start if it's <=15bp but has better coding/RBS
        for j in 0..2 {
            let mndx = maxndx[j];
            if mndx == -1 {
                continue;
            }

            let nm = &*nod.offset(mndx as isize);
            let nn_ref = &*nod.offset(ndx as isize);

            // Start of less common type but with better coding, rbs, and
            // upstream.  Must be 18 or more bases away from original.
            if nm.tscore < nn_ref.tscore
                && maxsc[j] - nm.tscore >= sc - nn_ref.tscore + (*tinf).st_wt
                && nm.rscore > nn_ref.rscore
                && nm.uscore > nn_ref.uscore
                && nm.cscore > nn_ref.cscore
                && (nm.ndx - nn_ref.ndx).abs() > 15
            {
                maxsc[j] += nn_ref.tscore - nm.tscore;
            }
            // Close starts.  Ignore coding and see if start has better rbs
            // and type.
            else if (nm.ndx - nn_ref.ndx).abs() <= 15
                && nm.rscore + nm.tscore > nn_ref.rscore + nn_ref.tscore
                && nn_ref.edge == 0
                && nm.edge == 0
            {
                if nn_ref.cscore > nm.cscore {
                    maxsc[j] += nn_ref.cscore - nm.cscore;
                }
                if nn_ref.uscore > nm.uscore {
                    maxsc[j] += nn_ref.uscore - nm.uscore;
                }
                if igm > maxigm[j] {
                    maxsc[j] += igm - maxigm[j];
                }
            } else {
                maxsc[j] = -1000.0;
            }
        }

        // Change the gene coordinates to the new maximum.
        let mut mndx: c_int = -1;
        for j in 0..2i32 {
            if maxndx[j as usize] == -1 {
                continue;
            }
            if mndx == -1 && maxsc[j as usize] + maxigm[j as usize] > sc + igm {
                mndx = j;
            } else if mndx >= 0
                && maxsc[j as usize] + maxigm[j as usize]
                    > maxsc[mndx as usize] + maxigm[mndx as usize]
            {
                mndx = j;
            }
        }
        if mndx != -1 && (*nod.offset(maxndx[mndx as usize] as isize)).strand == 1 {
            (*genes.offset(i as isize)).start_ndx = maxndx[mndx as usize];
            (*genes.offset(i as isize)).begin =
                (*nod.offset(maxndx[mndx as usize] as isize)).ndx + 1;
        } else if mndx != -1 && (*nod.offset(maxndx[mndx as usize] as isize)).strand == -1 {
            (*genes.offset(i as isize)).start_ndx = maxndx[mndx as usize];
            (*genes.offset(i as isize)).end = (*nod.offset(maxndx[mndx as usize] as isize)).ndx + 1;
        }
    }
}

/// Builds the `gene_data` and `score_data` annotation strings for each gene,
/// recording start type, partial-edge flags, RBS motif/spacer, GC content,
/// and the component scores (confidence, total, coding, start, RBS, upstream,
/// type) used in the final output.
pub unsafe fn record_gene_data(
    genes: *mut Gene,
    ng: c_int,
    nod: *mut Node,
    tinf: *const Training,
    sctr: c_int,
) {
    // SD motif string tables
    let sd_string: [&str; 28] = [
        "None",
        "GGA/GAG/AGG",
        "3Base/5BMM",
        "4Base/6BMM",
        "AGxAG",
        "AGxAG",
        "GGA/GAG/AGG",
        "GGxGG",
        "GGxGG",
        "AGxAG",
        "AGGAG(G)/GGAGG",
        "AGGA/GGAG/GAGG",
        "AGGA/GGAG/GAGG",
        "GGA/GAG/AGG",
        "GGxGG",
        "AGGA",
        "GGAG/GAGG",
        "AGxAGG/AGGxGG",
        "AGxAGG/AGGxGG",
        "AGxAGG/AGGxGG",
        "AGGAG/GGAGG",
        "AGGAG",
        "AGGAG",
        "GGAGG",
        "GGAGG",
        "AGGAGG",
        "AGGAGG",
        "AGGAGG",
    ];

    let sd_spacer: [&str; 28] = [
        "None", "3-4bp", "13-15bp", "13-15bp", "11-12bp", "3-4bp", "11-12bp", "11-12bp", "3-4bp",
        "5-10bp", "13-15bp", "3-4bp", "11-12bp", "5-10bp", "5-10bp", "5-10bp", "5-10bp", "11-12bp",
        "3-4bp", "5-10bp", "11-12bp", "3-4bp", "5-10bp", "3-4bp", "5-10bp", "11-12bp", "3-4bp",
        "5-10bp",
    ];

    let type_string: [&str; 4] = ["ATG", "GTG", "TTG", "Edge"];

    let mut qt: [c_char; 10] = [0; 10];

    for i in 0..ng {
        let gi = &mut *genes.offset(i as isize);
        let ndx = gi.start_ndx;
        let sndx = gi.stop_ndx;
        let n = &*nod.offset(ndx as isize);
        let sn = &*nod.offset(sndx as isize);

        // Record basic gene data
        let partial_left: c_int =
            if (n.edge == 1 && n.strand == 1) || (sn.edge == 1 && n.strand == -1) {
                1
            } else {
                0
            };
        let partial_right: c_int =
            if (sn.edge == 1 && n.strand == 1) || (n.edge == 1 && n.strand == -1) {
                1
            } else {
                0
            };
        let st_type: c_int = if n.edge == 1 { 3 } else { n.type_ };

        let s = format!(
            "ID={}_{};partial={}{};start_type={};",
            sctr,
            i + 1,
            partial_left,
            partial_right,
            type_string[st_type as usize]
        );
        copy_to_cbuf(&mut gi.gene_data, &s);

        // Record rbs data
        let rbs1 = (*tinf).rbs_wt[n.rbs[0] as usize] * (*tinf).st_wt;
        let rbs2 = (*tinf).rbs_wt[n.rbs[1] as usize] * (*tinf).st_wt;

        if (*tinf).uses_sd == 1 {
            if rbs1 > rbs2 {
                let s = format!(
                    "rbs_motif={};rbs_spacer={}",
                    sd_string[n.rbs[0] as usize], sd_spacer[n.rbs[0] as usize]
                );
                append_to_cbuf(&mut gi.gene_data, &s);
            } else {
                let s = format!(
                    "rbs_motif={};rbs_spacer={}",
                    sd_string[n.rbs[1] as usize], sd_spacer[n.rbs[1] as usize]
                );
                append_to_cbuf(&mut gi.gene_data, &s);
            }
        } else {
            mer_text(qt.as_mut_ptr(), n.mot.len, n.mot.ndx);
            let qt_str = cstr(qt.as_ptr());
            if (*tinf).no_mot > -0.5 && rbs1 > rbs2 && rbs1 > n.mot.score * (*tinf).st_wt {
                let s = format!(
                    "rbs_motif={};rbs_spacer={}",
                    sd_string[n.rbs[0] as usize], sd_spacer[n.rbs[0] as usize]
                );
                append_to_cbuf(&mut gi.gene_data, &s);
            } else if (*tinf).no_mot > -0.5 && rbs2 >= rbs1 && rbs2 > n.mot.score * (*tinf).st_wt {
                let s = format!(
                    "rbs_motif={};rbs_spacer={}",
                    sd_string[n.rbs[1] as usize], sd_spacer[n.rbs[1] as usize]
                );
                append_to_cbuf(&mut gi.gene_data, &s);
            } else if n.mot.len == 0 {
                append_to_cbuf(&mut gi.gene_data, "rbs_motif=None;rbs_spacer=None");
            } else {
                let s = format!("rbs_motif={};rbs_spacer={}bp", qt_str, n.mot.spacer);
                append_to_cbuf(&mut gi.gene_data, &s);
            }
        }
        let gc_s = format!(";gc_cont={:.3}", n.gc_cont);
        append_to_cbuf(&mut gi.gene_data, &gc_s);

        // Record score data
        let confidence = calculate_confidence(n.cscore + n.sscore, (*tinf).st_wt);
        let s = format!(
            "conf={:.2};score={:.2};cscore={:.2};sscore={:.2};rscore={:.2};uscore={:.2};",
            confidence,
            n.cscore + n.sscore,
            n.cscore,
            n.sscore,
            n.rscore,
            n.uscore
        );
        copy_to_cbuf(&mut gi.score_data, &s);

        let s = format!("tscore={:.2};", n.tscore);
        append_to_cbuf(&mut gi.score_data, &s);
    }
}

/// Print the genes. 'flag' indicates which format to use.
pub unsafe fn print_genes(
    fp: c_int,
    genes: *mut Gene,
    ng: c_int,
    nod: *mut Node,
    slen: c_int,
    flag: c_int,
    sctr: c_int,
    is_meta: c_int,
    mdesc: *mut c_char,
    tinf: *const Training,
    header: *mut c_char,
    short_hdr: *mut c_char,
    version: *mut c_char,
) {
    let header_str = cstr(header);
    let short_hdr_str = cstr(short_hdr);
    let version_str = cstr(version);

    // Initialize sequence data
    let seq_data = format!("seqnum={};seqlen={};seqhdr=\"{}\"", sctr, slen, header_str);

    // Initialize run data string
    let mut run_data;
    if is_meta == 0 {
        run_data = format!("version=Prodigal.v{};run_type=Single;", version_str);
        run_data.push_str("model=\"Ab initio\";");
    } else {
        let mdesc_str = cstr(mdesc);
        run_data = format!("version=Prodigal.v{};run_type=Metagenomic;", version_str);
        run_data.push_str(&format!("model=\"{}\";", mdesc_str));
    }
    run_data.push_str(&format!(
        "gc_cont={:.2};transl_table={};uses_sd={}",
        (*tinf).gc * 100.0,
        (*tinf).trans_table,
        (*tinf).uses_sd
    ));

    // Print the gff header once
    if flag == 3 && sctr == 1 {
        fprint(fp, "##gff-version  3\n");
    }

    // Print sequence/model information
    if flag == 0 {
        fprint(fp, &format!("DEFINITION  {};{}\n", seq_data, run_data));
        fprint(fp, "FEATURES             Location/Qualifiers\n");
    } else if flag != 1 {
        fprint(fp, &format!("# Sequence Data: {}\n", seq_data));
        fprint(fp, &format!("# Model Data: {}\n", run_data));
    }

    // Print the genes
    for i in 0..ng {
        let gi = &*genes.offset(i as isize);
        let ndx = gi.start_ndx;
        let sndx = gi.stop_ndx;
        let n = &*nod.offset(ndx as isize);
        let sn = &*nod.offset(sndx as isize);

        let gene_data_str = cstr(gi.gene_data.as_ptr());
        let score_data_str = cstr(gi.score_data.as_ptr());

        // Print the coordinates and data
        if n.strand == 1 {
            let left = if n.edge == 1 {
                format!("<{}", gi.begin)
            } else {
                format!("{}", gi.begin)
            };
            let right = if sn.edge == 1 {
                format!(">{}", gi.end)
            } else {
                format!("{}", gi.end)
            };

            if flag == 0 {
                fprint(fp, &format!("     CDS             {}..{}\n", left, right));
                fprint(fp, "                     ");
                fprint(
                    fp,
                    &format!("/note=\"{};{}\"\n", gene_data_str, score_data_str),
                );
            }
            if flag == 1 {
                fprint(
                    fp,
                    &format!(
                        "gene_prodigal={}|1|f|y|y|3|0|{}|{}|{}|{}|-1|-1|1.0\n",
                        i + 1,
                        gi.begin,
                        gi.end,
                        gi.begin,
                        gi.end
                    ),
                );
            }
            if flag == 2 {
                fprint(fp, &format!(">{}_{}_{}_{}\n", i + 1, gi.begin, gi.end, "+"));
            }
            if flag == 3 {
                fprint(
                    fp,
                    &format!(
                        "{}\tProdigal_v{}\tCDS\t{}\t{}\t{:.1}\t+\t0\t{};{}",
                        short_hdr_str,
                        version_str,
                        gi.begin,
                        gi.end,
                        n.cscore + n.sscore,
                        gene_data_str,
                        score_data_str
                    ),
                );
                fprint(fp, "\n");
            }
        } else {
            let left = if sn.edge == 1 {
                format!("<{}", gi.begin)
            } else {
                format!("{}", gi.begin)
            };
            let right = if n.edge == 1 {
                format!(">{}", gi.end)
            } else {
                format!("{}", gi.end)
            };

            if flag == 0 {
                fprint(
                    fp,
                    &format!("     CDS             complement({}..{})\n", left, right),
                );
                fprint(fp, "                     ");
                fprint(
                    fp,
                    &format!("/note=\"{};{}\"\n", gene_data_str, score_data_str),
                );
            }
            if flag == 1 {
                fprint(
                    fp,
                    &format!(
                        "gene_prodigal={}|1|r|y|y|3|0|{}|{}|{}|{}|-1|-1|1.0\n",
                        i + 1,
                        slen + 1 - gi.end,
                        slen + 1 - gi.begin,
                        slen + 1 - gi.end,
                        slen + 1 - gi.begin
                    ),
                );
            }
            if flag == 2 {
                fprint(fp, &format!(">{}_{}_{}_{}\n", i + 1, gi.begin, gi.end, "-"));
            }
            if flag == 3 {
                fprint(
                    fp,
                    &format!(
                        "{}\tProdigal_v{}\tCDS\t{}\t{}\t{:.1}\t-\t0\t{};{}",
                        short_hdr_str,
                        version_str,
                        gi.begin,
                        gi.end,
                        n.cscore + n.sscore,
                        gene_data_str,
                        score_data_str
                    ),
                );
                fprint(fp, "\n");
            }
        }
    }

    // Footer
    if flag == 0 {
        fprint(fp, "//\n");
    }
}

/// Print the gene translations.
pub unsafe fn write_translations(
    fh: c_int,
    genes: *mut Gene,
    ng: c_int,
    nod: *mut Node,
    seq: *mut u8,
    rseq: *mut u8,
    useq: *mut u8,
    slen: c_int,
    tinf: *const Training,
    _sctr: c_int,
    short_hdr: *mut c_char,
) {
    let short_hdr_str = cstr(short_hdr);

    for i in 0..ng {
        let gi = &*genes.offset(i as isize);
        let gene_data_str = cstr(gi.gene_data.as_ptr());
        if (*nod.offset(gi.start_ndx as isize)).strand == 1 {
            fprint(
                fh,
                &format!(
                    ">{}_{} # {} # {} # 1 # {}\n",
                    short_hdr_str,
                    i + 1,
                    gi.begin,
                    gi.end,
                    gene_data_str
                ),
            );
            let mut j = gi.begin;
            while j < gi.end {
                if is_n(useq, j - 1) == 1 || is_n(useq, j) == 1 || is_n(useq, j + 1) == 1 {
                    fprint(fh, "X");
                } else {
                    let is_init = if j == gi.begin { 1 } else { 0 };
                    let edge_val = 1 - (*nod.offset(gi.start_ndx as isize)).edge;
                    let ch = amino(seq, j - 1, tinf, is_init & edge_val) as u8 as char;
                    fprint(fh, &format!("{}", ch));
                }
                if (j - gi.begin) % 180 == 177 {
                    fprint(fh, "\n");
                }
                j += 3;
            }
            if (j - gi.begin) % 180 != 0 {
                fprint(fh, "\n");
            }
        } else {
            fprint(
                fh,
                &format!(
                    ">{}_{} # {} # {} # -1 # {}\n",
                    short_hdr_str,
                    i + 1,
                    gi.begin,
                    gi.end,
                    gene_data_str
                ),
            );
            let mut j = slen + 1 - gi.end;
            while j < slen + 1 - gi.begin {
                if is_n(useq, slen - j) == 1
                    || is_n(useq, slen - 1 - j) == 1
                    || is_n(useq, slen - 2 - j) == 1
                {
                    fprint(fh, "X");
                } else {
                    let is_init = if j == slen + 1 - gi.end { 1 } else { 0 };
                    let edge_val = 1 - (*nod.offset(gi.start_ndx as isize)).edge;
                    let ch = amino(rseq, j - 1, tinf, is_init & edge_val) as u8 as char;
                    fprint(fh, &format!("{}", ch));
                }
                if (j - slen - 1 + gi.end) % 180 == 177 {
                    fprint(fh, "\n");
                }
                j += 3;
            }
            if (j - slen - 1 + gi.end) % 180 != 0 {
                fprint(fh, "\n");
            }
        }
    }
}

/// Print the gene nucleotide sequences.
pub unsafe fn write_nucleotide_seqs(
    fh: c_int,
    genes: *mut Gene,
    ng: c_int,
    nod: *mut Node,
    seq: *mut u8,
    rseq: *mut u8,
    useq: *mut u8,
    slen: c_int,
    _tinf: *const Training,
    _sctr: c_int,
    short_hdr: *mut c_char,
) {
    let short_hdr_str = cstr(short_hdr);

    for i in 0..ng {
        let gi = &*genes.offset(i as isize);
        let gene_data_str = cstr(gi.gene_data.as_ptr());
        if (*nod.offset(gi.start_ndx as isize)).strand == 1 {
            fprint(
                fh,
                &format!(
                    ">{}_{} # {} # {} # 1 # {}\n",
                    short_hdr_str,
                    i + 1,
                    gi.begin,
                    gi.end,
                    gene_data_str
                ),
            );
            let mut j = gi.begin - 1;
            while j < gi.end {
                if is_a(seq, j) == 1 {
                    fprint(fh, "A");
                } else if is_t(seq, j) == 1 {
                    fprint(fh, "T");
                } else if is_g(seq, j) == 1 {
                    fprint(fh, "G");
                } else if is_c(seq, j) == 1 && is_n(useq, j) == 0 {
                    fprint(fh, "C");
                } else {
                    fprint(fh, "N");
                }
                if (j - gi.begin + 1) % 70 == 69 {
                    fprint(fh, "\n");
                }
                j += 1;
            }
            if (j - gi.begin + 1) % 70 != 0 {
                fprint(fh, "\n");
            }
        } else {
            fprint(
                fh,
                &format!(
                    ">{}_{} # {} # {} # -1 # {}\n",
                    short_hdr_str,
                    i + 1,
                    gi.begin,
                    gi.end,
                    gene_data_str
                ),
            );
            let mut j = slen - gi.end;
            while j < slen + 1 - gi.begin {
                if is_a(rseq, j) == 1 {
                    fprint(fh, "A");
                } else if is_t(rseq, j) == 1 {
                    fprint(fh, "T");
                } else if is_g(rseq, j) == 1 {
                    fprint(fh, "G");
                } else if is_c(rseq, j) == 1 && is_n(useq, slen - 1 - j) == 0 {
                    fprint(fh, "C");
                } else {
                    fprint(fh, "N");
                }
                if (j - slen + gi.end) % 70 == 69 {
                    fprint(fh, "\n");
                }
                j += 1;
            }
            if (j - slen + gi.end) % 70 != 0 {
                fprint(fh, "\n");
            }
        }
    }
}

/// Convert score to a percent confidence.
pub unsafe fn calculate_confidence(score: f64, start_weight: f64) -> f64 {
    let mut conf: f64;

    if score / start_weight < 41.0 {
        conf = (score / start_weight).exp();
        conf = (conf / (conf + 1.0)) * 100.0;
    } else {
        conf = 99.99;
    }
    if conf <= 50.00 {
        conf = 50.00;
    }
    conf
}

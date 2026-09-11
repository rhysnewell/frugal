/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

use std::ffi::CStr;
use std::os::raw::{c_char, c_int, c_void};

use crate::sequence::{
    is_atg, is_gc, is_gtg, is_start, is_stop, is_ttg, mer_ndx, shine_dalgarno_exact,
    shine_dalgarno_mm,
};
use crate::types::{
    Mask, Motif, Node, Training, ATG, EDGE_BONUS, EDGE_UPS, GTG, MAX_SAM_OVLP, META_PEN,
    MIN_EDGE_GENE, MIN_GENE, OPER_DIST, STOP, TTG,
};

use crate::sequence::{calc_mer_bg, max_fr, mer_text};
use crate::sequence::{shine_dalgarno_exact_mask, shine_dalgarno_mm_mask};

/// Write a Rust string to a file descriptor.
#[inline]
unsafe fn fprint(fp: c_int, s: &str) {
    let _ = crate::output::write_to_handle(fp, s);
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

#[inline]
unsafe fn clear_next_node(nodes: *mut Node, nn: c_int) -> *mut Node {
    let node = nodes.offset(nn as isize);
    std::ptr::write(node, std::mem::zeroed());
    node
}

/*******************************************************************************
  Adds nodes to the node list.  Genes must be >=90bp in length, unless they
  run off the edge, in which case they only have to be 50bp.
*******************************************************************************/

/// Adds nodes to the node list. Genes must be >=90bp in length, unless they
/// run off the edge, in which case they only have to be 50bp. Walks both the
/// forward and reverse strands, emitting start nodes (ATG/GTG/TTG) and stop
/// nodes in all three reading frames, respecting masked regions. Returns the
/// total number of nodes written.
pub unsafe fn add_nodes(
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nodes: *mut Node,
    closed: c_int,
    mlist: *mut Mask,
    nm: c_int,
    tinf: *mut Training,
) -> c_int {
    let mut nn: c_int = 0;
    let mut last: [c_int; 3] = [0; 3];
    let mut saw_start: [c_int; 3] = [0; 3];
    let mut min_dist: [c_int; 3] = [0; 3];
    let slmod = slen % 3;

    /* Forward strand nodes */
    for i in 0..3 {
        last[((i + slmod) % 3) as usize] = slen + i;
        saw_start[(i % 3) as usize] = 0;
        min_dist[(i % 3) as usize] = MIN_EDGE_GENE;
        if closed == 0 {
            while last[((i + slmod) % 3) as usize] + 2 > slen - 1 {
                last[((i + slmod) % 3) as usize] -= 3;
            }
        }
    }
    let mut i = slen - 3;
    while i >= 0 {
        if is_stop(seq, i, tinf) == 1 {
            if saw_start[(i % 3) as usize] == 1 {
                let node = clear_next_node(nodes, nn);
                if is_stop(seq, last[(i % 3) as usize], tinf) == 0 {
                    (*node).edge = 1;
                }
                (*node).ndx = last[(i % 3) as usize];
                (*node).type_ = STOP;
                (*node).strand = 1;
                (*node).stop_val = i;
                nn += 1;
            }
            min_dist[(i % 3) as usize] = MIN_GENE;
            last[(i % 3) as usize] = i;
            saw_start[(i % 3) as usize] = 0;
            i -= 1;
            continue;
        }
        if last[(i % 3) as usize] >= slen {
            i -= 1;
            continue;
        }

        if is_start(seq, i, tinf) == 1
            && is_atg(seq, i) == 1
            && (last[(i % 3) as usize] - i + 3) >= min_dist[(i % 3) as usize]
            && cross_mask(i, last[(i % 3) as usize], mlist, nm) == 0
        {
            let node = clear_next_node(nodes, nn);
            (*node).ndx = i;
            (*node).type_ = ATG;
            saw_start[(i % 3) as usize] = 1;
            (*node).stop_val = last[(i % 3) as usize];
            (*node).strand = 1;
            nn += 1;
        } else if is_start(seq, i, tinf) == 1
            && is_gtg(seq, i) == 1
            && (last[(i % 3) as usize] - i + 3) >= min_dist[(i % 3) as usize]
            && cross_mask(i, last[(i % 3) as usize], mlist, nm) == 0
        {
            let node = clear_next_node(nodes, nn);
            (*node).ndx = i;
            (*node).type_ = GTG;
            saw_start[(i % 3) as usize] = 1;
            (*node).stop_val = last[(i % 3) as usize];
            (*node).strand = 1;
            nn += 1;
        } else if is_start(seq, i, tinf) == 1
            && is_ttg(seq, i) == 1
            && (last[(i % 3) as usize] - i + 3) >= min_dist[(i % 3) as usize]
            && cross_mask(i, last[(i % 3) as usize], mlist, nm) == 0
        {
            let node = clear_next_node(nodes, nn);
            (*node).ndx = i;
            (*node).type_ = TTG;
            saw_start[(i % 3) as usize] = 1;
            (*node).stop_val = last[(i % 3) as usize];
            (*node).strand = 1;
            nn += 1;
        } else if i <= 2
            && closed == 0
            && (last[(i % 3) as usize] - i) > MIN_EDGE_GENE
            && cross_mask(i, last[(i % 3) as usize], mlist, nm) == 0
        {
            let node = clear_next_node(nodes, nn);
            (*node).ndx = i;
            (*node).type_ = ATG;
            saw_start[(i % 3) as usize] = 1;
            (*node).edge = 1;
            (*node).stop_val = last[(i % 3) as usize];
            (*node).strand = 1;
            nn += 1;
        }
        i -= 1;
    }
    for i in 0..3 {
        if saw_start[(i % 3) as usize] == 1 {
            let node = clear_next_node(nodes, nn);
            if is_stop(seq, last[(i % 3) as usize], tinf) == 0 {
                (*node).edge = 1;
            }
            (*node).ndx = last[(i % 3) as usize];
            (*node).type_ = STOP;
            (*node).strand = 1;
            (*node).stop_val = i - 6;
            nn += 1;
        }
    }

    /* Reverse strand nodes */
    for i in 0..3 {
        last[((i + slmod) % 3) as usize] = slen + i;
        saw_start[(i % 3) as usize] = 0;
        min_dist[(i % 3) as usize] = MIN_EDGE_GENE;
        if closed == 0 {
            while last[((i + slmod) % 3) as usize] + 2 > slen - 1 {
                last[((i + slmod) % 3) as usize] -= 3;
            }
        }
    }
    i = slen - 3;
    while i >= 0 {
        if is_stop(rseq, i, tinf) == 1 {
            if saw_start[(i % 3) as usize] == 1 {
                let node = clear_next_node(nodes, nn);
                if is_stop(rseq, last[(i % 3) as usize], tinf) == 0 {
                    (*node).edge = 1;
                }
                (*node).ndx = slen - last[(i % 3) as usize] - 1;
                (*node).type_ = STOP;
                (*node).strand = -1;
                (*node).stop_val = slen - i - 1;
                nn += 1;
            }
            min_dist[(i % 3) as usize] = MIN_GENE;
            last[(i % 3) as usize] = i;
            saw_start[(i % 3) as usize] = 0;
            i -= 1;
            continue;
        }
        if last[(i % 3) as usize] >= slen {
            i -= 1;
            continue;
        }

        if is_start(rseq, i, tinf) == 1
            && is_atg(rseq, i) == 1
            && (last[(i % 3) as usize] - i + 3) >= min_dist[(i % 3) as usize]
            && cross_mask(slen - last[(i % 3) as usize] - 1, slen - i - 1, mlist, nm) == 0
        {
            let node = clear_next_node(nodes, nn);
            (*node).ndx = slen - i - 1;
            (*node).type_ = ATG;
            saw_start[(i % 3) as usize] = 1;
            (*node).stop_val = slen - last[(i % 3) as usize] - 1;
            (*node).strand = -1;
            nn += 1;
        } else if is_start(rseq, i, tinf) == 1
            && is_gtg(rseq, i) == 1
            && (last[(i % 3) as usize] - i + 3) >= min_dist[(i % 3) as usize]
            && cross_mask(slen - last[(i % 3) as usize] - 1, slen - i - 1, mlist, nm) == 0
        {
            let node = clear_next_node(nodes, nn);
            (*node).ndx = slen - i - 1;
            (*node).type_ = GTG;
            saw_start[(i % 3) as usize] = 1;
            (*node).stop_val = slen - last[(i % 3) as usize] - 1;
            (*node).strand = -1;
            nn += 1;
        } else if is_start(rseq, i, tinf) == 1
            && is_ttg(rseq, i) == 1
            && (last[(i % 3) as usize] - i + 3) >= min_dist[(i % 3) as usize]
            && cross_mask(slen - last[(i % 3) as usize] - 1, slen - i - 1, mlist, nm) == 0
        {
            let node = clear_next_node(nodes, nn);
            (*node).ndx = slen - i - 1;
            (*node).type_ = TTG;
            saw_start[(i % 3) as usize] = 1;
            (*node).stop_val = slen - last[(i % 3) as usize] - 1;
            (*node).strand = -1;
            nn += 1;
        } else if i <= 2
            && closed == 0
            && (last[(i % 3) as usize] - i) > MIN_EDGE_GENE
            && cross_mask(slen - last[(i % 3) as usize] - 1, slen - i - 1, mlist, nm) == 0
        {
            let node = clear_next_node(nodes, nn);
            (*node).ndx = slen - i - 1;
            (*node).type_ = ATG;
            saw_start[(i % 3) as usize] = 1;
            (*node).edge = 1;
            (*node).stop_val = slen - last[(i % 3) as usize] - 1;
            (*node).strand = -1;
            nn += 1;
        }
        i -= 1;
    }
    for i in 0..3 {
        if saw_start[(i % 3) as usize] == 1 {
            let node = clear_next_node(nodes, nn);
            if is_stop(rseq, last[(i % 3) as usize], tinf) == 0 {
                (*node).edge = 1;
            }
            (*node).ndx = slen - last[(i % 3) as usize] - 1;
            (*node).type_ = STOP;
            (*node).strand = -1;
            (*node).stop_val = slen - i + 5;
            nn += 1;
        }
    }
    nn
}

/* Simple routine to zero out the node scores */

/// Simple routine to zero out the node scores: clears all per-node score,
/// trace, motif and RBS fields so a fresh scoring pass can be performed.
pub unsafe fn reset_node_scores(nod: *mut Node, nn: c_int) {
    for i in 0..nn as isize {
        for j in 0..3 {
            (*nod.offset(i)).star_ptr[j] = 0;
            (*nod.offset(i)).gc_score[j] = 0.0;
        }
        for j in 0..2 {
            (*nod.offset(i)).rbs[j] = 0;
        }
        (*nod.offset(i)).score = 0.0;
        (*nod.offset(i)).cscore = 0.0;
        (*nod.offset(i)).sscore = 0.0;
        (*nod.offset(i)).rscore = 0.0;
        (*nod.offset(i)).tscore = 0.0;
        (*nod.offset(i)).uscore = 0.0;
        (*nod.offset(i)).traceb = -1;
        (*nod.offset(i)).tracef = -1;
        (*nod.offset(i)).ov_mark = -1;
        (*nod.offset(i)).elim = 0;
        (*nod.offset(i)).gc_bias = 0;
        std::ptr::write_bytes(&mut (*nod.offset(i)).mot as *mut Motif, 0, 1);
    }
}

/*******************************************************************************
  Since dynamic programming can't go 'backwards', we have to record
  information about overlapping genes in order to build the models.
*******************************************************************************/

/// Since dynamic programming can't go 'backwards', we have to record
/// information about overlapping genes in order to build the models. So, for
/// example, in cases like 5'->3', 5'-3' overlapping on the same strand, we
/// record information about the 2nd 5' end under the first 3' end's
/// information. For every stop, we calculate and store all the best starts
/// that could be used in genes that overlap that 3' end.
pub unsafe fn record_overlapping_starts(
    nod: *mut Node,
    nn: c_int,
    tinf: *mut Training,
    flag: c_int,
) {
    let mut max_sc: f64;

    for i in 0..nn {
        for j in 0..3 {
            (*nod.offset(i as isize)).star_ptr[j] = -1;
        }
        if (*nod.offset(i as isize)).type_ != STOP || (*nod.offset(i as isize)).edge == 1 {
            continue;
        }
        if (*nod.offset(i as isize)).strand == 1 {
            max_sc = -100.0;
            let mut j = i + 3;
            loop {
                if j < 0 {
                    break;
                }
                if j >= nn || (*nod.offset(j as isize)).ndx > (*nod.offset(i as isize)).ndx + 2 {
                    j -= 1;
                    continue;
                }
                if (*nod.offset(j as isize)).ndx + MAX_SAM_OVLP < (*nod.offset(i as isize)).ndx {
                    break;
                }
                if (*nod.offset(j as isize)).strand == 1 && (*nod.offset(j as isize)).type_ != STOP
                {
                    if (*nod.offset(j as isize)).stop_val <= (*nod.offset(i as isize)).ndx {
                        j -= 1;
                        continue;
                    }
                    if flag == 0
                        && (*nod.offset(i as isize)).star_ptr
                            [((*nod.offset(j as isize)).ndx % 3) as usize]
                            == -1
                    {
                        (*nod.offset(i as isize)).star_ptr
                            [((*nod.offset(j as isize)).ndx % 3) as usize] = j;
                    } else if flag == 1
                        && ((*nod.offset(j as isize)).cscore
                            + (*nod.offset(j as isize)).sscore
                            + intergenic_mod(
                                &mut *nod.offset(i as isize),
                                &mut *nod.offset(j as isize),
                                tinf,
                            )
                            > max_sc)
                    {
                        (*nod.offset(i as isize)).star_ptr
                            [((*nod.offset(j as isize)).ndx % 3) as usize] = j;
                        max_sc = (*nod.offset(j as isize)).cscore
                            + (*nod.offset(j as isize)).sscore
                            + intergenic_mod(
                                &mut *nod.offset(i as isize),
                                &mut *nod.offset(j as isize),
                                tinf,
                            );
                    }
                }
                j -= 1;
            }
        } else {
            max_sc = -100.0;
            let mut j = i - 3;
            while j < nn {
                if j < 0 || (*nod.offset(j as isize)).ndx < (*nod.offset(i as isize)).ndx - 2 {
                    j += 1;
                    continue;
                }
                if (*nod.offset(j as isize)).ndx - MAX_SAM_OVLP > (*nod.offset(i as isize)).ndx {
                    break;
                }
                if (*nod.offset(j as isize)).strand == -1 && (*nod.offset(j as isize)).type_ != STOP
                {
                    if (*nod.offset(j as isize)).stop_val >= (*nod.offset(i as isize)).ndx {
                        j += 1;
                        continue;
                    }
                    if flag == 0
                        && (*nod.offset(i as isize)).star_ptr
                            [((*nod.offset(j as isize)).ndx % 3) as usize]
                            == -1
                    {
                        (*nod.offset(i as isize)).star_ptr
                            [((*nod.offset(j as isize)).ndx % 3) as usize] = j;
                    } else if flag == 1
                        && ((*nod.offset(j as isize)).cscore
                            + (*nod.offset(j as isize)).sscore
                            + intergenic_mod(
                                &mut *nod.offset(j as isize),
                                &mut *nod.offset(i as isize),
                                tinf,
                            )
                            > max_sc)
                    {
                        (*nod.offset(i as isize)).star_ptr
                            [((*nod.offset(j as isize)).ndx % 3) as usize] = j;
                        max_sc = (*nod.offset(j as isize)).cscore
                            + (*nod.offset(j as isize)).sscore
                            + intergenic_mod(
                                &mut *nod.offset(j as isize),
                                &mut *nod.offset(i as isize),
                                tinf,
                            );
                    }
                }
                j += 1;
            }
        }
    }
}

/*******************************************************************************
  This routine goes through all the ORFs and counts the relative frequency of
  the most common frame for G+C content.
*******************************************************************************/

/// This routine goes through all the ORFs and counts the relative frequency
/// of the most common frame for G+C content. In high GC genomes, this tends
/// to be the third position. In low GC genomes, this tends to be the first
/// position. Genes will be selected as a training set based on the nature of
/// this bias for this particular organism.
#[allow(unused_assignments)]
pub unsafe fn record_gc_bias(gc: *mut c_int, nod: *mut Node, nn: c_int, tinf: *mut Training) {
    let mut ctr: [[c_int; 3]; 3] = [[0; 3]; 3];
    let mut last: [c_int; 3] = [0; 3];
    let mut tot: f64 = 0.0;

    if nn == 0 {
        return;
    }
    for i in 0..3 {
        for j in 0..3 {
            ctr[i][j] = 0;
        }
    }
    let mut i = nn - 1;
    while i >= 0 {
        let fr = ((*nod.offset(i as isize)).ndx % 3) as usize;
        let frmod = 3 - fr as c_int;
        if (*nod.offset(i as isize)).strand == 1 && (*nod.offset(i as isize)).type_ == STOP {
            for j in 0..3 {
                ctr[fr][j] = 0;
            }
            last[fr] = (*nod.offset(i as isize)).ndx;
            ctr[fr][((*gc.offset((*nod.offset(i as isize)).ndx as isize) + frmod) % 3) as usize] =
                1;
        } else if (*nod.offset(i as isize)).strand == 1 {
            let mut j = last[fr] - 3;
            while j >= (*nod.offset(i as isize)).ndx {
                ctr[fr][((*gc.offset(j as isize) + frmod) % 3) as usize] += 1;
                j -= 3;
            }
            let mfr = max_fr(ctr[fr][0], ctr[fr][1], ctr[fr][2]);
            (*nod.offset(i as isize)).gc_bias = mfr;
            for j in 0..3 {
                (*nod.offset(i as isize)).gc_score[j] = 3.0 * ctr[fr][j] as f64;
                (*nod.offset(i as isize)).gc_score[j] /= 1.0
                    * ((*nod.offset(i as isize)).stop_val - (*nod.offset(i as isize)).ndx + 3)
                        as f64;
            }
            last[fr] = (*nod.offset(i as isize)).ndx;
        }
        i -= 1;
    }
    for i in 0..nn {
        let fr = ((*nod.offset(i as isize)).ndx % 3) as usize;
        let frmod = fr as c_int;
        if (*nod.offset(i as isize)).strand == -1 && (*nod.offset(i as isize)).type_ == STOP {
            for j in 0..3 {
                ctr[fr][j] = 0;
            }
            last[fr] = (*nod.offset(i as isize)).ndx;
            ctr[fr]
                [((3 - *gc.offset((*nod.offset(i as isize)).ndx as isize) + frmod) % 3) as usize] =
                1;
        } else if (*nod.offset(i as isize)).strand == -1 {
            let mut j = last[fr] + 3;
            while j <= (*nod.offset(i as isize)).ndx {
                ctr[fr][((3 - *gc.offset(j as isize) + frmod) % 3) as usize] += 1;
                j += 3;
            }
            let mfr = max_fr(ctr[fr][0], ctr[fr][1], ctr[fr][2]);
            (*nod.offset(i as isize)).gc_bias = mfr;
            for j in 0..3 {
                (*nod.offset(i as isize)).gc_score[j] = 3.0 * ctr[fr][j] as f64;
                (*nod.offset(i as isize)).gc_score[j] /= 1.0
                    * ((*nod.offset(i as isize)).ndx - (*nod.offset(i as isize)).stop_val + 3)
                        as f64;
            }
            last[fr] = (*nod.offset(i as isize)).ndx;
        }
    }

    for i in 0..3 {
        (*tinf).bias[i] = 0.0;
    }
    for i in 0..nn {
        if (*nod.offset(i as isize)).type_ != STOP {
            let len =
                ((*nod.offset(i as isize)).stop_val - (*nod.offset(i as isize)).ndx).abs() + 1;
            (*tinf).bias[(*nod.offset(i as isize)).gc_bias as usize] += (*nod.offset(i as isize))
                .gc_score[(*nod.offset(i as isize)).gc_bias as usize]
                * len as f64
                / 1000.0;
        }
    }
    tot = (*tinf).bias[0] + (*tinf).bias[1] + (*tinf).bias[2];
    for i in 0..3 {
        (*tinf).bias[i] *= 3.0 / tot;
    }
}

/*******************************************************************************
  Simple routine that calculates the dicodon frequency in genes and in the
  background, and then stores the log likelihood of each 6-mer relative to the
  background.
*******************************************************************************/

/// Simple routine that calculates the dicodon frequency in genes and in the
/// background, and then stores the log likelihood of each 6-mer relative to
/// the background. Walks the DP traceback starting from `dbeg` to accumulate
/// 6-mer counts inside genes, then writes log-odds values into
/// `tinf.gene_dc`, clamped to the range [-5.0, 5.0].
pub unsafe fn calc_dicodon_gene(
    tinf: *mut Training,
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nod: *mut Node,
    dbeg: c_int,
) {
    let mut counts: [c_int; 4096] = [0; 4096];
    let mut glob: c_int = 0;
    let mut prob: [f64; 4096] = [0.0; 4096];
    let mut bg: [f64; 4096] = [0.0; 4096];
    let mut left: c_int;
    let mut right: c_int;
    let mut in_gene: c_int;

    left = -1;
    right = -1;
    calc_mer_bg(6, seq, rseq, slen, bg.as_mut_ptr());
    let mut path = dbeg;
    in_gene = 0;
    while path != -1 {
        if (*nod.offset(path as isize)).strand == -1 && (*nod.offset(path as isize)).type_ != STOP {
            in_gene = -1;
            left = slen - (*nod.offset(path as isize)).ndx - 1;
        }
        if (*nod.offset(path as isize)).strand == 1 && (*nod.offset(path as isize)).type_ == STOP {
            in_gene = 1;
            right = (*nod.offset(path as isize)).ndx + 2;
        }
        if in_gene == -1
            && (*nod.offset(path as isize)).strand == -1
            && (*nod.offset(path as isize)).type_ == STOP
        {
            right = slen - (*nod.offset(path as isize)).ndx + 1;
            let mut i = left;
            while i < right - 5 {
                counts[mer_ndx(6, rseq, i) as usize] += 1;
                glob += 1;
                i += 3;
            }
            in_gene = 0;
        }
        if in_gene == 1
            && (*nod.offset(path as isize)).strand == 1
            && (*nod.offset(path as isize)).type_ != STOP
        {
            left = (*nod.offset(path as isize)).ndx;
            let mut i = left;
            while i < right - 5 {
                counts[mer_ndx(6, seq, i) as usize] += 1;
                glob += 1;
                i += 3;
            }
            in_gene = 0;
        }
        path = (*nod.offset(path as isize)).traceb;
    }
    for i in 0..4096 {
        prob[i] = (counts[i] as f64 * 1.0) / (glob as f64 * 1.0);
        if prob[i] == 0.0 && bg[i] != 0.0 {
            (*tinf).gene_dc[i] = -5.0;
        } else if bg[i] == 0.0 {
            (*tinf).gene_dc[i] = 0.0;
        } else {
            (*tinf).gene_dc[i] = (prob[i] / bg[i]).ln();
        }
        if (*tinf).gene_dc[i] > 5.0 {
            (*tinf).gene_dc[i] = 5.0;
        }
        if (*tinf).gene_dc[i] < -5.0 {
            (*tinf).gene_dc[i] = -5.0;
        }
    }
}

/*******************************************************************************
  calc_amino_bg - declared in header but not defined in original C code.
  Included as a stub for link compatibility.
*******************************************************************************/

/// Stub for an amino-acid-background routine declared in `node.h` but never
/// defined in the original C source. Kept here only for link compatibility;
/// performs no work.
pub unsafe fn calc_amino_bg(
    _tinf: *mut Training,
    _seq: *mut u8,
    _rseq: *mut u8,
    _slen: c_int,
    _nod: *mut Node,
    _nn: c_int,
) {
    /* Not implemented in original Prodigal */
}

/*******************************************************************************
  Scoring function for all the start nodes.
*******************************************************************************/

/// Scoring function for all the start nodes. This score has two factors:
/// (1) Coding, which is a composite of coding score and length, and
/// (2) Start score, which is a composite of RBS score and ATG/TTG/GTG.
/// Also applies edge-gene bonuses/penalties and a metagenomic coding
/// penalty when `is_meta` is set on short contigs.
pub unsafe fn score_nodes(
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nod: *mut Node,
    nn: c_int,
    tinf: *mut Training,
    closed: c_int,
    is_meta: c_int,
) {
    score_nodes_with_rbs(
        seq,
        rseq,
        slen,
        nod,
        nn,
        tinf,
        closed,
        is_meta,
        std::ptr::null(),
    )
}

#[allow(clippy::too_many_arguments)]
pub unsafe fn score_nodes_with_rbs(
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nod: *mut Node,
    nn: c_int,
    tinf: *mut Training,
    closed: c_int,
    is_meta: c_int,
    rbs_masks: *const u32,
) {
    let mut negf: f64;
    let mut posf: f64;
    let mut rbs1: f64;
    let mut rbs2: f64;
    let mut sd_score: f64;
    let mut edge_gene: f64;
    let mut min_meta_len: f64;

    /* Step 1: Calculate raw coding potential for every start-stop pair. */
    calc_orf_gc(seq, rseq, slen, nod, nn, tinf);
    raw_coding_score(seq, rseq, slen, nod, nn, tinf);

    /* Step 2: Calculate raw RBS Scores for every start node. */
    if (*tinf).uses_sd == 1 {
        if rbs_masks.is_null() {
            rbs_score(seq, rseq, slen, nod, nn, tinf);
        } else {
            rbs_score_from_masks(nod, nn, tinf, rbs_masks);
        }
    } else {
        for i in 0..nn {
            if (*nod.offset(i as isize)).type_ == STOP || (*nod.offset(i as isize)).edge == 1 {
                continue;
            }
            find_best_upstream_motif(tinf, seq, rseq, slen, &mut *nod.offset(i as isize), 2);
        }
    }

    /* Step 3: Score the start nodes */
    for i in 0..nn {
        if (*nod.offset(i as isize)).type_ == STOP {
            continue;
        }

        /* Does this gene run off the edge? */
        edge_gene = 0.0;
        if (*nod.offset(i as isize)).edge == 1 {
            edge_gene += 1.0;
        }
        if ((*nod.offset(i as isize)).strand == 1
            && is_stop(seq, (*nod.offset(i as isize)).stop_val, tinf) == 0)
            || ((*nod.offset(i as isize)).strand == -1
                && is_stop(rseq, slen - 1 - (*nod.offset(i as isize)).stop_val, tinf) == 0)
        {
            edge_gene += 1.0;
        }

        /* Edge Nodes : stops with no starts, give a small bonus */
        if (*nod.offset(i as isize)).edge == 1 {
            (*nod.offset(i as isize)).tscore = EDGE_BONUS * (*tinf).st_wt / edge_gene;
            (*nod.offset(i as isize)).uscore = 0.0;
            (*nod.offset(i as isize)).rscore = 0.0;
        } else {
            /* Type Score */
            (*nod.offset(i as isize)).tscore =
                (*tinf).type_wt[(*nod.offset(i as isize)).type_ as usize] * (*tinf).st_wt;

            /* RBS Motif Score */
            rbs1 = (*tinf).rbs_wt[(*nod.offset(i as isize)).rbs[0] as usize];
            rbs2 = (*tinf).rbs_wt[(*nod.offset(i as isize)).rbs[1] as usize];
            sd_score = dmax(rbs1, rbs2) * (*tinf).st_wt;
            if (*tinf).uses_sd == 1 {
                (*nod.offset(i as isize)).rscore = sd_score;
            } else {
                (*nod.offset(i as isize)).rscore =
                    (*tinf).st_wt * (*nod.offset(i as isize)).mot.score;
                if (*nod.offset(i as isize)).rscore < sd_score && (*tinf).no_mot > -0.5 {
                    (*nod.offset(i as isize)).rscore = sd_score;
                }
            }

            /* Upstream Score */
            if (*nod.offset(i as isize)).strand == 1 {
                score_upstream_composition(seq, slen, &mut *nod.offset(i as isize), tinf);
            } else {
                score_upstream_composition(rseq, slen, &mut *nod.offset(i as isize), tinf);
            }

            /* Penalize upstream score if choosing this start would stop
            the gene from running off the edge. */
            if closed == 0
                && (*nod.offset(i as isize)).ndx <= 2
                && (*nod.offset(i as isize)).strand == 1
            {
                (*nod.offset(i as isize)).uscore += EDGE_UPS * (*tinf).st_wt;
            } else if closed == 0
                && (*nod.offset(i as isize)).ndx >= slen - 3
                && (*nod.offset(i as isize)).strand == -1
            {
                (*nod.offset(i as isize)).uscore += EDGE_UPS * (*tinf).st_wt;
            } else if i < 500 && (*nod.offset(i as isize)).strand == 1 {
                let mut j = i - 1;
                while j >= 0 {
                    if (*nod.offset(j as isize)).edge == 1
                        && (*nod.offset(i as isize)).stop_val == (*nod.offset(j as isize)).stop_val
                    {
                        (*nod.offset(i as isize)).uscore += EDGE_UPS * (*tinf).st_wt;
                        break;
                    }
                    j -= 1;
                }
            } else if i >= nn - 500 && (*nod.offset(i as isize)).strand == -1 {
                let mut j = i + 1;
                while j < nn {
                    if (*nod.offset(j as isize)).edge == 1
                        && (*nod.offset(i as isize)).stop_val == (*nod.offset(j as isize)).stop_val
                    {
                        (*nod.offset(i as isize)).uscore += EDGE_UPS * (*tinf).st_wt;
                        break;
                    }
                    j += 1;
                }
            }
        }

        /* Convert starts at base 1 and slen to edge genes if closed = 0 */
        if (((*nod.offset(i as isize)).ndx <= 2 && (*nod.offset(i as isize)).strand == 1)
            || ((*nod.offset(i as isize)).ndx >= slen - 3
                && (*nod.offset(i as isize)).strand == -1))
            && (*nod.offset(i as isize)).edge == 0
            && closed == 0
        {
            edge_gene += 1.0;
            (*nod.offset(i as isize)).edge = 1;
            (*nod.offset(i as isize)).tscore = 0.0;
            (*nod.offset(i as isize)).uscore = EDGE_BONUS * (*tinf).st_wt / edge_gene;
            (*nod.offset(i as isize)).rscore = 0.0;
        }

        /* Penalize starts with no stop codon */
        if (*nod.offset(i as isize)).edge == 0 && edge_gene == 1.0 {
            (*nod.offset(i as isize)).uscore -= 0.5 * EDGE_BONUS * (*tinf).st_wt;
        }

        /* Penalize non-edge genes < 250bp */
        if edge_gene == 0.0
            && ((*nod.offset(i as isize)).ndx - (*nod.offset(i as isize)).stop_val).abs() < 250
        {
            negf = 250.0
                / ((*nod.offset(i as isize)).ndx - (*nod.offset(i as isize)).stop_val).abs() as f64;
            posf = ((*nod.offset(i as isize)).ndx - (*nod.offset(i as isize)).stop_val).abs()
                as f64
                / 250.0;
            if (*nod.offset(i as isize)).rscore < 0.0 {
                (*nod.offset(i as isize)).rscore *= negf;
            }
            if (*nod.offset(i as isize)).uscore < 0.0 {
                (*nod.offset(i as isize)).uscore *= negf;
            }
            if (*nod.offset(i as isize)).tscore < 0.0 {
                (*nod.offset(i as isize)).tscore *= negf;
            }
            if (*nod.offset(i as isize)).rscore > 0.0 {
                (*nod.offset(i as isize)).rscore *= posf;
            }
            if (*nod.offset(i as isize)).uscore > 0.0 {
                (*nod.offset(i as isize)).uscore *= posf;
            }
            if (*nod.offset(i as isize)).tscore > 0.0 {
                (*nod.offset(i as isize)).tscore *= posf;
            }
        }

        /* Coding Penalization in Metagenomic Fragments */
        if is_meta == 1
            && slen < 3000
            && edge_gene == 0.0
            && ((*nod.offset(i as isize)).cscore < 5.0
                || ((*nod.offset(i as isize)).ndx - (*nod.offset(i as isize)).stop_val).abs() < 120)
        {
            (*nod.offset(i as isize)).cscore -= META_PEN * dmax(0.0, (3000 - slen) as f64 / 2700.0);
        }

        /* Base Start Score */
        (*nod.offset(i as isize)).sscore = (*nod.offset(i as isize)).tscore
            + (*nod.offset(i as isize)).rscore
            + (*nod.offset(i as isize)).uscore;

        /* Penalize starts if coding is negative. */
        if (*nod.offset(i as isize)).cscore < 0.0 {
            if edge_gene > 0.0 && (*nod.offset(i as isize)).edge == 0 {
                if is_meta == 0 || slen > 1500 {
                    (*nod.offset(i as isize)).sscore -= (*tinf).st_wt;
                } else {
                    (*nod.offset(i as isize)).sscore -= 10.31 - 0.004 * slen as f64;
                }
            } else if is_meta == 1 && slen < 3000 && (*nod.offset(i as isize)).edge == 1 {
                min_meta_len = (slen as f64).sqrt() * 5.0;
                if ((*nod.offset(i as isize)).ndx - (*nod.offset(i as isize)).stop_val).abs() as f64
                    >= min_meta_len
                {
                    if (*nod.offset(i as isize)).cscore >= 0.0 {
                        (*nod.offset(i as isize)).cscore = -1.0;
                    }
                    (*nod.offset(i as isize)).sscore = 0.0;
                    (*nod.offset(i as isize)).uscore = 0.0;
                }
            } else {
                (*nod.offset(i as isize)).sscore -= 0.5;
            }
        } else if (*nod.offset(i as isize)).cscore < 5.0
            && is_meta == 1
            && ((*nod.offset(i as isize)).ndx - (*nod.offset(i as isize)).stop_val).abs() < 120
            && (*nod.offset(i as isize)).sscore < 0.0
        {
            (*nod.offset(i as isize)).sscore -= (*tinf).st_wt;
        }
    }
}

/* Calculate the GC Content for each start-stop pair */
/// Calculate the GC content for each start-stop pair and store it in
/// `nod[i].gc_cont`. Iterates each frame on both strands, accumulating
/// GC counts from the stop back to each start.
pub unsafe fn calc_orf_gc(
    seq: *mut u8,
    _rseq: *mut u8,
    _slen: c_int,
    nod: *mut Node,
    nn: c_int,
    _tinf: *mut Training,
) {
    let mut last: [c_int; 3] = [0; 3];
    let mut gc: [f64; 3] = [0.0; 3];
    let mut gsize: f64;

    /* Go through each start-stop pair and calculate the %GC of the gene */
    for k in 0..3 {
        gc[k] = 0.0;
    }
    let mut i = nn - 1;
    while i >= 0 {
        let fr = ((*nod.offset(i as isize)).ndx % 3) as usize;
        if (*nod.offset(i as isize)).strand == 1 && (*nod.offset(i as isize)).type_ == STOP {
            last[fr] = (*nod.offset(i as isize)).ndx;
            gc[fr] = is_gc(seq, (*nod.offset(i as isize)).ndx) as f64
                + is_gc(seq, (*nod.offset(i as isize)).ndx + 1) as f64
                + is_gc(seq, (*nod.offset(i as isize)).ndx + 2) as f64;
        } else if (*nod.offset(i as isize)).strand == 1 {
            let mut j = last[fr] - 3;
            while j >= (*nod.offset(i as isize)).ndx {
                gc[fr] +=
                    is_gc(seq, j) as f64 + is_gc(seq, j + 1) as f64 + is_gc(seq, j + 2) as f64;
                j -= 3;
            }
            gsize = ((*nod.offset(i as isize)).stop_val - (*nod.offset(i as isize)).ndx).abs()
                as f64
                + 3.0;
            (*nod.offset(i as isize)).gc_cont = gc[fr] / gsize;
            last[fr] = (*nod.offset(i as isize)).ndx;
        }
        i -= 1;
    }
    for k in 0..3 {
        gc[k] = 0.0;
    }
    for i in 0..nn {
        let fr = ((*nod.offset(i as isize)).ndx % 3) as usize;
        if (*nod.offset(i as isize)).strand == -1 && (*nod.offset(i as isize)).type_ == STOP {
            last[fr] = (*nod.offset(i as isize)).ndx;
            gc[fr] = is_gc(seq, (*nod.offset(i as isize)).ndx) as f64
                + is_gc(seq, (*nod.offset(i as isize)).ndx - 1) as f64
                + is_gc(seq, (*nod.offset(i as isize)).ndx - 2) as f64;
        } else if (*nod.offset(i as isize)).strand == -1 {
            let mut j = last[fr] + 3;
            while j <= (*nod.offset(i as isize)).ndx {
                gc[fr] +=
                    is_gc(seq, j) as f64 + is_gc(seq, j + 1) as f64 + is_gc(seq, j + 2) as f64;
                j += 3;
            }
            gsize = ((*nod.offset(i as isize)).stop_val - (*nod.offset(i as isize)).ndx).abs()
                as f64
                + 3.0;
            (*nod.offset(i as isize)).gc_cont = gc[fr] / gsize;
            last[fr] = (*nod.offset(i as isize)).ndx;
        }
    }
}

/*******************************************************************************
  Score each candidate's coding.
*******************************************************************************/

/// Score each candidate's coding. Also sharpens coding/noncoding thresholds
/// to prevent choosing interior starts when there is strong coding
/// continuing upstream. The routine makes three passes: an initial pass
/// summing the dicodon log-likelihoods from start to stop, a second pass
/// that penalizes start nodes with ascending coding to their left, and a
/// third pass that adds a length-based factor derived from the expected
/// no-stop probability.
pub unsafe fn raw_coding_score(
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nod: *mut Node,
    nn: c_int,
    tinf: *mut Training,
) {
    let mut last: [c_int; 3] = [0; 3];
    let mut score: [f64; 3] = [0.0; 3];
    let mut lfac: f64;
    let mut no_stop: f64;
    let mut gsize: f64;

    if (*tinf).trans_table != 11 {
        /* TGA or TAG is not a stop */
        no_stop = ((1.0 - (*tinf).gc) * (1.0 - (*tinf).gc) * (*tinf).gc) / 8.0;
        no_stop += ((1.0 - (*tinf).gc) * (1.0 - (*tinf).gc) * (1.0 - (*tinf).gc)) / 8.0;
        no_stop = 1.0 - no_stop;
    } else {
        no_stop = ((1.0 - (*tinf).gc) * (1.0 - (*tinf).gc) * (*tinf).gc) / 4.0;
        no_stop += ((1.0 - (*tinf).gc) * (1.0 - (*tinf).gc) * (1.0 - (*tinf).gc)) / 8.0;
        no_stop = 1.0 - no_stop;
    }

    /* Initial Pass: Score coding potential (start->stop) */
    for k in 0..3 {
        score[k] = 0.0;
    }
    let mut i = nn - 1;
    while i >= 0 {
        let fr = ((*nod.offset(i as isize)).ndx % 3) as usize;
        if (*nod.offset(i as isize)).strand == 1 && (*nod.offset(i as isize)).type_ == STOP {
            last[fr] = (*nod.offset(i as isize)).ndx;
            score[fr] = 0.0;
        } else if (*nod.offset(i as isize)).strand == 1 {
            let mut j = last[fr] - 3;
            while j >= (*nod.offset(i as isize)).ndx {
                score[fr] += (*tinf).gene_dc[mer_ndx(6, seq, j) as usize];
                j -= 3;
            }
            (*nod.offset(i as isize)).cscore = score[fr];
            last[fr] = (*nod.offset(i as isize)).ndx;
        }
        i -= 1;
    }
    for k in 0..3 {
        score[k] = 0.0;
    }
    for i in 0..nn {
        let fr = ((*nod.offset(i as isize)).ndx % 3) as usize;
        if (*nod.offset(i as isize)).strand == -1 && (*nod.offset(i as isize)).type_ == STOP {
            last[fr] = (*nod.offset(i as isize)).ndx;
            score[fr] = 0.0;
        } else if (*nod.offset(i as isize)).strand == -1 {
            let mut j = last[fr] + 3;
            while j <= (*nod.offset(i as isize)).ndx {
                score[fr] += (*tinf).gene_dc[mer_ndx(6, rseq, slen - j - 1) as usize];
                j += 3;
            }
            (*nod.offset(i as isize)).cscore = score[fr];
            last[fr] = (*nod.offset(i as isize)).ndx;
        }
    }

    /* Second Pass: Penalize start nodes with ascending coding to their left */
    for k in 0..3 {
        score[k] = -10000.0;
    }
    for i in 0..nn {
        let fr = ((*nod.offset(i as isize)).ndx % 3) as usize;
        if (*nod.offset(i as isize)).strand == 1 && (*nod.offset(i as isize)).type_ == STOP {
            score[fr] = -10000.0;
        } else if (*nod.offset(i as isize)).strand == 1 {
            if (*nod.offset(i as isize)).cscore > score[fr] {
                score[fr] = (*nod.offset(i as isize)).cscore;
            } else {
                (*nod.offset(i as isize)).cscore -= score[fr] - (*nod.offset(i as isize)).cscore;
            }
        }
    }
    for k in 0..3 {
        score[k] = -10000.0;
    }
    i = nn - 1;
    while i >= 0 {
        let fr = ((*nod.offset(i as isize)).ndx % 3) as usize;
        if (*nod.offset(i as isize)).strand == -1 && (*nod.offset(i as isize)).type_ == STOP {
            score[fr] = -10000.0;
        } else if (*nod.offset(i as isize)).strand == -1 {
            if (*nod.offset(i as isize)).cscore > score[fr] {
                score[fr] = (*nod.offset(i as isize)).cscore;
            } else {
                (*nod.offset(i as isize)).cscore -= score[fr] - (*nod.offset(i as isize)).cscore;
            }
        }
        i -= 1;
    }

    /* Third Pass: Add length-based factor to the score */
    let lfac_min = ((1.0 - no_stop.powf(80.0)) / no_stop.powf(80.0)).ln();
    let lfac_max = ((1.0 - no_stop.powf(1000.0)) / no_stop.powf(1000.0)).ln();
    for i in 0..nn {
        let fr = ((*nod.offset(i as isize)).ndx % 3) as usize;
        if (*nod.offset(i as isize)).strand == 1 && (*nod.offset(i as isize)).type_ == STOP {
            score[fr] = -10000.0;
        } else if (*nod.offset(i as isize)).strand == 1 {
            gsize = (((*nod.offset(i as isize)).stop_val - (*nod.offset(i as isize)).ndx).abs()
                as f64
                + 3.0)
                / 3.0;
            if gsize > 1000.0 {
                lfac = lfac_max;
                lfac -= lfac_min;
                lfac *= (gsize - 80.0) / 920.0;
            } else {
                let tmp = no_stop.powf(gsize);
                lfac = ((1.0 - tmp) / tmp).ln();
                lfac -= lfac_min;
            }
            if lfac > score[fr] {
                score[fr] = lfac;
            } else {
                lfac -= dmax(dmin(score[fr] - lfac, lfac), 0.0);
            }
            if lfac > 3.0 && (*nod.offset(i as isize)).cscore < 0.5 * lfac {
                (*nod.offset(i as isize)).cscore = 0.5 * lfac;
            }
            (*nod.offset(i as isize)).cscore += lfac;
        }
    }
    i = nn - 1;
    while i >= 0 {
        let fr = ((*nod.offset(i as isize)).ndx % 3) as usize;
        if (*nod.offset(i as isize)).strand == -1 && (*nod.offset(i as isize)).type_ == STOP {
            score[fr] = -10000.0;
        } else if (*nod.offset(i as isize)).strand == -1 {
            gsize = (((*nod.offset(i as isize)).stop_val - (*nod.offset(i as isize)).ndx).abs()
                as f64
                + 3.0)
                / 3.0;
            if gsize > 1000.0 {
                lfac = lfac_max;
                lfac -= lfac_min;
                lfac *= (gsize - 80.0) / 920.0;
            } else {
                let tmp = no_stop.powf(gsize);
                lfac = ((1.0 - tmp) / tmp).ln();
                lfac -= lfac_min;
            }
            if lfac > score[fr] {
                score[fr] = lfac;
            } else {
                lfac -= dmax(dmin(score[fr] - lfac, lfac), 0.0);
            }
            if lfac > 3.0 && (*nod.offset(i as isize)).cscore < 0.5 * lfac {
                (*nod.offset(i as isize)).cscore = 0.5 * lfac;
            }
            (*nod.offset(i as isize)).cscore += lfac;
        }
        i -= 1;
    }
}

/*******************************************************************************
  Examines the results of the SD motif search to determine if this organism
  uses an SD motif or not.
*******************************************************************************/

/// Examines the results of the SD motif search to determine if this organism
/// uses an SD motif or not. Some motif of 3-6bp has to be good or we set
/// `uses_sd` to 0, which will cause Prodigal to run the non-SD motif finder
/// for starts.
pub unsafe fn determine_sd_usage(tinf: *mut Training) {
    (*tinf).uses_sd = 1;
    if (*tinf).rbs_wt[0] >= 0.0 {
        (*tinf).uses_sd = 0;
    }
    if (*tinf).rbs_wt[16] < 1.0
        && (*tinf).rbs_wt[13] < 1.0
        && (*tinf).rbs_wt[15] < 1.0
        && ((*tinf).rbs_wt[0] >= -0.5
            || ((*tinf).rbs_wt[22] < 2.0 && (*tinf).rbs_wt[24] < 2.0 && (*tinf).rbs_wt[27] < 2.0))
    {
        (*tinf).uses_sd = 0;
    }
}

/*******************************************************************************
  RBS Scoring Function: Calculate the RBS motif and then multiply it by the
  appropriate weight for that motif.
*******************************************************************************/

/// RBS scoring function: calculate the RBS motif and then multiply it by the
/// appropriate weight for that motif (determined in the start training
/// function). For each non-edge start node, the upstream window is scanned
/// with both `shine_dalgarno_exact` and `shine_dalgarno_mm` and the best
/// scoring exact and mismatch indices are stored in `nod[i].rbs[0]` and
/// `nod[i].rbs[1]`.
pub const RBS_SLOTS: usize = 15;
pub const RBS_MASKS_PER_NODE: usize = RBS_SLOTS * 2;

/// The Shine-Dalgarno candidate set at a start depends only on the sequence; the model enters
/// solely through the `rbs_wt` argmax, so the candidates are recorded once and folded per model.
pub unsafe fn record_rbs_masks(
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nod: *mut Node,
    nn: c_int,
    masks: *mut u32,
) {
    for i in 0..nn {
        let slot_base = i as usize * RBS_MASKS_PER_NODE;
        for k in 0..RBS_MASKS_PER_NODE {
            *masks.add(slot_base + k) = 0;
        }
        if (*nod.offset(i as isize)).type_ == STOP || (*nod.offset(i as isize)).edge == 1 {
            continue;
        }
        let (wseq, first, start) = if (*nod.offset(i as isize)).strand == 1 {
            (
                seq,
                (*nod.offset(i as isize)).ndx - 20,
                (*nod.offset(i as isize)).ndx,
            )
        } else if (*nod.offset(i as isize)).strand == -1 {
            (
                rseq,
                slen - (*nod.offset(i as isize)).ndx - 21,
                slen - 1 - (*nod.offset(i as isize)).ndx,
            )
        } else {
            continue;
        };
        for slot in 0..RBS_SLOTS {
            let j = first + slot as c_int;
            if (*nod.offset(i as isize)).strand == 1 {
                if j < 0 {
                    continue;
                }
            } else if j > slen - 1 {
                continue;
            }
            *masks.add(slot_base + slot) = shine_dalgarno_exact_mask(wseq, j, start);
            *masks.add(slot_base + RBS_SLOTS + slot) = shine_dalgarno_mm_mask(wseq, j, start);
        }
    }
}

unsafe fn best_of_mask(mask: u32, rwt: *const f64) -> c_int {
    let mut best: c_int = 0;
    let mut rest = mask & !1u32;
    while rest != 0 {
        let value = rest.trailing_zeros() as c_int;
        rest &= rest - 1;
        let weight = *rwt.add(value as usize);
        let current = *rwt.add(best as usize);
        if weight > current || (weight == current && value > best) {
            best = value;
        }
    }
    best
}

unsafe fn rbs_score_from_masks(nod: *mut Node, nn: c_int, tinf: *mut Training, masks: *const u32) {
    let rwt = (*tinf).rbs_wt.as_ptr();
    for i in 0..nn {
        if (*nod.offset(i as isize)).type_ == STOP || (*nod.offset(i as isize)).edge == 1 {
            continue;
        }
        (*nod.offset(i as isize)).rbs[0] = 0;
        (*nod.offset(i as isize)).rbs[1] = 0;
        let slot_base = i as usize * RBS_MASKS_PER_NODE;
        for slot in 0..RBS_SLOTS {
            let exact = *masks.add(slot_base + slot);
            if exact != 0 {
                let value = best_of_mask(exact, rwt);
                if value > (*nod.offset(i as isize)).rbs[0] {
                    (*nod.offset(i as isize)).rbs[0] = value;
                }
            }
            let mm = *masks.add(slot_base + RBS_SLOTS + slot);
            if mm != 0 {
                let value = best_of_mask(mm, rwt);
                if value > (*nod.offset(i as isize)).rbs[1] {
                    (*nod.offset(i as isize)).rbs[1] = value;
                }
            }
        }
    }
}

pub unsafe fn rbs_score(
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nod: *mut Node,
    nn: c_int,
    tinf: *mut Training,
) {
    let mut cur_sc: [c_int; 2];

    /* Scan all starts looking for RBS's */
    for i in 0..nn {
        if (*nod.offset(i as isize)).type_ == STOP || (*nod.offset(i as isize)).edge == 1 {
            continue;
        }
        (*nod.offset(i as isize)).rbs[0] = 0;
        (*nod.offset(i as isize)).rbs[1] = 0;
        if (*nod.offset(i as isize)).strand == 1 {
            let mut j = (*nod.offset(i as isize)).ndx - 20;
            while j <= (*nod.offset(i as isize)).ndx - 6 {
                if j < 0 {
                    j += 1;
                    continue;
                }
                cur_sc = [0; 2];
                cur_sc[0] = shine_dalgarno_exact(
                    seq,
                    j,
                    (*nod.offset(i as isize)).ndx,
                    (*tinf).rbs_wt.as_ptr() as *mut f64,
                );
                cur_sc[1] = shine_dalgarno_mm(
                    seq,
                    j,
                    (*nod.offset(i as isize)).ndx,
                    (*tinf).rbs_wt.as_ptr() as *mut f64,
                );
                if cur_sc[0] > (*nod.offset(i as isize)).rbs[0] {
                    (*nod.offset(i as isize)).rbs[0] = cur_sc[0];
                }
                if cur_sc[1] > (*nod.offset(i as isize)).rbs[1] {
                    (*nod.offset(i as isize)).rbs[1] = cur_sc[1];
                }
                j += 1;
            }
        } else if (*nod.offset(i as isize)).strand == -1 {
            let mut j = slen - (*nod.offset(i as isize)).ndx - 21;
            while j <= slen - (*nod.offset(i as isize)).ndx - 7 {
                if j > slen - 1 {
                    j += 1;
                    continue;
                }
                cur_sc = [0; 2];
                cur_sc[0] = shine_dalgarno_exact(
                    rseq,
                    j,
                    slen - 1 - (*nod.offset(i as isize)).ndx,
                    (*tinf).rbs_wt.as_ptr() as *mut f64,
                );
                cur_sc[1] = shine_dalgarno_mm(
                    rseq,
                    j,
                    slen - 1 - (*nod.offset(i as isize)).ndx,
                    (*tinf).rbs_wt.as_ptr() as *mut f64,
                );
                if cur_sc[0] > (*nod.offset(i as isize)).rbs[0] {
                    (*nod.offset(i as isize)).rbs[0] = cur_sc[0];
                }
                if cur_sc[1] > (*nod.offset(i as isize)).rbs[1] {
                    (*nod.offset(i as isize)).rbs[1] = cur_sc[1];
                }
                j += 1;
            }
        }
    }
}

/*******************************************************************************
  Iterative Algorithm to train starts (Shine-Dalgarno motifs only).
*******************************************************************************/

/// Iterative algorithm to train starts. It begins with all the highest
/// coding starts in the model, scans for RBS/ATG-GTG-TTG usage, then starts
/// moving starts around attempting to match these discoveries. This start
/// trainer is for Shine-Dalgarno motifs only. Runs 10 iterations (a few more
/// than the typical 4-5 needed for convergence) and at the final iteration
/// also accumulates upstream base-composition counts, which are then
/// converted into per-position log scores in `tinf.ups_comp`.
pub unsafe fn train_starts_sd(
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nod: *mut Node,
    nn: c_int,
    tinf: *mut Training,
) {
    let mut fr: c_int;
    let mut rbs: [c_int; 3] = [0; 3];
    let mut type_: [c_int; 3] = [0; 3];
    let mut bndx: [c_int; 3] = [0; 3];
    let mut max_rb: c_int;
    let mut sum: f64;
    let wt: f64 = (*tinf).st_wt;
    let mut rbg: [f64; 28] = [0.0; 28];
    let mut rreal: [f64; 28] = [0.0; 28];
    let mut best: [f64; 3] = [0.0; 3];
    let mut sthresh: f64 = 35.0;
    let mut tbg: [f64; 3] = [0.0; 3];
    let mut treal: [f64; 3] = [0.0; 3];

    for j in 0..3 {
        (*tinf).type_wt[j] = 0.0;
    }
    for j in 0..28 {
        (*tinf).rbs_wt[j] = 0.0;
    }
    for i in 0..32 {
        for j in 0..4 {
            (*tinf).ups_comp[i][j] = 0.0;
        }
    }

    /* Build the background of random types */
    for i in 0..3 {
        tbg[i] = 0.0;
    }
    for i in 0..nn {
        if (*nod.offset(i as isize)).type_ == STOP {
            continue;
        }
        tbg[(*nod.offset(i as isize)).type_ as usize] += 1.0;
    }
    sum = 0.0;
    for i in 0..3 {
        sum += tbg[i];
    }
    for i in 0..3 {
        tbg[i] /= sum;
    }

    /* Iterate 10 times through the list of nodes */
    for i in 0..10 {
        /* Recalculate the RBS motif background */
        for j in 0..28 {
            rbg[j] = 0.0;
        }
        for j in 0..nn {
            if (*nod.offset(j as isize)).type_ == STOP || (*nod.offset(j as isize)).edge == 1 {
                continue;
            }
            if (*tinf).rbs_wt[(*nod.offset(j as isize)).rbs[0] as usize]
                > (*tinf).rbs_wt[(*nod.offset(j as isize)).rbs[1] as usize] + 1.0
                || (*nod.offset(j as isize)).rbs[1] == 0
            {
                max_rb = (*nod.offset(j as isize)).rbs[0];
            } else if (*tinf).rbs_wt[(*nod.offset(j as isize)).rbs[0] as usize]
                < (*tinf).rbs_wt[(*nod.offset(j as isize)).rbs[1] as usize] - 1.0
                || (*nod.offset(j as isize)).rbs[0] == 0
            {
                max_rb = (*nod.offset(j as isize)).rbs[1];
            } else {
                max_rb = dmax(
                    (*nod.offset(j as isize)).rbs[0] as f64,
                    (*nod.offset(j as isize)).rbs[1] as f64,
                ) as c_int;
            }
            rbg[max_rb as usize] += 1.0;
        }
        sum = 0.0;
        for j in 0..28 {
            sum += rbg[j];
        }
        for j in 0..28 {
            rbg[j] /= sum;
        }

        for j in 0..28 {
            rreal[j] = 0.0;
        }
        for j in 0..3 {
            treal[j] = 0.0;
        }

        /* Forward strand pass */
        for j in 0..3 {
            best[j] = 0.0;
            bndx[j] = -1;
            rbs[j] = 0;
            type_[j] = 0;
        }
        for j in 0..nn {
            if (*nod.offset(j as isize)).type_ != STOP && (*nod.offset(j as isize)).edge == 1 {
                continue;
            }
            fr = (*nod.offset(j as isize)).ndx % 3;
            if (*nod.offset(j as isize)).type_ == STOP && (*nod.offset(j as isize)).strand == 1 {
                if best[fr as usize] >= sthresh
                    && (*nod.offset(bndx[fr as usize] as isize)).ndx % 3 == fr
                {
                    rreal[rbs[fr as usize] as usize] += 1.0;
                    treal[type_[fr as usize] as usize] += 1.0;
                    if i == 9 {
                        count_upstream_composition(
                            seq,
                            slen,
                            1,
                            (*nod.offset(bndx[fr as usize] as isize)).ndx,
                            tinf,
                        );
                    }
                }
                best[fr as usize] = 0.0;
                bndx[fr as usize] = -1;
                rbs[fr as usize] = 0;
                type_[fr as usize] = 0;
            } else if (*nod.offset(j as isize)).strand == 1 {
                if (*tinf).rbs_wt[(*nod.offset(j as isize)).rbs[0] as usize]
                    > (*tinf).rbs_wt[(*nod.offset(j as isize)).rbs[1] as usize] + 1.0
                    || (*nod.offset(j as isize)).rbs[1] == 0
                {
                    max_rb = (*nod.offset(j as isize)).rbs[0];
                } else if (*tinf).rbs_wt[(*nod.offset(j as isize)).rbs[0] as usize]
                    < (*tinf).rbs_wt[(*nod.offset(j as isize)).rbs[1] as usize] - 1.0
                    || (*nod.offset(j as isize)).rbs[0] == 0
                {
                    max_rb = (*nod.offset(j as isize)).rbs[1];
                } else {
                    max_rb = dmax(
                        (*nod.offset(j as isize)).rbs[0] as f64,
                        (*nod.offset(j as isize)).rbs[1] as f64,
                    ) as c_int;
                }
                if (*nod.offset(j as isize)).cscore
                    + wt * (*tinf).rbs_wt[max_rb as usize]
                    + wt * (*tinf).type_wt[(*nod.offset(j as isize)).type_ as usize]
                    >= best[fr as usize]
                {
                    best[fr as usize] =
                        (*nod.offset(j as isize)).cscore + wt * (*tinf).rbs_wt[max_rb as usize];
                    best[fr as usize] +=
                        wt * (*tinf).type_wt[(*nod.offset(j as isize)).type_ as usize];
                    bndx[fr as usize] = j;
                    type_[fr as usize] = (*nod.offset(j as isize)).type_;
                    rbs[fr as usize] = max_rb;
                }
            }
        }

        /* Reverse strand pass */
        for j in 0..3 {
            best[j] = 0.0;
            bndx[j] = -1;
            rbs[j] = 0;
            type_[j] = 0;
        }
        let mut j = nn - 1;
        while j >= 0 {
            if (*nod.offset(j as isize)).type_ != STOP && (*nod.offset(j as isize)).edge == 1 {
                j -= 1;
                continue;
            }
            fr = (*nod.offset(j as isize)).ndx % 3;
            if (*nod.offset(j as isize)).type_ == STOP && (*nod.offset(j as isize)).strand == -1 {
                if best[fr as usize] >= sthresh
                    && (*nod.offset(bndx[fr as usize] as isize)).ndx % 3 == fr
                {
                    rreal[rbs[fr as usize] as usize] += 1.0;
                    treal[type_[fr as usize] as usize] += 1.0;
                    if i == 9 {
                        count_upstream_composition(
                            rseq,
                            slen,
                            -1,
                            (*nod.offset(bndx[fr as usize] as isize)).ndx,
                            tinf,
                        );
                    }
                }
                best[fr as usize] = 0.0;
                bndx[fr as usize] = -1;
                rbs[fr as usize] = 0;
                type_[fr as usize] = 0;
            } else if (*nod.offset(j as isize)).strand == -1 {
                if (*tinf).rbs_wt[(*nod.offset(j as isize)).rbs[0] as usize]
                    > (*tinf).rbs_wt[(*nod.offset(j as isize)).rbs[1] as usize] + 1.0
                    || (*nod.offset(j as isize)).rbs[1] == 0
                {
                    max_rb = (*nod.offset(j as isize)).rbs[0];
                } else if (*tinf).rbs_wt[(*nod.offset(j as isize)).rbs[0] as usize]
                    < (*tinf).rbs_wt[(*nod.offset(j as isize)).rbs[1] as usize] - 1.0
                    || (*nod.offset(j as isize)).rbs[0] == 0
                {
                    max_rb = (*nod.offset(j as isize)).rbs[1];
                } else {
                    max_rb = dmax(
                        (*nod.offset(j as isize)).rbs[0] as f64,
                        (*nod.offset(j as isize)).rbs[1] as f64,
                    ) as c_int;
                }
                if (*nod.offset(j as isize)).cscore
                    + wt * (*tinf).rbs_wt[max_rb as usize]
                    + wt * (*tinf).type_wt[(*nod.offset(j as isize)).type_ as usize]
                    >= best[fr as usize]
                {
                    best[fr as usize] =
                        (*nod.offset(j as isize)).cscore + wt * (*tinf).rbs_wt[max_rb as usize];
                    best[fr as usize] +=
                        wt * (*tinf).type_wt[(*nod.offset(j as isize)).type_ as usize];
                    bndx[fr as usize] = j;
                    type_[fr as usize] = (*nod.offset(j as isize)).type_;
                    rbs[fr as usize] = max_rb;
                }
            }
            j -= 1;
        }

        sum = 0.0;
        for j in 0..28 {
            sum += rreal[j];
        }
        if sum == 0.0 {
            for j in 0..28 {
                (*tinf).rbs_wt[j] = 0.0;
            }
        } else {
            for j in 0..28 {
                rreal[j] /= sum;
                if rbg[j] != 0.0 {
                    (*tinf).rbs_wt[j] = (rreal[j] / rbg[j]).ln();
                } else {
                    (*tinf).rbs_wt[j] = -4.0;
                }
                if (*tinf).rbs_wt[j] > 4.0 {
                    (*tinf).rbs_wt[j] = 4.0;
                }
                if (*tinf).rbs_wt[j] < -4.0 {
                    (*tinf).rbs_wt[j] = -4.0;
                }
            }
        }
        sum = 0.0;
        for j in 0..3 {
            sum += treal[j];
        }
        if sum == 0.0 {
            for j in 0..3 {
                (*tinf).type_wt[j] = 0.0;
            }
        } else {
            for j in 0..3 {
                treal[j] /= sum;
                if tbg[j] != 0.0 {
                    (*tinf).type_wt[j] = (treal[j] / tbg[j]).ln();
                } else {
                    (*tinf).type_wt[j] = -4.0;
                }
                if (*tinf).type_wt[j] > 4.0 {
                    (*tinf).type_wt[j] = 4.0;
                }
                if (*tinf).type_wt[j] < -4.0 {
                    (*tinf).type_wt[j] = -4.0;
                }
            }
        }
        if sum <= nn as f64 / 2000.0 {
            sthresh /= 2.0;
        }
    }

    /* Convert upstream base composition to a log score */
    for i in 0..32 {
        sum = 0.0;
        for j in 0..4 {
            sum += (*tinf).ups_comp[i][j];
        }
        if sum == 0.0 {
            for j in 0..4 {
                (*tinf).ups_comp[i][j] = 0.0;
            }
        } else {
            for j in 0..4 {
                (*tinf).ups_comp[i][j] /= sum;
                if (*tinf).gc > 0.1 && (*tinf).gc < 0.9 {
                    if j == 0 || j == 3 {
                        (*tinf).ups_comp[i][j] =
                            ((*tinf).ups_comp[i][j] * 2.0 / (1.0 - (*tinf).gc)).ln();
                    } else {
                        (*tinf).ups_comp[i][j] = ((*tinf).ups_comp[i][j] * 2.0 / (*tinf).gc).ln();
                    }
                } else if (*tinf).gc <= 0.1 {
                    if j == 0 || j == 3 {
                        (*tinf).ups_comp[i][j] = ((*tinf).ups_comp[i][j] * 2.0 / 0.90).ln();
                    } else {
                        (*tinf).ups_comp[i][j] = ((*tinf).ups_comp[i][j] * 2.0 / 0.10).ln();
                    }
                } else {
                    if j == 0 || j == 3 {
                        (*tinf).ups_comp[i][j] = ((*tinf).ups_comp[i][j] * 2.0 / 0.10).ln();
                    } else {
                        (*tinf).ups_comp[i][j] = ((*tinf).ups_comp[i][j] * 2.0 / 0.90).ln();
                    }
                }
                if (*tinf).ups_comp[i][j] > 4.0 {
                    (*tinf).ups_comp[i][j] = 4.0;
                }
                if (*tinf).ups_comp[i][j] < -4.0 {
                    (*tinf).ups_comp[i][j] = -4.0;
                }
            }
        }
    }
}

/*******************************************************************************
  Iterative Algorithm to train starts (non-SD version).
*******************************************************************************/

/// Iterative algorithm to train starts. It begins with all the highest
/// coding starts in the model, scans for RBS/ATG-GTG-TTG usage, then starts
/// moving starts around attempting to match these discoveries. Unlike the SD
/// algorithm, it allows for any popular motif to be discovered. Runs 20
/// iterations across three motif-finding stages (0: count all motifs, 1:
/// count best motif and its sub-motifs, 2: count only the best motif), and
/// at the final iteration also collects upstream base composition.
pub unsafe fn train_starts_nonsd(
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nod: *mut Node,
    nn: c_int,
    tinf: *mut Training,
) {
    let mut fr: c_int;
    let mut bndx: [c_int; 3] = [0; 3];
    let mut mgood: [[[c_int; 4096]; 4]; 4] = [[[0; 4096]; 4]; 4];
    let mut stage: c_int;
    let mut sum: f64;
    let mut ngenes: f64;
    let wt: f64 = (*tinf).st_wt;
    let mut best: [f64; 3] = [0.0; 3];
    let mut sthresh: f64 = 35.0;
    let mut tbg: [f64; 3] = [0.0; 3];
    let mut treal: [f64; 3] = [0.0; 3];
    let mut mbg: [[[f64; 4096]; 4]; 4] = [[[0.0; 4096]; 4]; 4];
    let mut mreal: [[[f64; 4096]; 4]; 4] = [[[0.0; 4096]; 4]; 4];
    let mut zbg: f64;
    let mut zreal: f64;

    for i in 0..32 {
        for j in 0..4 {
            (*tinf).ups_comp[i][j] = 0.0;
        }
    }

    /* Build the background of random types */
    for i in 0..3 {
        (*tinf).type_wt[i] = 0.0;
    }
    for i in 0..3 {
        tbg[i] = 0.0;
    }
    for i in 0..nn {
        if (*nod.offset(i as isize)).type_ == STOP {
            continue;
        }
        tbg[(*nod.offset(i as isize)).type_ as usize] += 1.0;
    }
    sum = 0.0;
    for i in 0..3 {
        sum += tbg[i];
    }
    for i in 0..3 {
        tbg[i] /= sum;
    }

    /* Iterate 20 times through the list of nodes */
    for i in 0..20 {
        /* Determine which stage of motif finding we're in */
        if i < 4 {
            stage = 0;
        } else if i < 12 {
            stage = 1;
        } else {
            stage = 2;
        }

        /* Recalculate the upstream motif background and set 'real' counts to 0 */
        for j in 0..4 {
            for k in 0..4 {
                for l in 0..4096 {
                    mbg[j][k][l] = 0.0;
                }
            }
        }
        zbg = 0.0;
        for j in 0..nn {
            if (*nod.offset(j as isize)).type_ == STOP || (*nod.offset(j as isize)).edge == 1 {
                continue;
            }
            find_best_upstream_motif(tinf, seq, rseq, slen, &mut *nod.offset(j as isize), stage);
            update_motif_counts(
                mbg.as_mut_ptr() as *mut [[f64; 4096]; 4],
                &mut zbg,
                seq,
                rseq,
                slen,
                &mut *nod.offset(j as isize),
                stage,
            );
        }
        sum = 0.0;
        for j in 0..4 {
            for k in 0..4 {
                for l in 0..4096 {
                    sum += mbg[j][k][l];
                }
            }
        }
        sum += zbg;
        for j in 0..4 {
            for k in 0..4 {
                for l in 0..4096 {
                    mbg[j][k][l] /= sum;
                }
            }
        }
        zbg /= sum;

        /* Reset counts of 'real' motifs/types to 0 */
        for j in 0..4 {
            for k in 0..4 {
                for l in 0..4096 {
                    mreal[j][k][l] = 0.0;
                }
            }
        }
        zreal = 0.0;
        for j in 0..3 {
            treal[j] = 0.0;
        }
        ngenes = 0.0;

        /* Forward strand pass */
        for j in 0..3 {
            best[j] = 0.0;
            bndx[j] = -1;
        }
        for j in 0..nn {
            if (*nod.offset(j as isize)).type_ != STOP && (*nod.offset(j as isize)).edge == 1 {
                continue;
            }
            fr = (*nod.offset(j as isize)).ndx % 3;
            if (*nod.offset(j as isize)).type_ == STOP && (*nod.offset(j as isize)).strand == 1 {
                if best[fr as usize] >= sthresh {
                    ngenes += 1.0;
                    treal[(*nod.offset(bndx[fr as usize] as isize)).type_ as usize] += 1.0;
                    update_motif_counts(
                        mreal.as_mut_ptr() as *mut [[f64; 4096]; 4],
                        &mut zreal,
                        seq,
                        rseq,
                        slen,
                        &mut *nod.offset(bndx[fr as usize] as isize),
                        stage,
                    );
                    if i == 19 {
                        count_upstream_composition(
                            seq,
                            slen,
                            1,
                            (*nod.offset(bndx[fr as usize] as isize)).ndx,
                            tinf,
                        );
                    }
                }
                best[fr as usize] = 0.0;
                bndx[fr as usize] = -1;
            } else if (*nod.offset(j as isize)).strand == 1 {
                if (*nod.offset(j as isize)).cscore
                    + wt * (*nod.offset(j as isize)).mot.score
                    + wt * (*tinf).type_wt[(*nod.offset(j as isize)).type_ as usize]
                    >= best[fr as usize]
                {
                    best[fr as usize] =
                        (*nod.offset(j as isize)).cscore + wt * (*nod.offset(j as isize)).mot.score;
                    best[fr as usize] +=
                        wt * (*tinf).type_wt[(*nod.offset(j as isize)).type_ as usize];
                    bndx[fr as usize] = j;
                }
            }
        }

        /* Reverse strand pass */
        for j in 0..3 {
            best[j] = 0.0;
            bndx[j] = -1;
        }
        let mut j = nn - 1;
        while j >= 0 {
            if (*nod.offset(j as isize)).type_ != STOP && (*nod.offset(j as isize)).edge == 1 {
                j -= 1;
                continue;
            }
            fr = (*nod.offset(j as isize)).ndx % 3;
            if (*nod.offset(j as isize)).type_ == STOP && (*nod.offset(j as isize)).strand == -1 {
                if best[fr as usize] >= sthresh {
                    ngenes += 1.0;
                    treal[(*nod.offset(bndx[fr as usize] as isize)).type_ as usize] += 1.0;
                    update_motif_counts(
                        mreal.as_mut_ptr() as *mut [[f64; 4096]; 4],
                        &mut zreal,
                        seq,
                        rseq,
                        slen,
                        &mut *nod.offset(bndx[fr as usize] as isize),
                        stage,
                    );
                    if i == 19 {
                        count_upstream_composition(
                            rseq,
                            slen,
                            -1,
                            (*nod.offset(bndx[fr as usize] as isize)).ndx,
                            tinf,
                        );
                    }
                }
                best[fr as usize] = 0.0;
                bndx[fr as usize] = -1;
            } else if (*nod.offset(j as isize)).strand == -1 {
                if (*nod.offset(j as isize)).cscore
                    + wt * (*nod.offset(j as isize)).mot.score
                    + wt * (*tinf).type_wt[(*nod.offset(j as isize)).type_ as usize]
                    >= best[fr as usize]
                {
                    best[fr as usize] =
                        (*nod.offset(j as isize)).cscore + wt * (*nod.offset(j as isize)).mot.score;
                    best[fr as usize] +=
                        wt * (*tinf).type_wt[(*nod.offset(j as isize)).type_ as usize];
                    bndx[fr as usize] = j;
                }
            }
            j -= 1;
        }

        /* Update the log likelihood weights for type and RBS motifs */
        if stage < 2 {
            build_coverage_map(
                mreal.as_mut_ptr() as *mut [[f64; 4096]; 4],
                mgood.as_mut_ptr() as *mut [[c_int; 4096]; 4],
                ngenes,
                stage,
            );
        }
        sum = 0.0;
        for j in 0..4 {
            for k in 0..4 {
                for l in 0..4096 {
                    sum += mreal[j][k][l];
                }
            }
        }
        sum += zreal;
        if sum == 0.0 {
            for j in 0..4 {
                for k in 0..4 {
                    for l in 0..4096 {
                        (*tinf).mot_wt[j][k][l] = 0.0;
                    }
                }
            }
            (*tinf).no_mot = 0.0;
        } else {
            for j in 0..4usize {
                for k in 0..4usize {
                    for l in 0..4096usize {
                        if mgood[j][k][l] == 0 {
                            zreal += mreal[j][k][l];
                            zbg += mreal[j][k][l];
                            mreal[j][k][l] = 0.0;
                            mbg[j][k][l] = 0.0;
                        }
                        mreal[j][k][l] /= sum;
                        if mbg[j][k][l] != 0.0 {
                            (*tinf).mot_wt[j][k][l] = (mreal[j][k][l] / mbg[j][k][l]).ln();
                        } else {
                            (*tinf).mot_wt[j][k][l] = -4.0;
                        }
                        if (*tinf).mot_wt[j][k][l] > 4.0 {
                            (*tinf).mot_wt[j][k][l] = 4.0;
                        }
                        if (*tinf).mot_wt[j][k][l] < -4.0 {
                            (*tinf).mot_wt[j][k][l] = -4.0;
                        }
                    }
                }
            }
        }
        zreal /= sum;
        if zbg != 0.0 {
            (*tinf).no_mot = (zreal / zbg).ln();
        } else {
            (*tinf).no_mot = -4.0;
        }
        if (*tinf).no_mot > 4.0 {
            (*tinf).no_mot = 4.0;
        }
        if (*tinf).no_mot < -4.0 {
            (*tinf).no_mot = -4.0;
        }
        sum = 0.0;
        for j in 0..3 {
            sum += treal[j];
        }
        if sum == 0.0 {
            for j in 0..3 {
                (*tinf).type_wt[j] = 0.0;
            }
        } else {
            for j in 0..3 {
                treal[j] /= sum;
                if tbg[j] != 0.0 {
                    (*tinf).type_wt[j] = (treal[j] / tbg[j]).ln();
                } else {
                    (*tinf).type_wt[j] = -4.0;
                }
                if (*tinf).type_wt[j] > 4.0 {
                    (*tinf).type_wt[j] = 4.0;
                }
                if (*tinf).type_wt[j] < -4.0 {
                    (*tinf).type_wt[j] = -4.0;
                }
            }
        }
        if sum <= nn as f64 / 2000.0 {
            sthresh /= 2.0;
        }
    }

    /* Convert upstream base composition to a log score */
    for i in 0..32 {
        sum = 0.0;
        for j in 0..4 {
            sum += (*tinf).ups_comp[i][j];
        }
        if sum == 0.0 {
            for j in 0..4 {
                (*tinf).ups_comp[i][j] = 0.0;
            }
        } else {
            for j in 0..4 {
                (*tinf).ups_comp[i][j] /= sum;
                if (*tinf).gc > 0.1 && (*tinf).gc < 0.9 {
                    if j == 0 || j == 3 {
                        (*tinf).ups_comp[i][j] =
                            ((*tinf).ups_comp[i][j] * 2.0 / (1.0 - (*tinf).gc)).ln();
                    } else {
                        (*tinf).ups_comp[i][j] = ((*tinf).ups_comp[i][j] * 2.0 / (*tinf).gc).ln();
                    }
                } else if (*tinf).gc <= 0.1 {
                    if j == 0 || j == 3 {
                        (*tinf).ups_comp[i][j] = ((*tinf).ups_comp[i][j] * 2.0 / 0.90).ln();
                    } else {
                        (*tinf).ups_comp[i][j] = ((*tinf).ups_comp[i][j] * 2.0 / 0.10).ln();
                    }
                } else {
                    if j == 0 || j == 3 {
                        (*tinf).ups_comp[i][j] = ((*tinf).ups_comp[i][j] * 2.0 / 0.10).ln();
                    } else {
                        (*tinf).ups_comp[i][j] = ((*tinf).ups_comp[i][j] * 2.0 / 0.90).ln();
                    }
                }
                if (*tinf).ups_comp[i][j] > 4.0 {
                    (*tinf).ups_comp[i][j] = 4.0;
                }
                if (*tinf).ups_comp[i][j] < -4.0 {
                    (*tinf).ups_comp[i][j] = -4.0;
                }
            }
        }
    }
}

/*******************************************************************************
  For a given start, record the base composition of the upstream region.
*******************************************************************************/

/// For a given start, record the base composition of the upstream region at
/// positions -1 and -2 and -15 to -44. This will be used to supplement the
/// SD (or other) motif finder with additional information.
pub unsafe fn count_upstream_composition(
    seq: *mut u8,
    slen: c_int,
    strand: c_int,
    pos: c_int,
    tinf: *mut Training,
) {
    let mut count: usize = 0;
    let start: c_int;
    if strand == 1 {
        start = pos;
    } else {
        start = slen - 1 - pos;
    }

    for i in 1..45 {
        if i > 2 && i < 15 {
            continue;
        }
        if start - i >= 0 {
            (*tinf).ups_comp[count][mer_ndx(1, seq, start - i) as usize] += 1.0;
        }
        count += 1;
    }
}

/*******************************************************************************
  For a given start, score the base composition of the upstream region.
*******************************************************************************/

/// For a given start, score the base composition of the upstream region at
/// positions -1 and -2 and -15 to -44. This will be used to supplement the
/// SD (or other) motif finder with additional information.
pub unsafe fn score_upstream_composition(
    seq: *mut u8,
    slen: c_int,
    nod: *mut Node,
    tinf: *mut Training,
) {
    let mut count: usize = 0;
    let start: c_int;
    if (*nod).strand == 1 {
        start = (*nod).ndx;
    } else {
        start = slen - 1 - (*nod).ndx;
    }

    (*nod).uscore = 0.0;
    for i in 1..45 {
        if i > 2 && i < 15 {
            continue;
        }
        if start - i < 0 {
            continue;
        }
        (*nod).uscore +=
            0.4 * (*tinf).st_wt * (*tinf).ups_comp[count][mer_ndx(1, seq, start - i) as usize];
        count += 1;
    }
}

/*******************************************************************************
  Given the weights for various motifs/distances from the training file,
  return the highest scoring mer/spacer combination.
*******************************************************************************/

/// Given the weights for various motifs/distances from the training file,
/// return the highest scoring mer/spacer combination of 3-6bp motifs with a
/// spacer ranging from 3bp to 15bp. In the final stage of start training
/// (`stage == 2`), only good scoring motifs are returned; otherwise a "no
/// motif" placeholder is stored in `nod.mot`.
pub unsafe fn find_best_upstream_motif(
    tinf: *mut Training,
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nod: *mut Node,
    stage: c_int,
) {
    let start: c_int;
    let mut spacer: c_int;
    let mut spacendx: c_int;
    let mut index: c_int;
    let mut max_spacer: c_int = 0;
    let mut max_spacendx: c_int = 0;
    let mut max_len: c_int = 0;
    let mut max_ndx: c_int = 0;
    let mut max_sc: f64 = -100.0;
    let mut score: f64;
    let wseq: *const u8;

    if (*nod).type_ == STOP || (*nod).edge == 1 {
        return;
    }
    if (*nod).strand == 1 {
        wseq = seq;
        start = (*nod).ndx;
    } else {
        wseq = rseq;
        start = slen - 1 - (*nod).ndx;
    }

    let mut i: c_int = 3;
    while i >= 0 {
        let mut j = start - 18 - i;
        while j <= start - 6 - i {
            if j < 0 {
                j += 1;
                continue;
            }
            spacer = start - j - i - 3;
            if j <= start - 16 - i {
                spacendx = 3;
            } else if j <= start - 14 - i {
                spacendx = 2;
            } else if j >= start - 7 - i {
                spacendx = 1;
            } else {
                spacendx = 0;
            }
            index = mer_ndx(i + 3, wseq as *mut u8, j);
            score = (*tinf).mot_wt[i as usize][spacendx as usize][index as usize];
            if score > max_sc {
                max_sc = score;
                max_spacendx = spacendx;
                max_spacer = spacer;
                max_ndx = index;
                max_len = i + 3;
            }
            j += 1;
        }
        i -= 1;
    }

    if stage == 2 && (max_sc == -4.0 || max_sc < (*tinf).no_mot + 0.69) {
        (*nod).mot.ndx = 0;
        (*nod).mot.len = 0;
        (*nod).mot.spacendx = 0;
        (*nod).mot.spacer = 0;
        (*nod).mot.score = (*tinf).no_mot;
    } else {
        (*nod).mot.ndx = max_ndx;
        (*nod).mot.len = max_len;
        (*nod).mot.spacendx = max_spacendx;
        (*nod).mot.spacer = max_spacer;
        (*nod).mot.score = max_sc;
    }
}

/*******************************************************************************
  Update the motif counts from a putative "real" start.
*******************************************************************************/

/// Update the motif counts from a putative "real" start. This is done in
/// three stages. In stage 0, all motifs sizes 3-6bp in the region with
/// spacer 3-15bp are counted. In stage 1, only the best motif and all its
/// subsets are counted (e.g. for AGGAG, we would count AGGAG, AGGA, GGAG,
/// AGG, GGA, and GAG). In stage 2, only the best single motif is counted.
#[allow(unused_assignments)]
pub unsafe fn update_motif_counts(
    mcnt: *mut [[f64; 4096]; 4],
    zero: *mut f64,
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nod: *mut Node,
    stage: c_int,
) {
    let start: c_int;
    let mut spacendx: c_int;
    let wseq: *const u8;
    let mot: *const Motif = &(*nod).mot;

    if (*nod).type_ == STOP || (*nod).edge == 1 {
        return;
    }
    if (*mot).len == 0 {
        *zero += 1.0;
        return;
    }

    if (*nod).strand == 1 {
        wseq = seq;
        start = (*nod).ndx;
    } else {
        wseq = rseq;
        start = slen - 1 - (*nod).ndx;
    }

    /* Stage 0: Count all motifs. */
    if stage == 0 {
        let mut i: c_int = 3;
        while i >= 0 {
            let mut j = start - 18 - i;
            while j <= start - 6 - i {
                if j < 0 {
                    j += 1;
                    continue;
                }
                if j <= start - 16 - i {
                    spacendx = 3;
                } else if j <= start - 14 - i {
                    spacendx = 2;
                } else if j >= start - 7 - i {
                    spacendx = 1;
                } else {
                    spacendx = 0;
                }
                let ndx = mer_ndx(i + 3, wseq as *mut u8, j);
                for k in 0..4 {
                    (*mcnt.offset(i as isize))[k][ndx as usize] += 1.0;
                }
                j += 1;
            }
            i -= 1;
        }
    }
    /* Stage 1: Count only the best motif, but also count all its sub-motifs. */
    else if stage == 1 {
        (*mcnt.offset(((*mot).len - 3) as isize))[(*mot).spacendx as usize][(*mot).ndx as usize] +=
            1.0;
        let mut i: c_int = 0;
        while i < (*mot).len - 3 {
            let mut j = start - (*mot).spacer - (*mot).len;
            while j <= start - (*mot).spacer - (i + 3) {
                if j < 0 {
                    j += 1;
                    continue;
                }
                if j <= start - 16 - i {
                    spacendx = 3;
                } else if j <= start - 14 - i {
                    spacendx = 2;
                } else if j >= start - 7 - i {
                    spacendx = 1;
                } else {
                    spacendx = 0;
                }
                let ndx = mer_ndx(i + 3, wseq as *mut u8, j);
                (*mcnt.offset(i as isize))[spacendx as usize][ndx as usize] += 1.0;
                j += 1;
            }
            i += 1;
        }
    }
    /* Stage 2: Only count the highest scoring motif. */
    else if stage == 2 {
        (*mcnt.offset(((*mot).len - 3) as isize))[(*mot).spacendx as usize][(*mot).ndx as usize] +=
            1.0;
    }
}

/*******************************************************************************
  Build coverage map for motifs.
*******************************************************************************/

/// In addition to log likelihood, we also require a motif to actually be
/// present a good portion of the time in an absolute sense across the
/// genome. The coverage map is just a numerical map of whether or not to
/// accept a putative motif as a real one (despite its log likelihood score).
/// A motif is considered "good" if it contains a 3-base subset of itself
/// that is present in at least 20% of the total genes. In the final stage
/// of iterative start training, all motifs are labeled good. 0 = bad,
/// 1 = good, 2 = good w/mismatch.
pub unsafe fn build_coverage_map(
    real: *mut [[f64; 4096]; 4],
    good: *mut [[c_int; 4096]; 4],
    ng: f64,
    _stage: c_int,
) {
    let mut decomp: [c_int; 3] = [0; 3];
    let thresh: f64 = 0.2;

    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4096 {
                (*good.offset(i as isize))[j][k] = 0;
            }
        }
    }

    /* 3-base motifs */
    for i in 0..4 {
        for j in 0..64 {
            if (*real.offset(0))[i][j] / ng >= thresh {
                for k in 0..4 {
                    (*good.offset(0))[k][j] = 1;
                }
            }
        }
    }

    /* 4-base motifs, must contain two valid 3-base motifs */
    for i in 0..4 {
        for j in 0..256 {
            decomp[0] = ((j & 252) >> 2) as c_int;
            decomp[1] = (j & 63) as c_int;
            if (*good.offset(0))[i][decomp[0] as usize] == 0
                || (*good.offset(0))[i][decomp[1] as usize] == 0
            {
                continue;
            }
            (*good.offset(1))[i][j] = 1;
        }
    }

    /* 5-base motifs */
    for i in 0..4 {
        for j in 0..1024 {
            decomp[0] = ((j & 1008) >> 4) as c_int;
            decomp[1] = ((j & 252) >> 2) as c_int;
            decomp[2] = (j & 63) as c_int;
            if (*good.offset(0))[i][decomp[0] as usize] == 0
                || (*good.offset(0))[i][decomp[1] as usize] == 0
                || (*good.offset(0))[i][decomp[2] as usize] == 0
            {
                continue;
            }
            (*good.offset(2))[i][j] = 1;
            let mut tmp = j as c_int;
            let mut k: c_int = 0;
            while k <= 16 {
                tmp = tmp ^ k;
                let mut l: c_int = 0;
                while l <= 32 {
                    tmp = tmp ^ l;
                    if (*good.offset(2))[i][tmp as usize] == 0 {
                        (*good.offset(2))[i][tmp as usize] = 2;
                    }
                    l += 32;
                }
                k += 16;
            }
        }
    }

    /* 6-base motifs, must contain two valid 5-base motifs */
    for i in 0..4 {
        for j in 0..4096 {
            decomp[0] = ((j & 4092) >> 2) as c_int;
            decomp[1] = (j & 1023) as c_int;
            if (*good.offset(2))[i][decomp[0] as usize] == 0
                || (*good.offset(2))[i][decomp[1] as usize] == 0
            {
                continue;
            }
            if (*good.offset(2))[i][decomp[0] as usize] == 1
                && (*good.offset(2))[i][decomp[1] as usize] == 1
            {
                (*good.offset(3))[i][j] = 1;
            } else {
                (*good.offset(3))[i][j] = 2;
            }
        }
    }
}

/*******************************************************************************
  Intergenic modifier for connecting two genes.
*******************************************************************************/

/// When connecting two genes, we add a bonus for the -1 and -4 base overlaps
/// on the same strand, which often signify an operon and negate the need for
/// an RBS for the second gene. In addition, we add a slight bonus when genes
/// are close and a slight penalty when switching strands or having a large
/// intergenic space.
#[inline(always)]
pub unsafe fn intergenic_mod(n1: *mut Node, n2: *mut Node, tinf: *mut Training) -> f64 {
    let mut rval: f64 = 0.0;
    let mut ovlp: f64 = 0.0;

    if ((*n1).strand == 1
        && (*n2).strand == 1
        && ((*n1).ndx + 2 == (*n2).ndx || (*n1).ndx - 1 == (*n2).ndx))
        || ((*n1).strand == -1
            && (*n2).strand == -1
            && ((*n1).ndx + 2 == (*n2).ndx || (*n1).ndx - 1 == (*n2).ndx))
    {
        if (*n1).strand == 1 && (*n2).rscore < 0.0 {
            rval -= (*n2).rscore;
        }
        if (*n1).strand == -1 && (*n1).rscore < 0.0 {
            rval -= (*n1).rscore;
        }
        if (*n1).strand == 1 && (*n2).uscore < 0.0 {
            rval -= (*n2).uscore;
        }
        if (*n1).strand == -1 && (*n1).uscore < 0.0 {
            rval -= (*n1).uscore;
        }
    }
    let dist = ((*n1).ndx - (*n2).ndx).abs();
    if (*n1).strand == 1 && (*n2).strand == 1 && (*n1).ndx + 2 >= (*n2).ndx {
        ovlp = 1.0;
    } else if (*n1).strand == -1 && (*n2).strand == -1 && (*n1).ndx >= (*n2).ndx + 2 {
        ovlp = 1.0;
    }
    if dist > 3 * OPER_DIST || (*n1).strand != (*n2).strand {
        rval -= 0.15 * (*tinf).st_wt;
    } else if (dist <= OPER_DIST && ovlp == 0.0) || (dist as f64) < 0.25 * OPER_DIST as f64 {
        rval += (2.0 - dist as f64 / OPER_DIST as f64) * 0.15 * (*tinf).st_wt;
    }
    rval
}

/*******************************************************************************
  Write detailed scoring information about every single possible gene.
*******************************************************************************/

/// Write detailed scoring information about every single possible gene.
/// Only done at the user's request. Emits a tab-separated table per
/// candidate start (one block per stop codon) including total/coding/start
/// scores, RBS motif description, spacer, upstream/type scores and GC
/// content; sorts nodes by stop first, then restores the original ordering
/// before returning.
pub unsafe fn write_start_file(
    fh: c_int,
    nod: *mut Node,
    nn: c_int,
    tinf: *mut Training,
    sctr: c_int,
    slen: c_int,
    is_meta: c_int,
    mdesc: *mut c_char,
    version: *mut c_char,
    header: *mut c_char,
) {
    let mut prev_stop: c_int = -1;
    let mut prev_strand: c_int = 0;
    let mut st_type: c_int;
    let mut rbs1: f64;
    let mut rbs2: f64;

    static SD_STRING: [&str; 28] = [
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
    static SD_SPACER: [&str; 28] = [
        "None", "3-4bp", "13-15bp", "13-15bp", "11-12bp", "3-4bp", "11-12bp", "11-12bp", "3-4bp",
        "5-10bp", "13-15bp", "3-4bp", "11-12bp", "5-10bp", "5-10bp", "5-10bp", "5-10bp", "11-12bp",
        "3-4bp", "5-10bp", "11-12bp", "3-4bp", "5-10bp", "3-4bp", "5-10bp", "11-12bp", "3-4bp",
        "5-10bp",
    ];
    static TYPE_STRING: [&str; 4] = ["ATG", "GTG", "TTG", "Edge"];

    let mut qt: [c_char; 10] = [0; 10];

    let header_str = cstr(header);
    let version_str = cstr(version);

    /* Initialize sequence data */
    let seq_data = format!("seqnum={};seqlen={};seqhdr=\"{}\"", sctr, slen, header_str);

    /* Initialize run data string */
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

    {
        let nodes_slice = std::slice::from_raw_parts_mut(nod, nn as usize);
        nodes_slice.sort_unstable_by(|a, b| {
            a.stop_val
                .cmp(&b.stop_val)
                .then(b.strand.cmp(&a.strand))
                .then(a.ndx.cmp(&b.ndx))
        });
    }

    fprint(fh, &format!("# Sequence Data: {}\n", seq_data));
    fprint(fh, &format!("# Run Data: {}\n\n", run_data));

    fprint(fh, "Beg\tEnd\tStd\tTotal\tCodPot\tStrtSc\tCodon\tRBSMot\t");
    fprint(fh, "Spacer\tRBSScr\tUpsScr\tTypeScr\tGCCont\n");

    for i in 0..nn {
        let ni = &*nod.offset(i as isize);
        if ni.type_ == STOP {
            continue;
        }
        if ni.edge == 1 {
            st_type = 3;
        } else {
            st_type = ni.type_;
        }
        if ni.stop_val != prev_stop || ni.strand != prev_strand {
            prev_stop = ni.stop_val;
            prev_strand = ni.strand;
            fprint(fh, "\n");
        }
        if ni.strand == 1 {
            fprint(
                fh,
                &format!(
                    "{}\t{}\t+\t{:.2}\t{:.2}\t{:.2}\t{}\t",
                    ni.ndx + 1,
                    ni.stop_val + 3,
                    ni.cscore + ni.sscore,
                    ni.cscore,
                    ni.sscore,
                    TYPE_STRING[st_type as usize],
                ),
            );
        }
        if ni.strand == -1 {
            fprint(
                fh,
                &format!(
                    "{}\t{}\t-\t{:.2}\t{:.2}\t{:.2}\t{}\t",
                    ni.stop_val - 1,
                    ni.ndx + 1,
                    ni.cscore + ni.sscore,
                    ni.cscore,
                    ni.sscore,
                    TYPE_STRING[st_type as usize],
                ),
            );
        }
        rbs1 = (*tinf).rbs_wt[ni.rbs[0] as usize] * (*tinf).st_wt;
        rbs2 = (*tinf).rbs_wt[ni.rbs[1] as usize] * (*tinf).st_wt;
        if (*tinf).uses_sd == 1 {
            if rbs1 > rbs2 {
                fprint(
                    fh,
                    &format!(
                        "{}\t{}\t{:.2}\t",
                        SD_STRING[ni.rbs[0] as usize], SD_SPACER[ni.rbs[0] as usize], ni.rscore,
                    ),
                );
            } else {
                fprint(
                    fh,
                    &format!(
                        "{}\t{}\t{:.2}\t",
                        SD_STRING[ni.rbs[1] as usize], SD_SPACER[ni.rbs[1] as usize], ni.rscore,
                    ),
                );
            }
        } else {
            mer_text(qt.as_mut_ptr(), ni.mot.len, ni.mot.ndx);
            let qt_str = cstr(qt.as_ptr());
            if (*tinf).no_mot > -0.5 && rbs1 > rbs2 && rbs1 > ni.mot.score * (*tinf).st_wt {
                fprint(
                    fh,
                    &format!(
                        "{}\t{}\t{:.2}\t",
                        SD_STRING[ni.rbs[0] as usize], SD_SPACER[ni.rbs[0] as usize], ni.rscore,
                    ),
                );
            } else if (*tinf).no_mot > -0.5 && rbs2 >= rbs1 && rbs2 > ni.mot.score * (*tinf).st_wt {
                fprint(
                    fh,
                    &format!(
                        "{}\t{}\t{:.2}\t",
                        SD_STRING[ni.rbs[1] as usize], SD_SPACER[ni.rbs[1] as usize], ni.rscore,
                    ),
                );
            } else {
                if ni.mot.len == 0 {
                    fprint(fh, &format!("None\tNone\t{:.2}\t", ni.rscore));
                } else {
                    fprint(
                        fh,
                        &format!("{}\t{}bp\t{:.2}\t", qt_str, ni.mot.spacer, ni.rscore,),
                    );
                }
            }
        }
        fprint(
            fh,
            &format!("{:.2}\t{:.2}\t{:.3}\n", ni.uscore, ni.tscore, ni.gc_cont,),
        );
    }
    fprint(fh, "\n");

    {
        let nodes_slice = std::slice::from_raw_parts_mut(nod, nn as usize);
        nodes_slice.sort_unstable_by(|a, b| a.ndx.cmp(&b.ndx).then(b.strand.cmp(&a.strand)));
    }
}

/* Checks to see if a node boundary crosses a mask */

/// Checks to see if a node boundary crosses a mask. Returns 1 if the
/// interval `[x, y]` overlaps any region in `mlist`, otherwise 0.
pub unsafe fn cross_mask(x: c_int, y: c_int, mlist: *mut Mask, nm: c_int) -> c_int {
    for i in 0..nm {
        if y < (*mlist.offset(i as isize)).begin || x > (*mlist.offset(i as isize)).end {
            continue;
        }
        return 1;
    }
    0
}

/* Return the minimum of two numbers */

/// Return the minimum of two numbers.
pub unsafe fn dmin(x: f64, y: f64) -> f64 {
    if x < y {
        x
    } else {
        y
    }
}

/* Return the maximum of two numbers */

/// Return the maximum of two numbers.
pub unsafe fn dmax(x: f64, y: f64) -> f64 {
    if x > y {
        x
    } else {
        y
    }
}

/* Sorting routine for nodes */

/// qsort-style comparator that orders nodes by ascending `ndx`, breaking
/// ties by descending `strand` (forward strand first).
pub unsafe fn compare_nodes(v1: *const c_void, v2: *const c_void) -> c_int {
    let n1 = v1 as *const Node;
    let n2 = v2 as *const Node;
    if (*n1).ndx < (*n2).ndx {
        return -1;
    }
    if (*n1).ndx > (*n2).ndx {
        return 1;
    }
    if (*n1).strand > (*n2).strand {
        return -1;
    }
    if (*n1).strand < (*n2).strand {
        return 1;
    }
    0
}

/* Sorts all nodes by common stop */

/// qsort-style comparator that sorts all nodes by common stop. Orders by
/// ascending `stop_val`, then descending `strand`, then ascending `ndx`.
pub unsafe fn stopcmp_nodes(v1: *const c_void, v2: *const c_void) -> c_int {
    let n1 = v1 as *const Node;
    let n2 = v2 as *const Node;
    if (*n1).stop_val < (*n2).stop_val {
        return -1;
    }
    if (*n1).stop_val > (*n2).stop_val {
        return 1;
    }
    if (*n1).strand > (*n2).strand {
        return -1;
    }
    if (*n1).strand < (*n2).strand {
        return 1;
    }
    if (*n1).ndx < (*n2).ndx {
        return -1;
    }
    if (*n1).ndx > (*n2).ndx {
        return 1;
    }
    0
}

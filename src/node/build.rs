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

use super::*;

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
    tinf: *const Training,
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
    tinf: *const Training,
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

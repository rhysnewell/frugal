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
    tinf: *const Training,
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
    tinf: *const Training,
    closed: c_int,
    is_meta: c_int,
    rbs_masks: *const u32,
) {
    /* Step 1: Calculate raw coding potential for every start-stop pair. */
    calc_orf_gc(seq, rseq, slen, nod, nn, tinf);
    raw_coding_score(seq, rseq, slen, nod, nn, tinf);
    let mer = |forward: bool, len: c_int, at: c_int| match forward {
        true => mer_ndx(len, seq, at),
        false => mer_ndx(len, rseq, at),
    };
    score_starts(
        seq, rseq, slen, nod, nn, tinf, closed, is_meta, rbs_masks, mer,
    );
}

#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(crate) unsafe fn score_starts(
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nod: *mut Node,
    nn: c_int,
    tinf: *const Training,
    closed: c_int,
    is_meta: c_int,
    rbs_masks: *const u32,
    mer: impl Fn(bool, c_int, c_int) -> c_int,
) {
    let mut negf: f64;
    let mut posf: f64;
    let mut rbs1: f64;
    let mut rbs2: f64;
    let mut sd_score: f64;
    let mut edge_gene: f64;
    let mut min_meta_len: f64;

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
            let node = &mut *nod.offset(i as isize);
            let forward = node.strand == 1;
            find_best_upstream_motif_with(tinf, slen, node, 2, |len, at| mer(forward, len, at));
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
            let forward = (*nod.offset(i as isize)).strand == 1;
            score_upstream_composition_with(slen, &mut *nod.offset(i as isize), tinf, |at| {
                mer(forward, 1, at)
            });

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

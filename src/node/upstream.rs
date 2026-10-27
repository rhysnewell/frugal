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
    tinf: *const Training,
) {
    score_upstream_composition_with(slen, nod, tinf, |at| mer_ndx(1, seq, at))
}

#[inline(always)]
pub(crate) unsafe fn score_upstream_composition_with(
    slen: c_int,
    nod: *mut Node,
    tinf: *const Training,
    base: impl Fn(c_int) -> c_int,
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
        (*nod).uscore += 0.4 * (*tinf).st_wt * (*tinf).ups_comp[count][base(start - i) as usize];
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
    tinf: *const Training,
    seq: *mut u8,
    rseq: *mut u8,
    slen: c_int,
    nod: *mut Node,
    stage: c_int,
) {
    let wseq = match (*nod).strand == 1 {
        true => seq,
        false => rseq,
    };
    find_best_upstream_motif_with(tinf, slen, nod, stage, |len, at| mer_ndx(len, wseq, at))
}

#[inline(always)]
pub(crate) unsafe fn find_best_upstream_motif_with(
    tinf: *const Training,
    slen: c_int,
    nod: *mut Node,
    stage: c_int,
    mer: impl Fn(c_int, c_int) -> c_int,
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

    if (*nod).type_ == STOP || (*nod).edge == 1 {
        return;
    }
    if (*nod).strand == 1 {
        start = (*nod).ndx;
    } else {
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
            index = mer(i + 3, j);
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

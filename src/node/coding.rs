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
    _tinf: *const Training,
    _seq: *mut u8,
    _rseq: *mut u8,
    _slen: c_int,
    _nod: *mut Node,
    _nn: c_int,
) {
    /* Not implemented in original Prodigal */
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
    _tinf: *const Training,
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
    tinf: *const Training,
) {
    raw_coding_score_with(
        nod,
        nn,
        tinf,
        |j| mer_ndx(6, seq, j) as usize,
        |j| mer_ndx(6, rseq, slen - j - 1) as usize,
    )
}

// A shorter mer's index is the low bits of the hexamer's at the same position, so one table per
// strand serves every model's coding, motif and upstream reads.
pub(crate) unsafe fn mer_table(seq: *mut u8, slen: c_int, out: &mut Vec<u16>) {
    out.clear();
    out.extend((0..slen).map(|at| mer_ndx(6, seq, at) as u16));
}

#[inline(always)]
pub(crate) fn mer_from(table: &[u16], len: c_int, at: c_int) -> c_int {
    c_int::from(table[at as usize] & ((1u16 << (2 * len)) - 1))
}

#[inline(always)]
pub(crate) unsafe fn raw_coding_score_with(
    nod: *mut Node,
    nn: c_int,
    tinf: *const Training,
    forward: impl Fn(c_int) -> usize,
    reverse: impl Fn(c_int) -> usize,
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
                score[fr] += (*tinf).gene_dc[forward(j)];
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
                score[fr] += (*tinf).gene_dc[reverse(j)];
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

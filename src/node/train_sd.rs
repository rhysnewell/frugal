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

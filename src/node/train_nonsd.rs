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

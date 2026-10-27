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

pub(super) unsafe fn rbs_score_from_masks(
    nod: *mut Node,
    nn: c_int,
    tinf: *const Training,
    masks: *const u32,
) {
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
    tinf: *const Training,
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

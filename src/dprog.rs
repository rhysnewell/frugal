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

use crate::types::{Node, Training, MAX_NODE_DIST, STOP};
use std::os::raw::c_int;

use crate::connection::{backward_start, backward_stop, forward_start, forward_stop};
use crate::connection_filter::{allowed_for, class_of};
use crate::node::intergenic_mod;


/// `flag` 0 scores on the GC frame plot alone, to build a training set; 1 scores on coding and
/// RBS for the final call.
pub unsafe fn dprog(nod: *mut Node, nn: c_int, tinf: *const Training, flag: c_int) -> c_int {
    let mut min: c_int;
    let mut max_ndx: c_int = -1;
    let mut max_sc: f64 = -1.0;
    let mut path: c_int;
    let mut nxt: c_int;
    let mut tmp: c_int;

    if nn == 0 {
        return -1;
    }
    for i in 0..nn {
        (*nod.offset(i as isize)).score = 0.0;
        (*nod.offset(i as isize)).traceb = -1;
        (*nod.offset(i as isize)).tracef = -1;
    }

    // Only classes[j] for j below the current i is ever read, and iteration j writes it, so
    // nothing here is observed before the loop sets it.
    let mut classes = vec![0u8; nn as usize];

    for i in 0..nn {
        /* Set up distance constraints for making connections, */
        /* but make exceptions for giant ORFs.                 */
        if i < MAX_NODE_DIST {
            min = 0;
        } else {
            min = i - MAX_NODE_DIST;
        }
        if (*nod.offset(i as isize)).strand == -1
            && (*nod.offset(i as isize)).type_ != STOP
            && (*nod.offset(min as isize)).ndx >= (*nod.offset(i as isize)).stop_val
        {
            while min >= 0 && (*nod.offset(i as isize)).ndx != (*nod.offset(i as isize)).stop_val {
                min -= 1;
            }
        }
        if (*nod.offset(i as isize)).strand == 1
            && (*nod.offset(i as isize)).type_ == STOP
            && (*nod.offset(min as isize)).ndx >= (*nod.offset(i as isize)).stop_val
        {
            while min >= 0 && (*nod.offset(i as isize)).ndx != (*nod.offset(i as isize)).stop_val {
                min -= 1;
            }
        }
        if min < MAX_NODE_DIST {
            min = 0;
        } else {
            min -= MAX_NODE_DIST;
        }
        let n2: *mut Node = nod.offset(i as isize);
        let table = allowed_for(n2);
        macro_rules! pass {
            ($scorer:path) => {
                for j in min..i {
                    if *table.get_unchecked(*classes.get_unchecked(j as usize) as usize) != 0 {
                        $scorer(nod, j, n2, tinf, flag);
                    }
                }
            };
        }
        match 2 * u8::from((*n2).strand != 1) + u8::from((*n2).type_ == STOP) {
            0 => pass!(forward_start),
            1 => pass!(forward_stop),
            2 => pass!(backward_start),
            _ => pass!(backward_stop),
        }
        classes[i as usize] = class_of(n2);
    }
    for i in (0..nn).rev() {
        if (*nod.offset(i as isize)).strand == 1 && (*nod.offset(i as isize)).type_ != STOP {
            continue;
        }
        if (*nod.offset(i as isize)).strand == -1 && (*nod.offset(i as isize)).type_ == STOP {
            continue;
        }
        if (*nod.offset(i as isize)).score > max_sc {
            max_sc = (*nod.offset(i as isize)).score;
            max_ndx = i;
        }
    }

    /* First Pass: untangle the triple overlaps */
    path = max_ndx;
    while (*nod.offset(path as isize)).traceb != -1 {
        nxt = (*nod.offset(path as isize)).traceb;
        if (*nod.offset(path as isize)).strand == -1
            && (*nod.offset(path as isize)).type_ == STOP
            && (*nod.offset(nxt as isize)).strand == 1
            && (*nod.offset(nxt as isize)).type_ == STOP
            && (*nod.offset(path as isize)).ov_mark != -1
            && (*nod.offset(path as isize)).ndx > (*nod.offset(nxt as isize)).ndx
        {
            tmp = (*nod.offset(path as isize)).star_ptr
                [(*nod.offset(path as isize)).ov_mark as usize];
            let mut ii = tmp;
            while (*nod.offset(ii as isize)).ndx != (*nod.offset(tmp as isize)).stop_val {
                ii -= 1;
            }
            (*nod.offset(path as isize)).traceb = tmp;
            (*nod.offset(tmp as isize)).traceb = ii;
            (*nod.offset(ii as isize)).ov_mark = -1;
            (*nod.offset(ii as isize)).traceb = nxt;
        }
        path = (*nod.offset(path as isize)).traceb;
    }

    /* Second Pass: Untangle the simple overlaps */
    path = max_ndx;
    while (*nod.offset(path as isize)).traceb != -1 {
        nxt = (*nod.offset(path as isize)).traceb;
        if (*nod.offset(path as isize)).strand == -1
            && (*nod.offset(path as isize)).type_ != STOP
            && (*nod.offset(nxt as isize)).strand == 1
            && (*nod.offset(nxt as isize)).type_ == STOP
        {
            let mut ii = path;
            while (*nod.offset(ii as isize)).ndx != (*nod.offset(path as isize)).stop_val {
                ii -= 1;
            }
            if ii < 0 {
                path = nxt;
                continue;
            }
            (*nod.offset(path as isize)).traceb = ii;
            (*nod.offset(ii as isize)).traceb = nxt;
        }
        if (*nod.offset(path as isize)).strand == 1
            && (*nod.offset(path as isize)).type_ == STOP
            && (*nod.offset(nxt as isize)).strand == 1
            && (*nod.offset(nxt as isize)).type_ == STOP
        {
            (*nod.offset(path as isize)).traceb = (*nod.offset(nxt as isize)).star_ptr
                [((*nod.offset(path as isize)).ndx % 3) as usize];
            (*nod.offset((*nod.offset(path as isize)).traceb as isize)).traceb = nxt;
        }
        if (*nod.offset(path as isize)).strand == -1
            && (*nod.offset(path as isize)).type_ == STOP
            && (*nod.offset(nxt as isize)).strand == -1
            && (*nod.offset(nxt as isize)).type_ == STOP
        {
            (*nod.offset(path as isize)).traceb = (*nod.offset(path as isize)).star_ptr
                [((*nod.offset(nxt as isize)).ndx % 3) as usize];
            (*nod.offset((*nod.offset(path as isize)).traceb as isize)).traceb = nxt;
        }
        path = (*nod.offset(path as isize)).traceb;
    }

    /* Mark forward pointers */
    path = max_ndx;
    while (*nod.offset(path as isize)).traceb != -1 {
        (*nod.offset((*nod.offset(path as isize)).traceb as isize)).tracef = path;
        path = (*nod.offset(path as isize)).traceb;
    }

    if (*nod.offset(max_ndx as isize)).traceb == -1 {
        return -1;
    } else {
        return max_ndx;
    }
}

/// Sometimes bad genes creep into the model due to the node distance constraint
/// in the dynamic programming routine.  This routine just does a sweep through
/// the genes and eliminates ones with negative scores.
pub unsafe fn eliminate_bad_genes(nod: *mut Node, dbeg: c_int, tinf: *const Training) {
    let mut path: c_int;

    if dbeg == -1 {
        return;
    }
    path = dbeg;
    while (*nod.offset(path as isize)).traceb != -1 {
        path = (*nod.offset(path as isize)).traceb;
    }
    while (*nod.offset(path as isize)).tracef != -1 {
        if (*nod.offset(path as isize)).strand == 1 && (*nod.offset(path as isize)).type_ == STOP {
            let tracef = (*nod.offset(path as isize)).tracef;
            (*nod.offset(tracef as isize)).sscore += intergenic_mod(
                &mut *nod.offset(path as isize),
                &mut *nod.offset(tracef as isize),
                tinf,
            );
        }
        if (*nod.offset(path as isize)).strand == -1 && (*nod.offset(path as isize)).type_ != STOP {
            let tracef = (*nod.offset(path as isize)).tracef;
            (*nod.offset(path as isize)).sscore += intergenic_mod(
                &mut *nod.offset(path as isize),
                &mut *nod.offset(tracef as isize),
                tinf,
            );
        }
        path = (*nod.offset(path as isize)).tracef;
    }

    path = dbeg;
    while (*nod.offset(path as isize)).traceb != -1 {
        path = (*nod.offset(path as isize)).traceb;
    }
    while (*nod.offset(path as isize)).tracef != -1 {
        if (*nod.offset(path as isize)).strand == 1
            && (*nod.offset(path as isize)).type_ != STOP
            && (*nod.offset(path as isize)).cscore + (*nod.offset(path as isize)).sscore < 0.0
        {
            (*nod.offset(path as isize)).elim = 1;
            let tracef = (*nod.offset(path as isize)).tracef;
            (*nod.offset(tracef as isize)).elim = 1;
        }
        if (*nod.offset(path as isize)).strand == -1 && (*nod.offset(path as isize)).type_ == STOP {
            let tracef = (*nod.offset(path as isize)).tracef;
            if (*nod.offset(tracef as isize)).cscore + (*nod.offset(tracef as isize)).sscore < 0.0 {
                (*nod.offset(path as isize)).elim = 1;
                (*nod.offset(tracef as isize)).elim = 1;
            }
        }
        path = (*nod.offset(path as isize)).tracef;
    }
}

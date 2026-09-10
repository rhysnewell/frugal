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

use crate::types::{Node, Training, MAX_NODE_DIST, MAX_OPP_OVLP, STOP};
use std::os::raw::c_int;

use crate::connection_filter::{allowed_for, class_of, keep};
use crate::node::intergenic_mod;

/// Basic dynamic programming routine for predicting genes.  The `flag` variable
/// is set to 0 for the initial dynamic programming routine based solely on GC
/// frame plot (used to construct a training set).  If the flag is set to 1, the
/// routine does the final dynamic programming based on coding, RBS scores, etc.
pub unsafe fn dprog(nod: *mut Node, nn: c_int, tinf: *mut Training, flag: c_int) -> c_int {
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

    let mut classes = Vec::with_capacity(nn as usize);
    for i in 0..nn {
        classes.push(class_of(nod.offset(i as isize)));
    }
    let mut reachable: Vec<c_int> = Vec::with_capacity(MAX_NODE_DIST as usize * 2);

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
        let table = allowed_for(nod.offset(i as isize));
        keep(
            &classes[min as usize..i as usize],
            &table,
            &mut reachable,
        );
        for offset in &reachable {
            score_connection(nod, min + *offset, i, tinf, flag);
        }
        classes[i as usize] = class_of(nod.offset(i as isize));
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

/// This routine scores the connection between two nodes, the most basic of which
/// is 5'fwd->3'fwd (gene) and 3'rev->5'rev (rev gene).  If the connection ending
/// at `n2` is the maximal scoring model, it updates the pointers in the dynamic
/// programming model.  `n3` is used to handle overlaps, i.e. cases where 5->3'
/// overlaps 5'->3' on the same strand.  In this case, 3' connects directly to 3',
/// and `n3` is used to untangle the 5' end of the second gene.
#[inline(always)]
pub unsafe fn score_connection(
    nod: *mut Node,
    p1: c_int,
    p2: c_int,
    tinf: *mut Training,
    flag: c_int,
) {
    let n1: *mut Node = &mut *nod.offset(p1 as isize);
    let n2: *mut Node = &mut *nod.offset(p2 as isize);
    let mut n3: *mut Node;
    let mut left: c_int = (*n1).ndx;
    let mut right: c_int = (*n2).ndx;
    let bnd: c_int;
    let mut ovlp: c_int = 0;
    let mut maxfr: c_int = -1;
    let mut score: f64 = 0.0;
    let mut scr_mod: f64 = 0.0;
    let mut maxval: f64;

    /***********************/
    /* Invalid Connections */
    /***********************/

    /* 5'fwd->5'fwd, 5'rev->5'rev */
    if (*n1).type_ != STOP && (*n2).type_ != STOP && (*n1).strand == (*n2).strand {
        return;
    }
    /* 5'fwd->5'rev, 5'fwd->3'rev */
    else if (*n1).strand == 1 && (*n1).type_ != STOP && (*n2).strand == -1 {
        return;
    }
    /* 3'rev->5'fwd, 3'rev->3'fwd) */
    else if (*n1).strand == -1 && (*n1).type_ == STOP && (*n2).strand == 1 {
        return;
    }
    /* 5'rev->3'fwd */
    else if (*n1).strand == -1 && (*n1).type_ != STOP && (*n2).strand == 1 && (*n2).type_ == STOP
    {
        return;
    }

    /******************/
    /* Edge Artifacts */
    /******************/
    if (*n1).traceb == -1 && (*n1).strand == 1 && (*n1).type_ == STOP {
        return;
    }
    if (*n1).traceb == -1 && (*n1).strand == -1 && (*n1).type_ != STOP {
        return;
    }
    /*********/
    /* Genes */
    /*********/

    /* 5'fwd->3'fwd */
    else if (*n1).strand == (*n2).strand
        && (*n1).strand == 1
        && (*n1).type_ != STOP
        && (*n2).type_ == STOP
    {
        if (*n2).stop_val >= (*n1).ndx {
            return;
        }
        if (*n1).ndx % 3 != (*n2).ndx % 3 {
            return;
        }
        right += 2;
        if flag == 0 {
            scr_mod = (*tinf).bias[0] * (*n1).gc_score[0]
                + (*tinf).bias[1] * (*n1).gc_score[1]
                + (*tinf).bias[2] * (*n1).gc_score[2];
        } else if flag == 1 {
            score = (*n1).cscore + (*n1).sscore;
        }
    }
    /* 3'rev->5'rev */
    else if (*n1).strand == (*n2).strand
        && (*n1).strand == -1
        && (*n1).type_ == STOP
        && (*n2).type_ != STOP
    {
        if (*n1).stop_val <= (*n2).ndx {
            return;
        }
        if (*n1).ndx % 3 != (*n2).ndx % 3 {
            return;
        }
        left -= 2;
        if flag == 0 {
            scr_mod = (*tinf).bias[0] * (*n2).gc_score[0]
                + (*tinf).bias[1] * (*n2).gc_score[1]
                + (*tinf).bias[2] * (*n2).gc_score[2];
        } else if flag == 1 {
            score = (*n2).cscore + (*n2).sscore;
        }
    }
    /********************************/
    /* Intergenic Space (Noncoding) */
    /********************************/

    /* 3'fwd->5'fwd */
    else if (*n1).strand == 1 && (*n1).type_ == STOP && (*n2).strand == 1 && (*n2).type_ != STOP {
        left += 2;
        if left >= right {
            return;
        }
        if flag == 1 {
            score = intergenic_mod(n1, n2, tinf);
        }
    }
    /* 3'fwd->3'rev */
    else if (*n1).strand == 1 && (*n1).type_ == STOP && (*n2).strand == -1 && (*n2).type_ == STOP
    {
        left += 2;
        right -= 2;
        if left >= right {
            return;
        }
        /* Overlapping Gene Case 2: Three consecutive overlapping genes f r r */
        maxfr = -1;
        maxval = 0.0;
        for i in 0..3 {
            if (*n2).star_ptr[i as usize] == -1 {
                continue;
            }
            n3 = &mut *nod.offset((*n2).star_ptr[i as usize] as isize);
            ovlp = left - (*n3).stop_val + 3;
            if ovlp <= 0 || ovlp >= MAX_OPP_OVLP {
                continue;
            }
            if ovlp >= (*n3).ndx - left {
                continue;
            }
            if (*n1).traceb == -1 {
                continue;
            }
            if ovlp >= (*n3).stop_val - (*nod.offset((*n1).traceb as isize)).ndx - 2 {
                continue;
            }
            if (flag == 1 && (*n3).cscore + (*n3).sscore + intergenic_mod(n3, n2, tinf) > maxval)
                || (flag == 0
                    && (*tinf).bias[0] * (*n3).gc_score[0]
                        + (*tinf).bias[1] * (*n3).gc_score[1]
                        + (*tinf).bias[2] * (*n3).gc_score[2]
                        > maxval)
            {
                maxfr = i;
                maxval = (*n3).cscore + (*n3).sscore + intergenic_mod(n3, n2, tinf);
            }
        }
        if maxfr != -1 {
            n3 = &mut *nod.offset((*n2).star_ptr[maxfr as usize] as isize);
            if flag == 0 {
                scr_mod = (*tinf).bias[0] * (*n3).gc_score[0]
                    + (*tinf).bias[1] * (*n3).gc_score[1]
                    + (*tinf).bias[2] * (*n3).gc_score[2];
            } else if flag == 1 {
                score = (*n3).cscore + (*n3).sscore + intergenic_mod(n3, n2, tinf);
            }
        } else if flag == 1 {
            score = intergenic_mod(n1, n2, tinf);
        }
    }
    /* 5'rev->3'rev */
    else if (*n1).strand == -1 && (*n1).type_ != STOP && (*n2).strand == -1 && (*n2).type_ == STOP
    {
        right -= 2;
        if left >= right {
            return;
        }
        if flag == 1 {
            score = intergenic_mod(n1, n2, tinf);
        }
    }
    /* 5'rev->5'fwd */
    else if (*n1).strand == -1 && (*n1).type_ != STOP && (*n2).strand == 1 && (*n2).type_ != STOP
    {
        if left >= right {
            return;
        }
        if flag == 1 {
            score = intergenic_mod(n1, n2, tinf);
        }
    }
    /********************/
    /* Possible Operons */
    /********************/

    /* 3'fwd->3'fwd, check for a start just to left of first 3' */
    else if (*n1).strand == 1 && (*n2).strand == 1 && (*n1).type_ == STOP && (*n2).type_ == STOP {
        if (*n2).stop_val >= (*n1).ndx {
            return;
        }
        if (*n1).star_ptr[((*n2).ndx % 3) as usize] == -1 {
            return;
        }
        n3 = &mut *nod.offset((*n1).star_ptr[((*n2).ndx % 3) as usize] as isize);
        left = (*n3).ndx;
        right += 2;
        if flag == 0 {
            scr_mod = (*tinf).bias[0] * (*n3).gc_score[0]
                + (*tinf).bias[1] * (*n3).gc_score[1]
                + (*tinf).bias[2] * (*n3).gc_score[2];
        } else if flag == 1 {
            score = (*n3).cscore + (*n3).sscore + intergenic_mod(n1, n3, tinf);
        }
    }
    /* 3'rev->3'rev, check for a start just to right of second 3' */
    else if (*n1).strand == -1 && (*n1).type_ == STOP && (*n2).strand == -1 && (*n2).type_ == STOP
    {
        if (*n1).stop_val <= (*n2).ndx {
            return;
        }
        if (*n2).star_ptr[((*n1).ndx % 3) as usize] == -1 {
            return;
        }
        n3 = &mut *nod.offset((*n2).star_ptr[((*n1).ndx % 3) as usize] as isize);
        left -= 2;
        right = (*n3).ndx;
        if flag == 0 {
            scr_mod = (*tinf).bias[0] * (*n3).gc_score[0]
                + (*tinf).bias[1] * (*n3).gc_score[1]
                + (*tinf).bias[2] * (*n3).gc_score[2];
        } else if flag == 1 {
            score = (*n3).cscore + (*n3).sscore + intergenic_mod(n3, n2, tinf);
        }
    }
    /***************************************/
    /* Overlapping Opposite Strand 3' Ends */
    /***************************************/

    /* 3'for->5'rev */
    else if (*n1).strand == 1 && (*n1).type_ == STOP && (*n2).strand == -1 && (*n2).type_ != STOP
    {
        if (*n2).stop_val - 2 >= (*n1).ndx + 2 {
            return;
        }
        ovlp = ((*n1).ndx + 2) - ((*n2).stop_val - 2) + 1;
        if ovlp >= MAX_OPP_OVLP {
            return;
        }
        if ((*n1).ndx + 2 - (*n2).stop_val - 2 + 1) >= ((*n2).ndx - (*n1).ndx + 3 + 1) {
            return;
        }
        if (*n1).traceb == -1 {
            bnd = 0;
        } else {
            bnd = (*nod.offset((*n1).traceb as isize)).ndx;
        }
        if ((*n1).ndx + 2 - (*n2).stop_val - 2 + 1) >= ((*n2).stop_val - 3 - bnd + 1) {
            return;
        }
        left = (*n2).stop_val - 2;
        if flag == 0 {
            scr_mod = (*tinf).bias[0] * (*n2).gc_score[0]
                + (*tinf).bias[1] * (*n2).gc_score[1]
                + (*tinf).bias[2] * (*n2).gc_score[2];
        } else if flag == 1 {
            score = (*n2).cscore + (*n2).sscore - 0.15 * (*tinf).st_wt;
        }
    }

    if flag == 0 {
        score = ((right - left + 1 - (ovlp * 2)) as f64) * scr_mod;
    }

    if (*n1).score + score >= (*n2).score {
        (*n2).score = (*n1).score + score;
        (*n2).traceb = p1;
        (*n2).ov_mark = maxfr;
    }
}

/// Sometimes bad genes creep into the model due to the node distance constraint
/// in the dynamic programming routine.  This routine just does a sweep through
/// the genes and eliminates ones with negative scores.
pub unsafe fn eliminate_bad_genes(nod: *mut Node, dbeg: c_int, tinf: *mut Training) {
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

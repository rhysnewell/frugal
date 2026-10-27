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
    ATG, EDGE_BONUS, EDGE_UPS, GTG, MAX_SAM_OVLP, META_PEN, MIN_EDGE_GENE, MIN_GENE, Mask, Motif,
    Node, OPER_DIST, STOP, TTG, Training,
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

mod build;
mod coding;
mod output;
mod rbs;
mod score;
mod train_nonsd;
mod train_sd;
mod upstream;

pub use build::*;
pub use coding::*;
pub(crate) use coding::{mer_from, mer_table, raw_coding_score_with};
pub use output::*;
pub use rbs::*;
pub(crate) use score::score_starts;
pub use score::*;
pub use train_nonsd::*;
pub use train_sd::*;
pub use upstream::*;
pub(crate) use upstream::{find_best_upstream_motif_with, score_upstream_composition_with};

/*******************************************************************************
  Intergenic modifier for connecting two genes.
*******************************************************************************/

/// When connecting two genes, we add a bonus for the -1 and -4 base overlaps
/// on the same strand, which often signify an operon and negate the need for
/// an RBS for the second gene. In addition, we add a slight bonus when genes
/// are close and a slight penalty when switching strands or having a large
/// intergenic space.
#[inline(always)]
pub unsafe fn intergenic_mod(n1: *mut Node, n2: *mut Node, tinf: *const Training) -> f64 {
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
    if x < y { x } else { y }
}

/* Return the maximum of two numbers */

/// Return the maximum of two numbers.
pub unsafe fn dmax(x: f64, y: f64) -> f64 {
    if x > y { x } else { y }
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

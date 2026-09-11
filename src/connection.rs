// The right node holds still for a whole pass, so splitting the scorer on its strand and type
// deletes every branch the prefilter has already ruled out.

use crate::types::{Node, Training, MAX_OPP_OVLP, OPER_DIST, STOP};
use std::os::raw::c_int;

#[inline(always)]
unsafe fn intergenic_mod_same(n1: *mut Node, n2: *mut Node, tinf: *const Training) -> f64 {
    let mut rval: f64 = 0.0;

    if (*n1).ndx + 2 == (*n2).ndx || (*n1).ndx - 1 == (*n2).ndx {
        if (*n1).strand == 1 {
            if (*n2).rscore < 0.0 {
                rval -= (*n2).rscore;
            }
            if (*n2).uscore < 0.0 {
                rval -= (*n2).uscore;
            }
        } else {
            if (*n1).rscore < 0.0 {
                rval -= (*n1).rscore;
            }
            if (*n1).uscore < 0.0 {
                rval -= (*n1).uscore;
            }
        }
    }
    let dist = ((*n1).ndx - (*n2).ndx).abs();
    let ovlp = if (*n1).strand == 1 {
        (*n1).ndx + 2 >= (*n2).ndx
    } else {
        (*n1).ndx >= (*n2).ndx + 2
    };
    if dist > 3 * OPER_DIST {
        rval -= 0.15 * (*tinf).st_wt;
    } else if (dist <= OPER_DIST && !ovlp) || (dist as f64) < 0.25 * OPER_DIST as f64 {
        rval += (2.0 - dist as f64 / OPER_DIST as f64) * 0.15 * (*tinf).st_wt;
    }
    rval
}

#[inline(always)]
unsafe fn intergenic_mod_diff(tinf: *const Training) -> f64 {
    -0.15 * (*tinf).st_wt
}

#[inline(always)]
unsafe fn frame_bias(node: *mut Node, tinf: *const Training) -> f64 {
    (*tinf).bias[0] * (*node).gc_score[0]
        + (*tinf).bias[1] * (*node).gc_score[1]
        + (*tinf).bias[2] * (*node).gc_score[2]
}

#[inline(always)]
unsafe fn commit(n1: *mut Node, n2: *mut Node, p1: c_int, score: f64, maxfr: c_int) {
    if (*n1).score + score >= (*n2).score {
        (*n2).score = (*n1).score + score;
        (*n2).traceb = p1;
        (*n2).ov_mark = maxfr;
    }
}

#[inline(always)]
pub unsafe fn forward_start(
    nod: *mut Node,
    p1: c_int,
    n2: *mut Node,
    tinf: *const Training,
    flag: c_int,
) {
    let n1: *mut Node = nod.offset(p1 as isize);
    let mut left = (*n1).ndx;
    let right = (*n2).ndx;
    let mut score: f64 = 0.0;

    if (*n1).type_ == STOP {
        left += 2;
        if left >= right {
            return;
        }
        if flag == 1 {
            score = intergenic_mod_same(n1, n2, tinf);
        }
    } else {
        if left >= right {
            return;
        }
        if flag == 1 {
            score = intergenic_mod_diff(tinf);
        }
    }
    commit(n1, n2, p1, score, -1);
}

#[inline(always)]
pub unsafe fn forward_stop(
    nod: *mut Node,
    p1: c_int,
    n2: *mut Node,
    tinf: *const Training,
    flag: c_int,
) {
    let n1: *mut Node = nod.offset(p1 as isize);
    let mut left = (*n1).ndx;
    let mut right = (*n2).ndx;
    let mut score: f64 = 0.0;
    let mut scr_mod: f64 = 0.0;

    if (*n2).stop_val >= (*n1).ndx {
        return;
    }
    if (*n1).type_ != STOP {
        if (*n1).ndx % 3 != (*n2).ndx % 3 {
            return;
        }
        right += 2;
        if flag == 0 {
            scr_mod = frame_bias(n1, tinf);
        } else {
            score = (*n1).cscore + (*n1).sscore;
        }
    } else {
        let star = (*n1).star_ptr[((*n2).ndx % 3) as usize];
        if star == -1 {
            return;
        }
        let n3: *mut Node = nod.offset(star as isize);
        left = (*n3).ndx;
        right += 2;
        if flag == 0 {
            scr_mod = frame_bias(n3, tinf);
        } else {
            score = (*n3).cscore + (*n3).sscore + intergenic_mod_same(n1, n3, tinf);
        }
    }
    if flag == 0 {
        score = ((right - left + 1) as f64) * scr_mod;
    }
    commit(n1, n2, p1, score, -1);
}

#[inline(always)]
pub unsafe fn backward_start(
    nod: *mut Node,
    p1: c_int,
    n2: *mut Node,
    tinf: *const Training,
    flag: c_int,
) {
    let n1: *mut Node = nod.offset(p1 as isize);
    let mut left = (*n1).ndx;
    let right = (*n2).ndx;
    let mut score: f64 = 0.0;
    let mut scr_mod: f64 = 0.0;
    let mut ovlp: c_int = 0;

    if (*n1).strand == -1 {
        if (*n1).stop_val <= (*n2).ndx {
            return;
        }
        if (*n1).ndx % 3 != (*n2).ndx % 3 {
            return;
        }
        left -= 2;
        if flag == 0 {
            scr_mod = frame_bias(n2, tinf);
        } else {
            score = (*n2).cscore + (*n2).sscore;
        }
    } else {
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
        let bnd = if (*n1).traceb == -1 {
            0
        } else {
            (*nod.offset((*n1).traceb as isize)).ndx
        };
        if ((*n1).ndx + 2 - (*n2).stop_val - 2 + 1) >= ((*n2).stop_val - 3 - bnd + 1) {
            return;
        }
        left = (*n2).stop_val - 2;
        if flag == 0 {
            scr_mod = frame_bias(n2, tinf);
        } else {
            score = (*n2).cscore + (*n2).sscore - 0.15 * (*tinf).st_wt;
        }
    }
    if flag == 0 {
        score = ((right - left + 1 - (ovlp * 2)) as f64) * scr_mod;
    }
    commit(n1, n2, p1, score, -1);
}

#[inline(always)]
pub unsafe fn backward_stop(
    nod: *mut Node,
    p1: c_int,
    n2: *mut Node,
    tinf: *const Training,
    flag: c_int,
) {
    let n1: *mut Node = nod.offset(p1 as isize);
    let mut left = (*n1).ndx;
    let mut right = (*n2).ndx;
    let mut score: f64 = 0.0;
    let mut scr_mod: f64 = 0.0;
    let mut ovlp: c_int = 0;
    let mut maxfr: c_int = -1;

    if (*n1).strand == 1 {
        left += 2;
        right -= 2;
        if left >= right {
            return;
        }
        let mut maxval: f64 = 0.0;
        for i in 0..3 {
            if (*n2).star_ptr[i as usize] == -1 {
                continue;
            }
            let n3: *mut Node = nod.offset((*n2).star_ptr[i as usize] as isize);
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
            let better = if flag == 1 {
                (*n3).cscore + (*n3).sscore + intergenic_mod_same(n3, n2, tinf) > maxval
            } else {
                frame_bias(n3, tinf) > maxval
            };
            if better {
                maxfr = i;
                maxval = (*n3).cscore + (*n3).sscore + intergenic_mod_same(n3, n2, tinf);
            }
        }
        if maxfr != -1 {
            let n3: *mut Node = nod.offset((*n2).star_ptr[maxfr as usize] as isize);
            if flag == 0 {
                scr_mod = frame_bias(n3, tinf);
            } else {
                score = (*n3).cscore + (*n3).sscore + intergenic_mod_same(n3, n2, tinf);
            }
        } else if flag == 1 {
            score = intergenic_mod_diff(tinf);
        }
    } else if (*n1).type_ != STOP {
        right -= 2;
        if left >= right {
            return;
        }
        if flag == 1 {
            score = intergenic_mod_same(n1, n2, tinf);
        }
    } else {
        if (*n1).stop_val <= (*n2).ndx {
            return;
        }
        let star = (*n2).star_ptr[((*n1).ndx % 3) as usize];
        if star == -1 {
            return;
        }
        let n3: *mut Node = nod.offset(star as isize);
        left -= 2;
        right = (*n3).ndx;
        if flag == 0 {
            scr_mod = frame_bias(n3, tinf);
        } else {
            score = (*n3).cscore + (*n3).sscore + intergenic_mod_same(n3, n2, tinf);
        }
    }
    if flag == 0 {
        score = ((right - left + 1 - (ovlp * 2)) as f64) * scr_mod;
    }
    commit(n1, n2, p1, score, maxfr);
}

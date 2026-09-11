use prodigal_rs::connection::{backward_start, backward_stop, forward_start, forward_stop};
use prodigal_rs::connection_filter::{allowed_for, class_of};
use prodigal_rs::node::intergenic_mod;
use prodigal_rs::types::{Node, Training, MAX_OPP_OVLP, STOP};
use std::os::raw::c_int;

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 11
    }

    fn below(&mut self, bound: u64) -> u64 {
        self.next() % bound
    }

    fn real(&mut self, span: f64) -> f64 {
        (self.below(20001) as f64 / 10000.0 - 1.0) * span
    }
}

// `record_overlapping_starts` only ever points a stop node at a start on its own strand, and the
// scorers lean on that, so the generator has to honour it too.
fn nodes(seed: u64, count: usize) -> Vec<Node> {
    let mut rng = Lcg(seed);
    let mut built: Vec<Node> = Vec::with_capacity(count);
    for i in 0..count {
        let mut node: Node = unsafe { std::mem::zeroed() };
        node.ndx = rng.below(3000) as c_int;
        node.stop_val = rng.below(3000) as c_int;
        node.strand = if rng.below(2) == 0 { 1 } else { -1 };
        node.type_ = rng.below(4) as c_int;
        node.edge = rng.below(2) as c_int;
        node.traceb = match i == 0 || rng.below(4) == 0 {
            true => -1,
            false => rng.below(i as u64) as c_int,
        };
        node.score = rng.real(50.0);
        node.cscore = rng.real(20.0);
        node.sscore = rng.real(20.0);
        node.rscore = rng.real(5.0);
        node.uscore = rng.real(5.0);
        for slot in 0..3 {
            node.gc_score[slot] = rng.real(2.0);
        }
        node.ov_mark = -1;
        node.star_ptr = [-1; 3];
        built.push(node);
    }
    let starts = |strand: c_int, built: &[Node]| {
        built
            .iter()
            .enumerate()
            .filter(|(_, node)| node.strand == strand && node.type_ != STOP)
            .map(|(index, _)| index as c_int)
            .collect::<Vec<_>>()
    };
    let forward = starts(1, &built);
    let reverse = starts(-1, &built);
    for i in 0..count {
        if built[i].type_ != STOP || built[i].edge == 1 {
            continue;
        }
        let pool = match built[i].strand == 1 {
            true => &forward,
            false => &reverse,
        };
        for slot in 0..3 {
            if pool.is_empty() || rng.below(3) == 0 {
                continue;
            }
            built[i].star_ptr[slot] = pool[rng.below(pool.len() as u64) as usize];
        }
    }
    built
}

fn training(seed: u64) -> Training {
    let mut rng = Lcg(seed);
    let mut tinf: Training = unsafe { std::mem::zeroed() };
    tinf.st_wt = 4.35;
    for slot in 0..3 {
        tinf.bias[slot] = rng.real(1.5);
    }
    tinf
}

fn landed(node: &Node) -> (u64, c_int, c_int) {
    (node.score.to_bits(), node.traceb, node.ov_mark)
}

unsafe fn dispatch(nod: *mut Node, p1: c_int, p2: c_int, tinf: *mut Training, flag: c_int) {
    let n2 = nod.offset(p2 as isize);
    match 2 * u8::from((*n2).strand != 1) + u8::from((*n2).type_ == STOP) {
        0 => forward_start(nod, p1, n2, tinf, flag),
        1 => forward_stop(nod, p1, n2, tinf, flag),
        2 => backward_start(nod, p1, n2, tinf, flag),
        _ => backward_stop(nod, p1, n2, tinf, flag),
    }
}

unsafe fn reference(
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

#[test]
fn the_split_scorers_match_the_original_on_every_allowed_pair() {
    let count = 140;
    let mut checked = [[0usize; 16]; 2];
    let mut rejected = 0usize;
    for flag in [0, 1] {
        for seed in 0..24u64 {
            let base = nodes(seed * 7 + 1, count);
            let mut tinf = training(seed + 3);
            for right in 1..count {
                let table = unsafe { allowed_for(&base[right]) };
                let kind = 2 * usize::from(base[right].strand != 1)
                    + usize::from(base[right].type_ == STOP);
                for left in 0..right {
                    let class = unsafe { class_of(&base[left]) };
                    let mut want = base.clone();
                    unsafe {
                        reference(
                            want.as_mut_ptr(),
                            left as c_int,
                            right as c_int,
                            &mut tinf,
                            flag,
                        )
                    };
                    if table[class as usize] == 0 {
                        assert_eq!(
                            landed(&want[right]),
                            landed(&base[right]),
                            "the filter rejected a pair the original scores"
                        );
                        rejected += 1;
                        continue;
                    }
                    let mut got = base.clone();
                    unsafe {
                        dispatch(
                            got.as_mut_ptr(),
                            left as c_int,
                            right as c_int,
                            &mut tinf,
                            flag,
                        )
                    };
                    assert_eq!(
                        landed(&got[right]),
                        landed(&want[right]),
                        "flag {flag} n1 {:?} n2 {:?}",
                        (base[left].strand, base[left].type_, base[left].ndx, base[left].stop_val, base[left].traceb),
                        (base[right].strand, base[right].type_, base[right].ndx, base[right].stop_val, base[right].star_ptr)
                    );
                    let from = 2 * usize::from(base[left].strand != 1)
                        + usize::from(base[left].type_ == STOP);
                    checked[flag as usize][kind * 4 + from] += 1;
                }
            }
        }
    }
    assert!(rejected > 0);
    let reachable = [1, 2, 4, 5, 9, 11, 13, 14, 15];
    for flag in 0..2 {
        for cell in reachable {
            assert!(
                checked[flag][cell] > 200,
                "flag {flag} pairing {cell} only reached {} times",
                checked[flag][cell]
            );
        }
        let stray = (0..16).filter(|c| !reachable.contains(c)).find(|c| checked[flag][*c] > 0);
        assert_eq!(stray, None, "the filter let an invalid pairing through");
    }
}

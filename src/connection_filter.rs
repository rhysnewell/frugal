/*******************************************************************************
    Connection validity prefilter for the dynamic programming routine.

    `score_connection` opens with six tests that reject a pair outright, and all
    six read only `strand`, `type_` and `traceb` of the left node and `strand`
    and `type_` of the right one. `score_connection` writes only to the right
    node, so during one pass of the inner loop every left node is immutable and
    the whole range can be classified at once.
*******************************************************************************/

use crate::types::{Node, STOP};
use std::os::raw::c_int;

pub const CLASSES: usize = 8;

const FORWARD: u8 = 0b100;
const IS_STOP: u8 = 0b010;
const NO_TRACE: u8 = 0b001;

#[inline]
pub unsafe fn class_of(node: *const Node) -> u8 {
    let mut held = 0u8;
    if (*node).strand == 1 {
        held |= FORWARD;
    }
    if (*node).type_ == STOP {
        held |= IS_STOP;
    }
    if (*node).traceb == -1 {
        held |= NO_TRACE;
    }
    held
}

/// 0xFF where a left node of that class can still reach the scoring body.
pub unsafe fn allowed_for(right: *const Node) -> [u8; 16] {
    let forward = (*right).strand == 1;
    let stop = (*right).type_ == STOP;

    let mut table = [0u8; 16];
    for (class, cell) in table.iter_mut().enumerate().take(CLASSES) {
        let class = class as u8;
        let left_forward = class & FORWARD != 0;
        let left_stop = class & IS_STOP != 0;
        let untraced = class & NO_TRACE != 0;

        let invalid = (!left_stop && !stop && left_forward == forward)
            || (left_forward && !left_stop && !forward)
            || (!left_forward && left_stop && forward)
            || (!left_forward && !left_stop && forward && stop)
            || (untraced && left_forward && left_stop)
            || (untraced && !left_forward && !left_stop);

        *cell = if invalid { 0x00 } else { 0xFF };
    }
    table
}

#[cfg(target_arch = "aarch64")]
#[inline]
fn survivors(classes: &[u8], table: &[u8; 16], out: &mut Vec<c_int>) {
    use std::arch::aarch64::*;
    let mut at = 0;
    unsafe {
        let lookup = vld1q_u8(table.as_ptr());
        while at + 16 <= classes.len() {
            let kept = vqtbl1q_u8(lookup, vld1q_u8(classes.as_ptr().add(at)));
            if vmaxvq_u8(kept) != 0 {
                let mut lanes = [0u8; 16];
                vst1q_u8(lanes.as_mut_ptr(), kept);
                for (lane, held) in lanes.iter().enumerate() {
                    if *held != 0 {
                        out.push((at + lane) as c_int);
                    }
                }
            }
            at += 16;
        }
    }
    tail(classes, table, at, out);
}

#[cfg(not(target_arch = "aarch64"))]
#[inline]
fn survivors(classes: &[u8], table: &[u8; 16], out: &mut Vec<c_int>) {
    tail(classes, table, 0, out);
}

#[inline]
fn tail(classes: &[u8], table: &[u8; 16], from: usize, out: &mut Vec<c_int>) {
    for (at, class) in classes.iter().enumerate().skip(from) {
        if table[*class as usize] != 0 {
            out.push(at as c_int);
        }
    }
}

/// Offsets into `classes` whose node can still reach the scoring body.
pub fn keep(classes: &[u8], table: &[u8; 16], out: &mut Vec<c_int>) {
    out.clear();
    survivors(classes, table, out);
}

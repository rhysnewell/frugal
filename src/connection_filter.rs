// `score_connection` opens with six tests that reject two thirds of the pairs, and all six read
// only three fields of the left node and two of the right. It writes only to the right node, so
// the left ones hold still for a whole pass and their three fields pack into one byte.

use crate::types::{Node, STOP};

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

pub unsafe fn allowed_for(right: *const Node) -> [u8; CLASSES] {
    let forward = (*right).strand == 1;
    let stop = (*right).type_ == STOP;

    let mut table = [0u8; CLASSES];
    for (class, cell) in table.iter_mut().enumerate() {
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

        *cell = u8::from(!invalid);
    }
    table
}

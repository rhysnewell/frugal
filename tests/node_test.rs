use std::mem::MaybeUninit;

use frugal::node::node_room;
use frugal::types::{Node, STT_NOD};

#[test]
fn a_long_contig_gets_room_for_every_node_it_may_hold() {
    let mut nodes: Vec<MaybeUninit<Node>> = Vec::with_capacity(STT_NOD);
    for slen in [1_200_000, 1_700_000, 7_900_000] {
        node_room(&mut nodes, slen);
        assert!(
            nodes.capacity() >= slen as usize / 8,
            "{slen} {}",
            nodes.capacity()
        );
    }
}

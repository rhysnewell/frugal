use std::os::raw::c_int;

use crate::bitmap::set;
use crate::sequence::rcom_seq;
use crate::types::{Gene, Mask, Node, MASK_SIZE, MAX_GENES, MAX_MASKS, MAX_SEQ, STT_NOD};

/// Internal scratch buffers for the prediction pipeline.
pub(crate) struct SequenceBuffer {
    pub seq: Vec<u8>,
    pub rseq: Vec<u8>,
    pub useq: Vec<u8>,
    pub nodes: Vec<Node>,
    pub genes: Vec<Gene>,
    pub masks: Vec<Mask>,
    pub nmask: c_int,
    prev_len: usize,
    prev_nn: usize,
}

impl SequenceBuffer {
    /// Allocate a fresh scratch buffer with capacity for `MAX_SEQ` bases and
    /// the default node/gene/mask counts.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self::with_base_capacity(MAX_SEQ)
    }

    /// Allocate a fresh scratch buffer sized for a known input upper bound.
    pub fn with_base_capacity(max_bases: usize) -> Self {
        let max_bases = max_bases.min(MAX_SEQ);
        let seq_bytes = (max_bases / 4 + 2).max(1);
        let useq_bytes = (max_bases / 8 + 2).max(1);
        SequenceBuffer {
            seq: vec![0u8; seq_bytes],
            rseq: vec![0u8; seq_bytes],
            useq: vec![0u8; useq_bytes],
            nodes: vec![unsafe { std::mem::zeroed() }; STT_NOD],
            genes: vec![unsafe { std::mem::zeroed() }; MAX_GENES],
            masks: vec![unsafe { std::mem::zeroed() }; MAX_MASKS],
            nmask: 0,
            prev_len: 0,
            prev_nn: 0,
        }
    }

    /// Node and gene buffers keep their full size because `add_nodes` and `add_genes` write
    /// through raw pointers with no bound check.
    pub fn reusable() -> Self {
        Self::with_base_capacity(0)
    }

    /// Record how much of the node buffer the caller dirtied, so the next sequence can clear
    /// exactly that prefix.
    pub fn mark_nodes(&mut self, nn: c_int) {
        self.prev_nn = (nn.max(0) as usize).min(self.nodes.len());
    }

    fn prepare(&mut self, dna_len: usize) {
        let max_bases = dna_len.min(MAX_SEQ);
        let seq_bytes = (max_bases / 4 + 2).max(1);
        let useq_bytes = (max_bases / 8 + 2).max(1);
        if self.seq.len() < seq_bytes {
            self.seq.resize(seq_bytes, 0);
            self.rseq.resize(seq_bytes, 0);
        }
        if self.useq.len() < useq_bytes {
            self.useq.resize(useq_bytes, 0);
        }
        self.clear_seq(self.prev_len);
        for i in 0..self.prev_nn {
            self.nodes[i] = unsafe { std::mem::zeroed() };
        }
        self.prev_nn = 0;
    }

    /// Clear sequence buffers for a new sequence.
    fn clear_seq(&mut self, prev_len: usize) {
        let seq_bytes = (prev_len / 4 + 1).min(self.seq.len());
        let useq_bytes = (prev_len / 8 + 1).min(self.useq.len());
        self.seq[..seq_bytes].fill(0);
        self.rseq[..seq_bytes].fill(0);
        self.useq[..useq_bytes].fill(0);
        self.nmask = 0;
    }

    /// Encode a single DNA sequence into the internal 2-bit bitmap format.
    /// Returns (sequence_length, gc_content).
    pub unsafe fn encode(&mut self, dna: &[u8], do_mask: bool) -> (c_int, f64) {
        self.prepare(dna.len());

        let mut bctr: c_int = 0;
        let mut len: c_int = 0;
        let mut gc_count: c_int = 0;
        let mut mask_beg: c_int = -1;

        for &base in dna {
            if base < b'A' || base > b'z' {
                continue;
            }

            // Handle mask tracking
            if do_mask && mask_beg != -1 && base != b'N' && base != b'n' {
                if len - mask_beg >= MASK_SIZE as c_int {
                    if (self.nmask as usize) < MAX_MASKS {
                        self.masks[self.nmask as usize].begin = mask_beg;
                        self.masks[self.nmask as usize].end = len - 1;
                        self.nmask += 1;
                    }
                }
                mask_beg = -1;
            }
            if do_mask && mask_beg == -1 && (base == b'N' || base == b'n') {
                mask_beg = len;
            }

            match base {
                b'g' | b'G' => {
                    set(self.seq.as_mut_ptr(), bctr);
                    gc_count += 1;
                }
                b't' | b'T' => {
                    set(self.seq.as_mut_ptr(), bctr);
                    set(self.seq.as_mut_ptr(), bctr + 1);
                }
                b'c' | b'C' => {
                    set(self.seq.as_mut_ptr(), bctr + 1);
                    gc_count += 1;
                }
                b'a' | b'A' => { /* bits 00, nothing to set */ }
                _ => {
                    // Unknown base (N, etc.)
                    set(self.seq.as_mut_ptr(), bctr + 1);
                    set(self.useq.as_mut_ptr(), len);
                }
            }
            bctr += 2;
            len += 1;
        }

        if len > 0 {
            rcom_seq(
                self.seq.as_mut_ptr(),
                self.rseq.as_mut_ptr(),
                self.useq.as_mut_ptr(),
                len,
            );
        }

        let gc = if len > 0 {
            gc_count as f64 / len as f64
        } else {
            0.0
        };
        self.prev_len = len as usize;
        (len, gc)
    }

    /// Encode multiple sequences concatenated with TTAATTAATTAA spacers.
    /// Used for training (mirrors read_seq_training behavior).
    /// Returns (total_length, gc_content).
    pub unsafe fn encode_training(&mut self, sequences: &[&[u8]], do_mask: bool) -> (c_int, f64) {
        self.clear_seq(self.seq.len() * 4);

        let mut bctr: c_int = 0;
        let mut len: c_int = 0;
        let mut gc_count: c_int = 0;
        let mut mask_beg: c_int = -1;
        let mut seq_count = 0;

        for &seq_data in sequences {
            // Insert TTAATTAATTAA spacer between sequences
            if seq_count > 0 {
                // 12 bases: T T A A T T A A T T A A
                for i in 0..12 {
                    if i % 4 == 0 || i % 4 == 1 {
                        set(self.seq.as_mut_ptr(), bctr);
                        set(self.seq.as_mut_ptr(), bctr + 1);
                    }
                    bctr += 2;
                    len += 1;
                }
            }
            seq_count += 1;

            for &base in seq_data {
                if base < b'A' || base > b'z' {
                    continue;
                }

                if do_mask && mask_beg != -1 && base != b'N' && base != b'n' {
                    if len - mask_beg >= MASK_SIZE as c_int {
                        if (self.nmask as usize) < MAX_MASKS {
                            self.masks[self.nmask as usize].begin = mask_beg;
                            self.masks[self.nmask as usize].end = len - 1;
                            self.nmask += 1;
                        }
                    }
                    mask_beg = -1;
                }
                if do_mask && mask_beg == -1 && (base == b'N' || base == b'n') {
                    mask_beg = len;
                }

                match base {
                    b'g' | b'G' => {
                        set(self.seq.as_mut_ptr(), bctr);
                        gc_count += 1;
                    }
                    b't' | b'T' => {
                        set(self.seq.as_mut_ptr(), bctr);
                        set(self.seq.as_mut_ptr(), bctr + 1);
                    }
                    b'c' | b'C' => {
                        set(self.seq.as_mut_ptr(), bctr + 1);
                        gc_count += 1;
                    }
                    b'a' | b'A' => {}
                    _ => {
                        set(self.seq.as_mut_ptr(), bctr + 1);
                        set(self.useq.as_mut_ptr(), len);
                    }
                }
                bctr += 2;
                len += 1;

                if (len as usize) + 10000 >= MAX_SEQ {
                    break;
                }
            }
        }

        // Add trailing spacer if multiple sequences
        if seq_count > 1 {
            for i in 0..12 {
                if i % 4 == 0 || i % 4 == 1 {
                    set(self.seq.as_mut_ptr(), bctr);
                    set(self.seq.as_mut_ptr(), bctr + 1);
                }
                bctr += 2;
                len += 1;
            }
        }

        if len > 0 {
            rcom_seq(
                self.seq.as_mut_ptr(),
                self.rseq.as_mut_ptr(),
                self.useq.as_mut_ptr(),
                len,
            );
        }

        let gc = if len > 0 {
            gc_count as f64 / len as f64
        } else {
            0.0
        };
        self.prev_len = len as usize;
        (len, gc)
    }

    /// Ensure node buffer is large enough for the given sequence length.
    pub fn ensure_node_capacity(&mut self, slen: c_int) {
        let needed = (slen as usize) / 8;
        if needed > self.nodes.len() {
            self.nodes.resize(needed, unsafe { std::mem::zeroed() });
        }
    }

    /// Clear node buffer.
    pub fn clear_nodes(&mut self, nn: c_int) {
        let dirty = (nn.max(0) as usize).max(self.prev_nn).min(self.nodes.len());
        for i in 0..dirty {
            self.nodes[i] = unsafe { std::mem::zeroed() };
        }
        self.prev_nn = 0;
    }
}

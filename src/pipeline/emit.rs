use std::os::raw::{c_char, c_int};

use super::VERSION_CSTR;
use super::handles::Handles;
use crate::gene::{print_genes, write_nucleotide_seqs, write_translations};
use crate::types::{Gene, Node, Training};

// A null bin description marks a single genome run, as print_genes reads it.
pub(super) unsafe fn write_genes(
    handles: &Handles,
    genes: *mut Gene,
    ng: c_int,
    nodes: *mut Node,
    seq: *mut u8,
    rseq: *mut u8,
    useq: *mut u8,
    slen: c_int,
    format: c_int,
    num_seq: c_int,
    desc: *mut c_char,
    tinf: *const Training,
    cur_header: *mut c_char,
    short_header: *mut c_char,
) {
    print_genes(
        handles.output,
        genes,
        ng,
        nodes,
        slen,
        format,
        num_seq,
        (!desc.is_null()) as c_int,
        desc,
        tinf,
        cur_header,
        short_header,
        VERSION_CSTR.as_ptr() as *mut c_char,
    );
    if handles.trans != handles.stdout {
        write_translations(
            handles.trans,
            genes,
            ng,
            nodes,
            seq,
            rseq,
            useq,
            slen,
            tinf,
            num_seq,
            short_header,
        );
    }
    if handles.nuc != handles.stdout {
        write_nucleotide_seqs(
            handles.nuc,
            genes,
            ng,
            nodes,
            seq,
            rseq,
            useq,
            slen,
            tinf,
            num_seq,
            short_header,
        );
    }
}

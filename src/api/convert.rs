use super::types::{PredictedGene, StartCodon, Strand};
use crate::gene::calculate_confidence;
use crate::types::{Gene, Node, Training};

// SD motif/spacer lookup tables (same as gene.rs record_gene_data)
const SD_STRING: [&str; 28] = [
    "None",
    "GGA/GAG/AGG",
    "3Base/5BMM",
    "4Base/6BMM",
    "AGxAG",
    "AGxAG",
    "GGA/GAG/AGG",
    "GGxGG",
    "GGxGG",
    "AGxAG",
    "AGGAG(G)/GGAGG",
    "AGGA/GGAG/GAGG",
    "AGGA/GGAG/GAGG",
    "GGA/GAG/AGG",
    "GGxGG",
    "AGGA",
    "GGAG/GAGG",
    "AGxAGG/AGGxGG",
    "AGxAGG/AGGxGG",
    "AGxAGG/AGGxGG",
    "AGGAG/GGAGG",
    "AGGAG",
    "AGGAG",
    "GGAGG",
    "GGAGG",
    "AGGAGG",
    "AGGAGG",
    "AGGAGG",
];

const SD_SPACER: [&str; 28] = [
    "None", "3-4bp", "13-15bp", "13-15bp", "11-12bp", "3-4bp", "11-12bp", "11-12bp", "3-4bp",
    "5-10bp", "13-15bp", "3-4bp", "11-12bp", "5-10bp", "5-10bp", "5-10bp", "5-10bp", "11-12bp",
    "3-4bp", "5-10bp", "11-12bp", "3-4bp", "5-10bp", "3-4bp", "5-10bp", "11-12bp", "3-4bp",
    "5-10bp",
];

use crate::sequence::mer_text;

/// Convert internal Gene + Node data to a PredictedGene.
pub(crate) unsafe fn gene_to_predicted(
    gene: &Gene,
    nodes: *const Node,
    tinf: &Training,
    _seq_len: usize,
) -> PredictedGene {
    let n = &*nodes.offset(gene.start_ndx as isize);
    let sn = &*nodes.offset(gene.stop_ndx as isize);

    let strand = if n.strand == 1 {
        Strand::Forward
    } else {
        Strand::Reverse
    };

    let partial_left = (n.edge == 1 && n.strand == 1) || (sn.edge == 1 && n.strand == -1);
    let partial_right = (sn.edge == 1 && n.strand == 1) || (n.edge == 1 && n.strand == -1);
    let begin = gene.begin as usize;
    let end = gene.end as usize;

    let start_codon = if n.edge == 1 {
        StartCodon::Edge
    } else {
        match n.type_ {
            0 => StartCodon::ATG,
            1 => StartCodon::GTG,
            2 => StartCodon::TTG,
            _ => StartCodon::Edge,
        }
    };

    // Determine RBS motif and spacer
    let rbs1_score = tinf.rbs_wt[n.rbs[0] as usize] * tinf.st_wt;
    let rbs2_score = tinf.rbs_wt[n.rbs[1] as usize] * tinf.st_wt;

    let (rbs_motif, rbs_spacer) = if tinf.uses_sd == 1 {
        if rbs1_score > rbs2_score {
            (SD_STRING[n.rbs[0] as usize], SD_SPACER[n.rbs[0] as usize])
        } else {
            (SD_STRING[n.rbs[1] as usize], SD_SPACER[n.rbs[1] as usize])
        }
    } else {
        // Non-SD mode: check if SD motif beats upstream motif
        if tinf.no_mot > -0.5 && rbs1_score > rbs2_score && rbs1_score > n.mot.score * tinf.st_wt {
            (SD_STRING[n.rbs[0] as usize], SD_SPACER[n.rbs[0] as usize])
        } else if tinf.no_mot > -0.5
            && rbs2_score >= rbs1_score
            && rbs2_score > n.mot.score * tinf.st_wt
        {
            (SD_STRING[n.rbs[1] as usize], SD_SPACER[n.rbs[1] as usize])
        } else if n.mot.len == 0 {
            ("None", "None")
        } else {
            // Upstream motif: convert mer index to text
            let mut qt: [std::os::raw::c_char; 10] = [0; 10];
            mer_text(qt.as_mut_ptr(), n.mot.len, n.mot.ndx);
            let motif = std::ffi::CStr::from_ptr(qt.as_ptr())
                .to_str()
                .unwrap_or("None")
                .to_string();
            let spacer = format!("{}bp", n.mot.spacer);
            return PredictedGene {
                begin,
                end,
                strand,
                start_codon,
                translation_table: tinf.trans_table as u8,
                partial: (partial_left, partial_right),
                rbs_motif: motif,
                rbs_spacer: spacer,
                gc_content: n.gc_cont,
                confidence: calculate_confidence(n.cscore + n.sscore, tinf.st_wt),
                score: n.cscore + n.sscore,
                cscore: n.cscore,
                sscore: n.sscore,
                rscore: n.rscore,
                uscore: n.uscore,
                tscore: n.tscore,
            };
        }
    };

    PredictedGene {
        begin,
        end,
        strand,
        start_codon,
        translation_table: tinf.trans_table as u8,
        partial: (partial_left, partial_right),
        rbs_motif: rbs_motif.to_string(),
        rbs_spacer: rbs_spacer.to_string(),
        gc_content: n.gc_cont,
        confidence: calculate_confidence(n.cscore + n.sscore, tinf.st_wt),
        score: n.cscore + n.sscore,
        cscore: n.cscore,
        sscore: n.sscore,
        rscore: n.rscore,
        uscore: n.uscore,
        tscore: n.tscore,
    }
}

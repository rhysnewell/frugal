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

use super::*;

/*******************************************************************************
  Write detailed scoring information about every single possible gene.
*******************************************************************************/

/// Write detailed scoring information about every single possible gene.
/// Only done at the user's request. Emits a tab-separated table per
/// candidate start (one block per stop codon) including total/coding/start
/// scores, RBS motif description, spacer, upstream/type scores and GC
/// content; sorts nodes by stop first, then restores the original ordering
/// before returning.
pub unsafe fn write_start_file(
    fh: c_int,
    nod: *mut Node,
    nn: c_int,
    tinf: *const Training,
    sctr: c_int,
    slen: c_int,
    is_meta: c_int,
    mdesc: *mut c_char,
    version: *mut c_char,
    header: *mut c_char,
) {
    let mut prev_stop: c_int = -1;
    let mut prev_strand: c_int = 0;
    let mut st_type: c_int;
    let mut rbs1: f64;
    let mut rbs2: f64;

    static SD_STRING: [&str; 28] = [
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
    static SD_SPACER: [&str; 28] = [
        "None", "3-4bp", "13-15bp", "13-15bp", "11-12bp", "3-4bp", "11-12bp", "11-12bp", "3-4bp",
        "5-10bp", "13-15bp", "3-4bp", "11-12bp", "5-10bp", "5-10bp", "5-10bp", "5-10bp", "11-12bp",
        "3-4bp", "5-10bp", "11-12bp", "3-4bp", "5-10bp", "3-4bp", "5-10bp", "11-12bp", "3-4bp",
        "5-10bp",
    ];
    static TYPE_STRING: [&str; 4] = ["ATG", "GTG", "TTG", "Edge"];

    let mut qt: [c_char; 10] = [0; 10];

    let header_str = cstr(header);
    let version_str = cstr(version);

    /* Initialize sequence data */
    let seq_data = format!("seqnum={};seqlen={};seqhdr=\"{}\"", sctr, slen, header_str);

    /* Initialize run data string */
    let mut run_data;
    if is_meta == 0 {
        run_data = format!("version=Prodigal.v{};run_type=Single;", version_str);
        run_data.push_str("model=\"Ab initio\";");
    } else {
        let mdesc_str = cstr(mdesc);
        run_data = format!("version=Prodigal.v{};run_type=Metagenomic;", version_str);
        run_data.push_str(&format!("model=\"{}\";", mdesc_str));
    }
    run_data.push_str(&format!(
        "gc_cont={:.2};transl_table={};uses_sd={}",
        (*tinf).gc * 100.0,
        (*tinf).trans_table,
        (*tinf).uses_sd
    ));

    {
        let nodes_slice = std::slice::from_raw_parts_mut(nod, nn as usize);
        nodes_slice.sort_unstable_by(|a, b| {
            a.stop_val
                .cmp(&b.stop_val)
                .then(b.strand.cmp(&a.strand))
                .then(a.ndx.cmp(&b.ndx))
        });
    }

    fprint(fh, &format!("# Sequence Data: {}\n", seq_data));
    fprint(fh, &format!("# Run Data: {}\n\n", run_data));

    fprint(fh, "Beg\tEnd\tStd\tTotal\tCodPot\tStrtSc\tCodon\tRBSMot\t");
    fprint(fh, "Spacer\tRBSScr\tUpsScr\tTypeScr\tGCCont\n");

    for i in 0..nn {
        let ni = &*nod.offset(i as isize);
        if ni.type_ == STOP {
            continue;
        }
        if ni.edge == 1 {
            st_type = 3;
        } else {
            st_type = ni.type_;
        }
        if ni.stop_val != prev_stop || ni.strand != prev_strand {
            prev_stop = ni.stop_val;
            prev_strand = ni.strand;
            fprint(fh, "\n");
        }
        if ni.strand == 1 {
            fprint(
                fh,
                &format!(
                    "{}\t{}\t+\t{:.2}\t{:.2}\t{:.2}\t{}\t",
                    ni.ndx + 1,
                    ni.stop_val + 3,
                    ni.cscore + ni.sscore,
                    ni.cscore,
                    ni.sscore,
                    TYPE_STRING[st_type as usize],
                ),
            );
        }
        if ni.strand == -1 {
            fprint(
                fh,
                &format!(
                    "{}\t{}\t-\t{:.2}\t{:.2}\t{:.2}\t{}\t",
                    ni.stop_val - 1,
                    ni.ndx + 1,
                    ni.cscore + ni.sscore,
                    ni.cscore,
                    ni.sscore,
                    TYPE_STRING[st_type as usize],
                ),
            );
        }
        rbs1 = (*tinf).rbs_wt[ni.rbs[0] as usize] * (*tinf).st_wt;
        rbs2 = (*tinf).rbs_wt[ni.rbs[1] as usize] * (*tinf).st_wt;
        if (*tinf).uses_sd == 1 {
            if rbs1 > rbs2 {
                fprint(
                    fh,
                    &format!(
                        "{}\t{}\t{:.2}\t",
                        SD_STRING[ni.rbs[0] as usize], SD_SPACER[ni.rbs[0] as usize], ni.rscore,
                    ),
                );
            } else {
                fprint(
                    fh,
                    &format!(
                        "{}\t{}\t{:.2}\t",
                        SD_STRING[ni.rbs[1] as usize], SD_SPACER[ni.rbs[1] as usize], ni.rscore,
                    ),
                );
            }
        } else {
            mer_text(qt.as_mut_ptr(), ni.mot.len, ni.mot.ndx);
            let qt_str = cstr(qt.as_ptr());
            if (*tinf).no_mot > -0.5 && rbs1 > rbs2 && rbs1 > ni.mot.score * (*tinf).st_wt {
                fprint(
                    fh,
                    &format!(
                        "{}\t{}\t{:.2}\t",
                        SD_STRING[ni.rbs[0] as usize], SD_SPACER[ni.rbs[0] as usize], ni.rscore,
                    ),
                );
            } else if (*tinf).no_mot > -0.5 && rbs2 >= rbs1 && rbs2 > ni.mot.score * (*tinf).st_wt {
                fprint(
                    fh,
                    &format!(
                        "{}\t{}\t{:.2}\t",
                        SD_STRING[ni.rbs[1] as usize], SD_SPACER[ni.rbs[1] as usize], ni.rscore,
                    ),
                );
            } else {
                if ni.mot.len == 0 {
                    fprint(fh, &format!("None\tNone\t{:.2}\t", ni.rscore));
                } else {
                    fprint(
                        fh,
                        &format!("{}\t{}bp\t{:.2}\t", qt_str, ni.mot.spacer, ni.rscore,),
                    );
                }
            }
        }
        fprint(
            fh,
            &format!("{:.2}\t{:.2}\t{:.3}\n", ni.uscore, ni.tscore, ni.gc_cont,),
        );
    }
    fprint(fh, "\n");

    {
        let nodes_slice = std::slice::from_raw_parts_mut(nod, nn as usize);
        nodes_slice.sort_unstable_by(|a, b| a.ndx.cmp(&b.ndx).then(b.strand.cmp(&a.strand)));
    }
}

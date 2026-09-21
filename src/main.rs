#![allow(unsafe_op_in_unsafe_fn)]
#![allow(clippy::manual_range_contains, clippy::unnecessary_map_or)]

use clap::Parser;
use std::process;

/// Prokaryotic gene prediction using dynamic programming.
///
/// Prodigal v2.6.3 — fully rewritten in Rust.
#[derive(Parser, Debug)]
#[command(name = "frugal", version = "2.6.3")]
struct Cli {
    /// Write protein translations to the selected file
    #[arg(short = 'a')]
    trans_file: Option<String>,

    /// Closed ends — do not allow genes to run off edges
    #[arg(short = 'c')]
    closed: bool,

    /// Write nucleotide sequences of genes to the selected file
    #[arg(short = 'd')]
    nuc_file: Option<String>,

    /// Output format: gbk, gff, sco, or gca (default gbk)
    #[arg(short = 'f')]
    output_format: Option<String>,

    /// Translation table (default 11)
    #[arg(short = 'g')]
    trans_table: Option<i32>,

    /// Input FASTA/Genbank file (default stdin)
    #[arg(short = 'i')]
    input_file: Option<String>,

    /// Treat runs of N as masked sequence
    #[arg(short = 'm')]
    mask: bool,

    /// Bypass Shine-Dalgarno trainer and force full motif scan
    #[arg(short = 'n')]
    force_nonsd: bool,

    /// Output file (default stdout)
    #[arg(short = 'o')]
    output_file: Option<String>,

    /// Procedure: single or meta (default single)
    #[arg(short = 'p')]
    mode: Option<String>,

    /// Run quietly (suppress normal stderr output)
    #[arg(short = 'q')]
    quiet: bool,

    /// Write all potential genes (with scores) to the selected file
    #[arg(short = 's')]
    start_file: Option<String>,

    /// Training file (read if exists, write if not)
    #[arg(short = 't')]
    train_file: Option<String>,
}

fn normalize_c_short_options() -> Vec<String> {
    std::env::args()
        .map(|arg| match arg.as_str() {
            "-A" => "-a".to_string(),
            "-C" => "-c".to_string(),
            "-D" => "-d".to_string(),
            "-F" => "-f".to_string(),
            "-G" => "-g".to_string(),
            "-H" => "-h".to_string(),
            "-I" => "-i".to_string(),
            "-M" => "-m".to_string(),
            "-N" => "-n".to_string(),
            "-O" => "-o".to_string(),
            "-P" => "-p".to_string(),
            "-Q" => "-q".to_string(),
            "-S" => "-s".to_string(),
            "-T" => "-t".to_string(),
            "-v" => "--version".to_string(),
            "-V" => "--version".to_string(),
            _ => arg,
        })
        .collect()
}

/// CLI entry point: parses command-line arguments, validates the output format,
/// translation table, and meta/single mode, then dispatches to `run_pipeline`.
fn main() {
    let cli = Cli::parse_from(normalize_c_short_options());

    let output_format = match cli.output_format.as_deref() {
        None => 0,
        Some(s) => match s.to_lowercase().as_str() {
            "0" | "gbk" => 0,
            "1" | "gca" => 1,
            "2" | "sco" => 2,
            "3" | "gff" => 3,
            _ => {
                eprintln!("\nInvalid output format specified.");
                process::exit(15);
            }
        },
    };

    let trans_table = match cli.trans_table {
        None => 0,
        Some(tt) => {
            if tt < 1 || tt > 25 || tt == 7 || tt == 8 || (tt >= 17 && tt <= 20) {
                eprintln!("\nInvalid translation table specified.");
                process::exit(15);
            }
            tt
        }
    };

    let is_meta = match cli.mode.as_deref() {
        None => false,
        Some(s) => match s.as_bytes().first() {
            Some(b'0') | Some(b's') | Some(b'S') => false,
            Some(b'1') | Some(b'm') | Some(b'M') => true,
            _ => {
                eprintln!("\nInvalid meta/single genome type specified.");
                process::exit(15);
            }
        },
    };

    let config = frugal::pipeline::PipelineConfig {
        input_file: cli.input_file,
        output_file: cli.output_file,
        trans_file: cli.trans_file,
        nuc_file: cli.nuc_file,
        start_file: cli.start_file,
        train_file: cli.train_file,
        output_format,
        trans_table,
        closed: cli.closed,
        do_mask: cli.mask,
        force_nonsd: cli.force_nonsd,
        is_meta,
        quiet: cli.quiet,
    };

    let rc = unsafe { frugal::pipeline::run_pipeline(&config) };

    process::exit(rc);
}

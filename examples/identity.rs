use std::io::{BufWriter, Write};
use std::sync::Arc;

use frugal::api::{META_PREDICTOR_STACK_SIZE, MetaPredictor, ProdigalConfig, Strand};

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("fasta path");
    let min_length: usize = args.next().map(|v| v.parse().unwrap()).unwrap_or(1500);
    let threads: usize = args.next().map(|v| v.parse().unwrap()).unwrap_or(4);
    let depth: usize = args.next().map(|v| v.parse().unwrap()).unwrap_or(0);

    let text = std::fs::read_to_string(&path).expect("read fasta");
    let mut contigs: Vec<Vec<u8>> = Vec::new();
    let mut current: Vec<u8> = Vec::new();
    for line in text.lines() {
        if line.starts_with('>') {
            if !current.is_empty() {
                contigs.push(std::mem::take(&mut current));
            }
        } else {
            current.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !current.is_empty() {
        contigs.push(current);
    }
    contigs.retain(|c| c.len() >= min_length);

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(META_PREDICTOR_STACK_SIZE)
        .build()
        .unwrap();
    let config = ProdigalConfig {
        model_depth: depth,
        ..ProdigalConfig::default()
    };
    let predictor = MetaPredictor::with_config_and_thread_pool(config, Arc::new(pool)).unwrap();

    let slices = contigs.iter().map(|c| c.as_slice()).collect::<Vec<_>>();
    let started = std::time::Instant::now();
    let batches = predictor.predict_batch(&slices).unwrap();
    let elapsed = started.elapsed();

    let mut total = 0usize;
    let mut sink = BufWriter::new(std::io::stdout().lock());
    for (index, genes) in batches.iter().enumerate() {
        for gene in genes {
            total += 1;
            let strand = match gene.strand {
                Strand::Forward => '+',
                Strand::Reverse => '-',
            };
            writeln!(
                sink,
                "{index}\t{}\t{}\t{strand}\t{}\t{}\t{}\t{:.6}\t{:.6}",
                gene.begin,
                gene.end,
                gene.partial.0,
                gene.partial.1,
                gene.translation_table,
                gene.score,
                gene.confidence
            )
            .unwrap();
        }
    }
    sink.flush().unwrap();
    eprintln!(
        "contigs {} genes {} wall {:.2}s",
        contigs.len(),
        total,
        elapsed.as_secs_f64()
    );
}

#![allow(clippy::collapsible_else_if, clippy::manual_strip)]

use std::sync::Arc;

use frugal::{
    predict, predict_meta, train, MetaPredictor, ProdigalConfig, StartCodon, Strand,
    META_PREDICTOR_STACK_SIZE,
};

/// Read the sample FASTA and extract raw sequence bytes (stripping headers/newlines).
fn load_sample_sequences() -> Vec<(String, Vec<u8>)> {
    let fasta = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("Prodigal/anthus_aco.fas"),
    )
    .unwrap();

    let mut seqs = Vec::new();
    let mut name = String::new();
    let mut seq = Vec::new();

    for line in fasta.lines() {
        if line.starts_with('>') {
            if !seq.is_empty() {
                seqs.push((name.clone(), seq.clone()));
                seq.clear();
            }
            name = line[1..].to_string();
        } else {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !seq.is_empty() {
        seqs.push((name, seq));
    }
    seqs
}

/// Concatenate all sequences into one (for training).
fn concat_sequences(seqs: &[(String, Vec<u8>)]) -> Vec<u8> {
    let mut all = Vec::new();
    for (_, s) in seqs {
        all.extend_from_slice(s);
    }
    all
}

#[test]
fn test_predict_meta_returns_genes() {
    let seqs = load_sample_sequences();
    // Pick a sequence that has genes
    for (name, seq) in &seqs {
        let genes = predict_meta(seq).unwrap();
        if genes.is_empty() {
            continue;
        }
        eprintln!(
            "[meta] {} ({} bp): {} genes predicted",
            name,
            seq.len(),
            genes.len()
        );

        // Verify basic properties
        for g in &genes {
            assert!(g.begin >= 1);
            assert!(g.end >= 1);
            assert!(g.begin != g.end);
            assert!(g.confidence >= 50.0);
            assert!(g.confidence <= 100.0);
            assert!(g.gc_content >= 0.0 && g.gc_content <= 1.0);
        }
        return;
    }
    panic!("No sequence produced genes");
}

#[test]
fn test_predict_meta_all_sequences() {
    let seqs = load_sample_sequences();
    let mut total_genes = 0;

    for (name, seq) in &seqs {
        let genes = predict_meta(seq).unwrap();
        total_genes += genes.len();
        eprintln!("  {} ({} bp): {} genes", name, seq.len(), genes.len());
    }
    eprintln!(
        "[meta all] total: {} genes across {} sequences",
        total_genes,
        seqs.len()
    );
    assert!(total_genes > 0);
}

#[test]
fn test_meta_predictor_accepts_shared_thread_pool() {
    let pool = Arc::new(
        rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .stack_size(META_PREDICTOR_STACK_SIZE)
            .build()
            .unwrap(),
    );
    let predictor_a = MetaPredictor::with_thread_pool(pool.clone()).unwrap();
    let predictor_b =
        MetaPredictor::with_config_and_thread_pool(ProdigalConfig::default(), pool).unwrap();

    let seqs = load_sample_sequences();
    let seq = &seqs[0].1;
    let genes_a = predictor_a.predict(seq).unwrap();
    let genes_b = predictor_b.predict(seq).unwrap();
    let direct = predict_meta(seq).unwrap();

    assert_eq!(genes_a.len(), genes_b.len());
    assert_eq!(genes_a.len(), direct.len());
    for (a, b) in genes_a.iter().zip(&genes_b) {
        assert_eq!((a.begin, a.end, a.strand), (b.begin, b.end, b.strand));
    }
    for (a, b) in genes_a.iter().zip(&direct) {
        assert_eq!((a.begin, a.end, a.strand), (b.begin, b.end, b.strand));
    }
}

#[test]
fn test_meta_predictor_preserves_reverse_edge_starts() {
    let cases = [
        (
            b"ATCGACTGATGTTCATATTCGCTGCCGATATAAACCTGCGCCGCCACCGCGCAACTGTTCAGGCGCACGGCGTCATCCATCGATAACGCCACGG".as_slice(),
            3,
            92,
        ),
        (
            b"CCAGAGAAGTCAGTGGGTTACGCAGAGAATCGTAACGGGATACCAGATCAAGAGCCTGTCCAATGTGCATAAAAAAATCCGGAAACAA".as_slice(),
            2,
            88,
        ),
        (
            b"ACCCGGCAACGGCAAAGGCCGCGAGAATACCGACGACTGTTGCGATAAGCAGACGGTGGAACATAGTGCGCAGATCCGGATAGATGTGTAGATGATGCATG".as_slice(),
            2,
            100,
        ),
    ];

    let predictor = MetaPredictor::new().unwrap();
    for (seq, begin, end) in cases {
        for genes in [predict_meta(seq).unwrap(), predictor.predict(seq).unwrap()] {
            assert_eq!(genes.len(), 1);
            assert_eq!((genes[0].begin, genes[0].end), (begin, end));
            assert_eq!(genes[0].strand, Strand::Reverse);
            assert_eq!(genes[0].start_codon, StartCodon::Edge);
            assert_eq!(genes[0].partial, (true, true));
        }
    }
}

#[test]
fn test_meta_predictor_preserves_reverse_internal_starts() {
    let seq = b"GGCTTTCGGAAAGCATTCCATGAGCGCGGGGTAGCCACCCGTATCCAGCGCTGAGTTAGCGATACCTACGCACACTGCCAGACCGTAGGCGAGAGTTAAATTCGGGCAAGCGGGAATACCAAAGAAGAATAGCAGATACATTATTACTGCCATTAATATCACCGCCCGACGACCAAACTTATCGGAGATCACACCGAAGAATAAAATACTGATCAATCGCCCCAAACCGATACCGGAAATTAAGTAGGCAATGCCCGCGTTGTCAGTGGAAAACTTTTCCGCCAGAGATGACATATTTTGGGCAAGCGTAATAACACTAATGCCGTGCAGGAAGTAGCTGAAGTAAATACAAAGAACAGCCAGGATAAATGGCGTGCTGAAAGCCTTATTTTGAGACATTTTCACAGCTCCTTGCTGAAAGATAATTTCGTTAC";
    let predictor = MetaPredictor::new().unwrap();
    for genes in [predict_meta(seq).unwrap(), predictor.predict(seq).unwrap()] {
        assert_eq!(genes.len(), 1);
        assert_eq!((genes[0].begin, genes[0].end), (1, 399));
        assert_eq!(genes[0].strand, Strand::Reverse);
        assert_eq!(genes[0].start_codon, StartCodon::ATG);
        assert_eq!(genes[0].partial, (true, false));
    }
}

#[test]
fn test_meta_predictor_reports_selected_translation_table() {
    let seq = b"ATTATTCAGCCTGGTGATAATGTGACGGCAGAGCAAATGATAGCTTCATGAAAAGCTGACTTCGCAATAAGCGCGCAAGAAAATGTAACGTTGGCTCGTGTGGA";
    let genes = MetaPredictor::new().unwrap().predict(seq).unwrap();

    assert_eq!(genes.len(), 1);
    assert_eq!((genes[0].begin, genes[0].end), (1, 102));
    assert_eq!(genes[0].strand, Strand::Forward);
    assert_eq!(genes[0].start_codon, StartCodon::Edge);
    assert_eq!(genes[0].partial, (true, true));
    assert_eq!(genes[0].translation_table, 4);
}

#[test]
fn test_train_and_predict() {
    let seqs = load_sample_sequences();
    let all = concat_sequences(&seqs);

    // Train on the full dataset
    let training = train(&all).unwrap();
    eprintln!(
        "[train] GC={:.2}%, trans_table={}, uses_sd={}",
        training.gc() * 100.0,
        training.translation_table(),
        training.uses_sd()
    );

    // Predict on each sequence
    let mut total_genes = 0;
    for (name, seq) in &seqs {
        let genes = predict(seq, &training).unwrap();
        total_genes += genes.len();
        eprintln!("  {} ({} bp): {} genes", name, seq.len(), genes.len());
    }
    eprintln!("[single] total: {} genes", total_genes);
    assert!(total_genes > 0);
}

#[test]
fn test_training_save_load_roundtrip() {
    let seqs = load_sample_sequences();
    let all = concat_sequences(&seqs);

    let training = train(&all).unwrap();

    let dir = tempfile::TempDir::new().unwrap();
    let path = dir.path().join("test.trn");

    training.save(&path).unwrap();
    let loaded = frugal::TrainingData::load(&path).unwrap();

    assert_eq!(training.gc(), loaded.gc());
    assert_eq!(training.translation_table(), loaded.translation_table());
    assert_eq!(training.uses_sd(), loaded.uses_sd());

    // Predictions should be identical
    let seq = &seqs[0].1;
    let genes1 = predict(seq, &training).unwrap();
    let genes2 = predict(seq, &loaded).unwrap();

    assert_eq!(genes1.len(), genes2.len());
    for (g1, g2) in genes1.iter().zip(genes2.iter()) {
        assert_eq!(g1.begin, g2.begin);
        assert_eq!(g1.end, g2.end);
    }
}

#[test]
fn test_empty_sequence_error() {
    let result = predict_meta(b"");
    assert!(result.is_err());
}

#[test]
fn test_short_sequence_training_error() {
    let result = train(b"ATGCGATCG");
    assert!(result.is_err());
}

#[test]
fn test_gene_properties() {
    let seqs = load_sample_sequences();
    for (_, seq) in &seqs {
        let genes = predict_meta(seq).unwrap();
        if genes.is_empty() {
            continue;
        }
        let g = &genes[0];

        // Check strand
        assert!(g.strand == Strand::Forward || g.strand == Strand::Reverse);

        // Check RBS motif is non-empty
        assert!(!g.rbs_motif.is_empty());

        // Check scores are finite
        assert!(g.score.is_finite());
        assert!(g.cscore.is_finite());
        assert!(g.sscore.is_finite());
        return;
    }
}

#[test]
fn test_custom_config() {
    let seqs = load_sample_sequences();
    let seq = &seqs[0].1; // use first sequence

    let config = ProdigalConfig {
        closed_ends: true,
        ..Default::default()
    };
    let genes = frugal::predict_meta_with(seq, &config).unwrap();
    // With closed ends, edge genes are suppressed
    for g in &genes {
        assert!(
            !g.partial.0 && !g.partial.1,
            "closed_ends should prevent partial genes"
        );
    }
}

// ── Real genome tests ────────────────────────────────────────────────────────

/// Helper: resolve a path relative to CARGO_MANIFEST_DIR's parent.
fn sibling_repo_path(relative: &str) -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join(relative)
}

/// Load a FASTA file, truncating each contig to at most `max_bp` bases.
/// Use `0` for no limit.
fn load_fasta_contigs(path: &std::path::Path, max_bp: usize) -> Vec<(String, Vec<u8>)> {
    let content = std::fs::read_to_string(path).unwrap();
    let mut contigs = Vec::new();
    let mut name = String::new();
    let mut seq = Vec::new();

    for line in content.lines() {
        if line.starts_with('>') {
            if !seq.is_empty() {
                contigs.push((name.clone(), seq.clone()));
                seq.clear();
            }
            name = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
        } else {
            if max_bp == 0 || seq.len() < max_bp {
                let remaining = if max_bp == 0 {
                    usize::MAX
                } else {
                    max_bp - seq.len()
                };
                let bytes = line.trim().as_bytes();
                seq.extend_from_slice(&bytes[..bytes.len().min(remaining)]);
            }
        }
    }
    if !seq.is_empty() {
        contigs.push((name, seq));
    }
    contigs
}

/// Size limit for real genome tests: 100 KB of sequence, enough for ~100 genes
/// but runs in seconds instead of minutes.
const TEST_SEQ_LIMIT: usize = 100_000;

/// Load a FASTA file into a single sequence (concatenating all contigs).
fn load_fasta_sequence(path: &std::path::Path) -> Vec<u8> {
    let content = std::fs::read_to_string(path).unwrap();
    content
        .lines()
        .filter(|l| !l.starts_with('>'))
        .flat_map(|l| l.trim().bytes())
        .collect()
}

/// Test that meta mode produces identical gene coordinates to
/// the C Prodigal binary on a real bacterial plasmid (E. faecium AUS0004_p1).
///
/// Reference coordinates from: `prodigal -c -m -g 11 -p meta -f sco`
/// on prokka/test/plasmid.fna (56520 bp, 63 genes).
#[test]
fn test_prokka_plasmid_meta_coordinates() {
    let plasmid_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../prokka-rs/prokka/test/plasmid.fna");
    if !plasmid_path.exists() {
        eprintln!(
            "Skipping: prokka plasmid test data not found at {}",
            plasmid_path.display()
        );
        return;
    }

    let seq = load_fasta_sequence(&plasmid_path);
    assert_eq!(seq.len(), 56520, "Expected plasmid length 56520 bp");

    let config = ProdigalConfig {
        translation_table: 11,
        closed_ends: true,
        mask_n_runs: true,
        ..Default::default()
    };
    let genes = frugal::predict_meta_with(&seq, &config).unwrap();

    // Reference: 63 genes from C Prodigal v2.6.3 meta mode
    // Format: (begin, end, strand)
    let reference: Vec<(usize, usize, Strand)> = vec![
        (1, 1041, Strand::Forward),
        (1878, 2210, Strand::Forward),
        (2801, 3130, Strand::Forward),
        (3099, 3458, Strand::Forward),
        (3572, 3811, Strand::Forward),
        (3878, 4117, Strand::Forward),
        (4244, 4933, Strand::Forward),
        (4939, 6105, Strand::Forward),
        (6145, 6858, Strand::Forward),
        (6911, 8887, Strand::Forward),
        (8933, 9652, Strand::Forward),
        (9668, 10429, Strand::Forward),
        (10442, 10702, Strand::Forward),
        (10699, 12789, Strand::Forward),
        (12825, 14927, Strand::Reverse),
        (15026, 15217, Strand::Reverse),
        (15420, 15860, Strand::Forward),
        (15857, 16099, Strand::Forward),
        (16117, 16419, Strand::Forward),
        (16488, 18542, Strand::Forward),
        (18542, 21154, Strand::Forward),
        (21151, 23529, Strand::Forward),
        (23590, 24132, Strand::Forward),
        (24210, 24542, Strand::Forward),
        (24543, 25136, Strand::Forward),
        (25158, 27182, Strand::Forward),
        (27182, 28381, Strand::Forward),
        (28397, 29041, Strand::Forward),
        (28998, 29858, Strand::Forward),
        (29880, 30110, Strand::Forward),
        (30127, 30216, Strand::Forward),
        (30677, 30955, Strand::Forward),
        (30973, 31395, Strand::Forward),
        (31414, 33006, Strand::Forward),
        (33183, 33377, Strand::Forward),
        (33390, 33857, Strand::Forward),
        (33796, 35763, Strand::Forward),
        (35839, 37284, Strand::Forward),
        (37292, 39499, Strand::Forward),
        (39515, 39664, Strand::Forward),
        (39697, 40341, Strand::Forward),
        (40354, 41253, Strand::Forward),
        (41256, 41741, Strand::Forward),
        (42107, 42367, Strand::Forward),
        (43019, 43456, Strand::Forward),
        (43449, 44156, Strand::Forward),
        (44537, 44716, Strand::Forward),
        (44891, 45295, Strand::Forward),
        (45303, 46460, Strand::Forward),
        (46707, 46976, Strand::Forward),
        (46969, 47226, Strand::Forward),
        (47789, 49192, Strand::Forward),
        (49189, 49719, Strand::Forward),
        (49714, 50319, Strand::Reverse),
        (50633, 51055, Strand::Forward),
        (51065, 51406, Strand::Forward),
        (51469, 52062, Strand::Forward),
        (52361, 53758, Strand::Forward),
        (53724, 54101, Strand::Forward),
        (54098, 54277, Strand::Forward),
        (54413, 54703, Strand::Forward),
        (55056, 55859, Strand::Forward),
        (55846, 56172, Strand::Forward),
    ];

    assert_eq!(
        genes.len(),
        reference.len(),
        "Gene count mismatch: got {} expected {} (C Prodigal reference)",
        genes.len(),
        reference.len()
    );

    let mut mismatches = 0;
    for (i, (gene, (ref_begin, ref_end, ref_strand))) in
        genes.iter().zip(reference.iter()).enumerate()
    {
        if gene.begin != *ref_begin || gene.end != *ref_end || gene.strand != *ref_strand {
            eprintln!(
                "Gene {}: got {}..{} {} expected {}..{} {}",
                i + 1,
                gene.begin,
                gene.end,
                gene.strand,
                ref_begin,
                ref_end,
                ref_strand
            );
            mismatches += 1;
        }
    }

    assert_eq!(
        mismatches, 0,
        "{} genes have different coordinates vs C Prodigal reference",
        mismatches
    );
}

// ── Real bacterial genome: prokka test genome (NC_011770, ~6.5 MB) ──────────

#[test]
fn test_prokka_genome_meta_mode() {
    let path = sibling_repo_path("prokka-rs/prokka/test/genome.fna");
    if !path.exists() {
        eprintln!("Skipping: prokka genome not found at {}", path.display());
        return;
    }

    let contigs = load_fasta_contigs(&path, TEST_SEQ_LIMIT);
    assert!(!contigs.is_empty(), "No contigs loaded from genome.fna");

    let mut total_genes = 0;
    let mut total_bp = 0;
    for (name, seq) in &contigs {
        let genes = predict_meta(seq).unwrap();
        total_genes += genes.len();
        total_bp += seq.len();
        eprintln!(
            "  [meta] {} ({} bp): {} genes",
            name,
            seq.len(),
            genes.len()
        );

        // Basic sanity on every gene
        for g in &genes {
            assert!(g.begin >= 1 && g.end >= 1);
            assert!(g.begin != g.end);
            assert!(g.confidence >= 0.0 && g.confidence <= 100.0);
            assert!(g.gc_content >= 0.0 && g.gc_content <= 1.0);
            assert!(g.score.is_finite());
        }
    }
    eprintln!(
        "[prokka genome meta] {} genes across {} contigs ({} bp total)",
        total_genes,
        contigs.len(),
        total_bp
    );
    assert!(
        total_genes > 20,
        "Expected >20 genes from prokka genome subset, got {}",
        total_genes
    );
}

#[test]
fn test_prokka_genome_single_mode() {
    let path = sibling_repo_path("prokka-rs/prokka/test/genome.fna");
    if !path.exists() {
        eprintln!("Skipping: prokka genome not found at {}", path.display());
        return;
    }

    let contigs = load_fasta_contigs(&path, TEST_SEQ_LIMIT);
    let all: Vec<u8> = contigs
        .iter()
        .flat_map(|(_, s)| s.iter().copied())
        .collect();

    let training = train(&all).unwrap();
    eprintln!(
        "[prokka genome train] GC={:.2}%, trans_table={}, uses_sd={}",
        training.gc() * 100.0,
        training.translation_table(),
        training.uses_sd()
    );

    // GC content sanity: should be between 20% and 80% for any real genome
    assert!(training.gc() > 0.20 && training.gc() < 0.80);

    let mut total_genes = 0;
    for (name, seq) in &contigs {
        let genes = predict(seq, &training).unwrap();
        total_genes += genes.len();
        eprintln!(
            "  [single] {} ({} bp): {} genes",
            name,
            seq.len(),
            genes.len()
        );
    }
    eprintln!("[prokka genome single] {} total genes", total_genes);
    assert!(
        total_genes > 20,
        "Expected >20 genes from prokka genome subset in single mode, got {}",
        total_genes
    );
}

// ── Real bacterial genome: Priestia megaterium CP157504.1 (~5.2 MB) ─────────

#[test]
fn test_priestia_megaterium_meta_mode() {
    let path = sibling_repo_path("gecco-rs/data/CP157504.1.fna");
    if !path.exists() {
        eprintln!(
            "Skipping: Priestia megaterium genome not found at {}",
            path.display()
        );
        return;
    }

    let contigs = load_fasta_contigs(&path, TEST_SEQ_LIMIT);
    assert!(!contigs.is_empty());

    let mut total_genes = 0;
    for (name, seq) in &contigs {
        let genes = predict_meta(seq).unwrap();
        total_genes += genes.len();
        eprintln!(
            "  [meta] {} ({} bp): {} genes",
            name,
            seq.len(),
            genes.len()
        );

        for g in &genes {
            assert!(g.begin >= 1 && g.end >= 1);
            assert!(g.begin != g.end);
            assert!(g.score.is_finite());
        }
    }
    eprintln!(
        "[priestia meta] {} genes across {} contigs",
        total_genes,
        contigs.len()
    );
    assert!(
        total_genes > 20,
        "Expected >20 genes from Priestia megaterium subset, got {}",
        total_genes
    );
}

#[test]
fn test_priestia_megaterium_single_mode() {
    let path = sibling_repo_path("gecco-rs/data/CP157504.1.fna");
    if !path.exists() {
        eprintln!(
            "Skipping: Priestia megaterium genome not found at {}",
            path.display()
        );
        return;
    }

    let contigs = load_fasta_contigs(&path, TEST_SEQ_LIMIT);
    let all: Vec<u8> = contigs
        .iter()
        .flat_map(|(_, s)| s.iter().copied())
        .collect();

    let training = train(&all).unwrap();
    eprintln!(
        "[priestia train] GC={:.2}%, trans_table={}, uses_sd={}",
        training.gc() * 100.0,
        training.translation_table(),
        training.uses_sd()
    );
    assert!(training.gc() > 0.20 && training.gc() < 0.80);

    let mut total_genes = 0;
    for (name, seq) in &contigs {
        let genes = predict(seq, &training).unwrap();
        total_genes += genes.len();
        eprintln!(
            "  [single] {} ({} bp): {} genes",
            name,
            seq.len(),
            genes.len()
        );
    }
    eprintln!("[priestia single] {} total genes", total_genes);
    assert!(
        total_genes > 20,
        "Expected >20 genes from Priestia megaterium subset single mode, got {}",
        total_genes
    );
}

// ── Cross-mode consistency: meta vs single should broadly agree ──────────────

#[test]
fn test_meta_vs_single_gene_count_consistency() {
    // Use real genomes for this comparison — short sequences (< 100 kb)
    // don't have enough genes for a meaningful ratio comparison.
    let test_files: Vec<(&str, std::path::PathBuf)> = vec![
        (
            "prokka-genome",
            sibling_repo_path("prokka-rs/prokka/test/genome.fna"),
        ),
        (
            "priestia",
            sibling_repo_path("gecco-rs/data/CP157504.1.fna"),
        ),
    ];

    let mut tested = 0;
    for (label, path) in &test_files {
        if !path.exists() {
            eprintln!("Skipping {}: file not found", label);
            continue;
        }
        let contigs = load_fasta_contigs(path, TEST_SEQ_LIMIT);
        let all: Vec<u8> = contigs
            .iter()
            .flat_map(|(_, s)| s.iter().copied())
            .collect();
        let training = train(&all).unwrap();

        for (name, seq) in &contigs {
            let meta_genes = predict_meta(seq).unwrap();
            let single_genes = predict(seq, &training).unwrap();

            if single_genes.is_empty() {
                continue;
            }

            let ratio = meta_genes.len() as f64 / single_genes.len() as f64;
            eprintln!(
                "  [{}] {} meta={} single={} ratio={:.2}",
                label,
                name,
                meta_genes.len(),
                single_genes.len(),
                ratio
            );

            // On large genomes, meta and single should agree within 10%
            assert!(
                ratio > 0.8 && ratio < 1.2,
                "Gene count ratio {:.2} for {} is too divergent (meta={}, single={})",
                ratio,
                name,
                meta_genes.len(),
                single_genes.len()
            );
            tested += 1;
        }
    }
    if tested == 0 {
        eprintln!("Skipping: no real genomes available for cross-mode comparison");
    }
}

// ── Gene coordinate validity across all modes ────────────────────────────────

#[test]
fn test_gene_coordinates_within_sequence_bounds() {
    let seqs = load_sample_sequences();
    let all = concat_sequences(&seqs);
    let training = train(&all).unwrap();

    for (name, seq) in &seqs {
        // Test both modes
        for (mode, genes) in [
            ("meta", predict_meta(seq).unwrap()),
            ("single", predict(seq, &training).unwrap()),
        ] {
            for (i, g) in genes.iter().enumerate() {
                let (start, end) = if g.begin < g.end {
                    (g.begin, g.end)
                } else {
                    (g.end, g.begin)
                };
                assert!(
                    end <= seq.len(),
                    "[{}] {} gene {} end ({}) exceeds sequence length ({})",
                    mode,
                    name,
                    i,
                    end,
                    seq.len()
                );
                assert!(
                    start >= 1,
                    "[{}] {} gene {} start ({}) is less than 1",
                    mode,
                    name,
                    i,
                    start
                );
                // Gene length should be a multiple of 3 (codon-aligned)
                let len = end - start + 1;
                assert!(
                    len % 3 == 0,
                    "[{}] {} gene {} length {} is not divisible by 3 ({}..{})",
                    mode,
                    name,
                    i,
                    len,
                    start,
                    end
                );
            }
        }
    }
}

// ── Non-overlapping gene validation ──────────────────────────────────────────

#[test]
fn test_genes_are_mostly_non_overlapping() {
    let seqs = load_sample_sequences();

    for (name, seq) in &seqs {
        let genes = predict_meta(seq).unwrap();
        if genes.len() < 2 {
            continue;
        }

        // Check that genes on the same strand don't have major overlaps
        let mut forward: Vec<(usize, usize)> = Vec::new();
        let mut reverse: Vec<(usize, usize)> = Vec::new();

        for g in &genes {
            let (s, e) = if g.begin < g.end {
                (g.begin, g.end)
            } else {
                (g.end, g.begin)
            };
            match g.strand {
                Strand::Forward => forward.push((s, e)),
                Strand::Reverse => reverse.push((s, e)),
            }
        }

        for (strand_name, mut strand_genes) in [("fwd", forward), ("rev", reverse)] {
            strand_genes.sort();
            let mut big_overlaps = 0;
            for w in strand_genes.windows(2) {
                if w[1].0 < w[0].1 {
                    let overlap = w[0].1 - w[1].0;
                    // Small overlaps (< 60bp) are normal in prokaryotes
                    if overlap > 60 {
                        big_overlaps += 1;
                    }
                }
            }
            assert!(
                big_overlaps == 0,
                "[{}] {} strand has {} large overlaps (>60bp)",
                name,
                strand_name,
                big_overlaps
            );
        }
    }
}

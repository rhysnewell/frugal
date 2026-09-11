/// Strand of a predicted gene.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}

/// Start codon type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StartCodon {
    ATG,
    GTG,
    TTG,
    Edge,
}

impl std::fmt::Display for StartCodon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            StartCodon::ATG => write!(f, "ATG"),
            StartCodon::GTG => write!(f, "GTG"),
            StartCodon::TTG => write!(f, "TTG"),
            StartCodon::Edge => write!(f, "Edge"),
        }
    }
}

/// A single predicted gene.
#[derive(Debug, Clone)]
pub struct PredictedGene {
    /// 1-indexed start position in the input sequence.
    pub begin: usize,
    /// 1-indexed end position (inclusive).
    pub end: usize,
    /// Strand.
    pub strand: Strand,
    /// Start codon type.
    pub start_codon: StartCodon,
    /// Translation table used for this gene.
    pub translation_table: u8,
    /// Whether left/right ends are partial (run off sequence edges).
    pub partial: (bool, bool),
    /// RBS motif name (e.g. "AGGAG" or "None").
    pub rbs_motif: String,
    /// RBS spacer distance (e.g. "3-4bp" or "None").
    pub rbs_spacer: String,
    /// GC content of this gene.
    pub gc_content: f64,
    /// Confidence score (50.0 - 99.99).
    pub confidence: f64,
    /// Total score (coding + start).
    pub score: f64,
    /// Coding score.
    pub cscore: f64,
    /// Start score (rscore + uscore + tscore).
    pub sscore: f64,
    /// RBS score.
    pub rscore: f64,
    /// Upstream composition score.
    pub uscore: f64,
    /// Start codon type score.
    pub tscore: f64,
}

/// Configuration for gene prediction.
#[derive(Debug, Clone)]
pub struct ProdigalConfig {
    /// NCBI translation table (1-25, excluding 7,8,17-20). Default: 11.
    pub translation_table: u8,
    /// Do not allow genes to run off sequence edges.
    pub closed_ends: bool,
    /// Treat runs of N as masked regions.
    pub mask_n_runs: bool,
    /// Bypass Shine-Dalgarno trainer, force full motif scan.
    pub force_non_sd: bool,
    /// In meta mode, run the full path only for this many models, ranked by their best node
    /// score. 0 runs every model in the GC window.
    pub model_depth: usize,
}

impl Default for ProdigalConfig {
    fn default() -> Self {
        ProdigalConfig {
            translation_table: 11,
            closed_ends: false,
            mask_n_runs: false,
            force_non_sd: false,
            model_depth: 0,
        }
    }
}

/// Errors from gene prediction.
#[derive(Debug)]
pub enum ProdigalError {
    EmptySequence,
    SequenceTooShort { length: usize, min: usize },
    SequenceTooLong { length: usize, max: usize },
    InvalidTranslationTable(u8),
    Io(std::io::Error),
}

impl std::fmt::Display for ProdigalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ProdigalError::EmptySequence => write!(f, "sequence is empty"),
            ProdigalError::SequenceTooShort { length, min } => {
                write!(f, "sequence too short ({length} bp, need >= {min} bp)")
            }
            ProdigalError::SequenceTooLong { length, max } => {
                write!(f, "sequence too long ({length} bp, max {max} bp)")
            }
            ProdigalError::InvalidTranslationTable(t) => {
                write!(f, "invalid translation table: {t}")
            }
            ProdigalError::Io(e) => write!(f, "I/O error: {e}"),
        }
    }
}

impl std::error::Error for ProdigalError {}

impl From<std::io::Error> for ProdigalError {
    fn from(e: std::io::Error) -> Self {
        ProdigalError::Io(e)
    }
}

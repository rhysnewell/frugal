#![allow(
    clippy::assign_op_pattern,
    clippy::collapsible_else_if,
    clippy::collapsible_if,
    clippy::cmp_null,
    clippy::if_same_then_else,
    clippy::manual_clamp,
    clippy::manual_range_contains,
    clippy::missing_safety_doc,
    clippy::needless_late_init,
    clippy::needless_range_loop,
    clippy::needless_return,
    clippy::not_unsafe_ptr_arg_deref,
    clippy::print_with_newline,
    clippy::too_many_arguments,
    clippy::unnecessary_cast
)]
// The low-level modules intentionally preserve the original Prodigal C
// function boundaries and control flow for auditability against `Prodigal/`.

pub mod api;
pub mod bitmap;
pub mod connection;
pub mod connection_filter;
pub mod dprog;
pub mod gene;
pub mod metagenomic;
pub mod node;
pub mod output;
pub mod pipeline;
pub mod reader;
pub mod sequence;
pub mod training;
pub mod training_data;
pub mod types;

// Re-export the high-level API at crate root
pub use api::{
    predict, predict_meta, predict_meta_with, predict_with, train, train_with, MetaPredictor,
    PredictedGene, ProdigalConfig, ProdigalError, StartCodon, Strand, TrainingData,
    META_PREDICTOR_STACK_SIZE,
};

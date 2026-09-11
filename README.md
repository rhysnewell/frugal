# frugal

Fast rosella prodigal. A fork of [prodigal-rs](https://github.com/henriksson-lab/prodigal-rs)
tuned for calling genes over a metagenome assembly, where the work is many short contigs rather
than one long genome.

Lineage: [Prodigal](https://github.com/hyattpd/Prodigal) v2.6.3 (Hyatt et al., 2010), translated
to Rust by henriksson-lab as `prodigal-rs` 0.3.5, forked here.

**This is a modified version of prodigal-rs.** Every change is listed below.

## What changed

Metagenomic mode is about three times faster than 0.3.5 and the genes are byte identical.
Measured on 23,985 contigs of a CAMI II urogenital sample, 157,065 genes, 4 threads.

| change | effect |
|---|---|
| Reuse the per-sequence scratch buffer across a batch | 0.3.5 allocated and zeroed about 180 MB for every contig, whatever its length |
| Cache the Shine-Dalgarno candidates across the 50 models | the candidate set is sequence-only, so it is recorded once and folded per model. That stage went 70.8 s to 13.3 s |
| Split connection scoring on the right node's strand and type | the right node is constant for a whole pass, so it picks one of four scorers once instead of being retested on every pair |
| Build the node list only for bins inside the GC window | the node set depends on the translation table alone |
| Hoist the length factor bounds out of the coding score loop | two `pow` and two `log` per node per model |

All of the above are exact. One optional change is not:

`ProdigalConfig::model_depth` runs the full path search for only the highest ranked models
instead of every model the GC window admits, ranking them on their best node score. Default 0
runs them all, which is the exact path. Depth 3 was 1.85x faster and left 98.7 per cent of gene
coordinates unchanged.

## Carried over from prodigal-rs

* **This is an LLM-mediated translation, not the original code.** Try the original first unless
  you have a reason not to. It has more features and more bug fixing behind it.
* **Do not send bug reports to the original developers**, of Prodigal or of prodigal-rs. Use this
  repository's issues.
* The translation aims to reproduce the original's behaviour, bugs included, so results stay
  comparable across studies.
* Single genome training mode is inherited untouched. The changes here are all on the
  metagenomic path.

## Use

```toml
[dependencies]
frugal = "0.4"
```

Calling genes over a batch of contigs, which is what this fork is for:

```rust
use std::sync::Arc;
use frugal::{META_PREDICTOR_STACK_SIZE, MetaPredictor, ProdigalConfig};

let pool = Arc::new(
    rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .stack_size(META_PREDICTOR_STACK_SIZE)
        .build()?,
);
let predictor = MetaPredictor::with_config_and_thread_pool(ProdigalConfig::default(), pool)?;
let genes = predictor.predict_batch(&contigs)?;
```

`predict_batch` holds one scratch buffer and one set of models per worker thread, so feeding it
the whole assembly at once costs no more memory than feeding it one contig at a time. Feed the
long contigs first if you can; the batch runs to its slowest member.

Everything else in the prodigal-rs API is unchanged: `predict_meta`, `train`, `predict`,
`TrainingData`, and the `prodigal` compatible CLI, now installed as `frugal`.

## Licence

GPL-3.0-only, as Prodigal and prodigal-rs are. The original copyright notices are intact in
every file they were in.

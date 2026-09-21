use std::path::{Path, PathBuf};
use std::process::Command;

fn manifest(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

fn run_mode(mode: &str) -> String {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.sco");
    let status = Command::new(env!("CARGO_BIN_EXE_frugal"))
        .args([
            "-i",
            manifest("tests/data/anthus_aco.fas").to_str().unwrap(),
            "-p",
            mode,
            "-f",
            "sco",
            "-q",
            "-o",
            out.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "frugal -p {mode} exited with {status}");
    std::fs::read_to_string(&out).unwrap()
}

#[test]
fn meta_mode_output_is_unchanged() {
    let expected = std::fs::read_to_string(manifest("tests/data/anthus_aco.meta.sco")).unwrap();
    assert_eq!(run_mode("meta"), expected);
}

#[test]
fn single_mode_output_is_unchanged() {
    let expected = std::fs::read_to_string(manifest("tests/data/anthus_aco.single.sco")).unwrap();
    assert_eq!(run_mode("single"), expected);
}

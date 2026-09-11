//! Verify that Rust `#[repr(C)]` struct sizes match their C counterparts.
//!
//! This test compiles and runs a small C program that prints sizeof() for each
//! struct, then compares against std::mem::size_of::<T>().

use std::process::Command;

fn get_c_sizes() -> Vec<(String, usize)> {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let c_src = format!("{}/Prodigal", manifest_dir);
    let out_bin = format!(
        "{}/target/struct_size_check{}",
        manifest_dir,
        std::env::consts::EXE_SUFFIX
    );

    // Compile the C size-check program
    let status = Command::new("gcc")
        .args([
            "-o",
            &out_bin,
            "-xc",
            "-", // read from stdin
            &format!("-I{}", c_src),
            "-DSUPPORT_GZIP_COMPRESSED",
        ])
        .stdin(std::process::Stdio::piped())
        .spawn()
        .and_then(|mut child| {
            use std::io::Write;
            let stdin = child.stdin.as_mut().unwrap();
            stdin.write_all(
                br#"
                #include <stdio.h>
                #include "training.h"
                #include "node.h"
                #include "gene.h"
                #include "metagenomic.h"
                #include "sequence.h"
                int main() {
                    printf("_mask %zu\n", sizeof(mask));
                    printf("_training %zu\n", sizeof(struct _training));
                    printf("_motif %zu\n", sizeof(struct _motif));
                    printf("_node %zu\n", sizeof(struct _node));
                    printf("_gene %zu\n", sizeof(struct _gene));
                    printf("_metagenomic_bin %zu\n", sizeof(struct _metagenomic_bin));
                    return 0;
                }
            "#,
            )?;
            child.wait()
        })
        .expect("Failed to compile C size-check program");
    assert!(status.success(), "C size-check program failed to compile");

    let output = Command::new(&out_bin)
        .output()
        .expect("Failed to run C size-check program");
    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout).unwrap();
    stdout
        .lines()
        .map(|line| {
            let mut parts = line.split_whitespace();
            let name = parts.next().unwrap().to_string();
            let size: usize = parts.next().unwrap().parse().unwrap();
            (name, size)
        })
        .collect()
}

#[test]
fn test_struct_sizes_match_c() {
    use frugal::types::*;

    let c_sizes = get_c_sizes();

    let rust_sizes: Vec<(&str, usize)> = vec![
        ("_mask", std::mem::size_of::<Mask>()),
        ("_training", std::mem::size_of::<Training>()),
        ("_motif", std::mem::size_of::<Motif>()),
        ("_node", std::mem::size_of::<Node>()),
        ("_gene", std::mem::size_of::<Gene>()),
        ("_metagenomic_bin", std::mem::size_of::<MetagenomicBin>()),
    ];

    for (c_name, c_size) in &c_sizes {
        let rust_entry = rust_sizes
            .iter()
            .find(|(name, _)| name == c_name)
            .unwrap_or_else(|| panic!("No Rust struct found for C struct {c_name}"));

        assert_eq!(
            *c_size, rust_entry.1,
            "Size mismatch for {c_name}: C={c_size}, Rust={}",
            rust_entry.1
        );
    }

    eprintln!("All struct sizes match:");
    for (name, size) in &c_sizes {
        eprintln!("  {name}: {size} bytes");
    }
}

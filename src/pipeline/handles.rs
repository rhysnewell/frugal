use std::ffi::CString;
use std::os::raw::{c_int, c_void};

use super::PipelineConfig;
use crate::output::{close_handle, create_file, stdout_handle};
use crate::reader::{seq_reader_close, seq_reader_open, seq_reader_open_stdin};

pub(super) struct Handles {
    pub input: *mut c_void,
    pub output: c_int,
    pub start: c_int,
    pub trans: c_int,
    pub nuc: c_int,
    pub stdout: c_int,
}

impl Handles {
    pub unsafe fn open(config: &PipelineConfig) -> Result<Self, i32> {
        let input = match &config.input_file {
            Some(path) => {
                let path_c = CString::new(path.as_str()).unwrap();
                let input = seq_reader_open(path_c.as_ptr());
                if input.is_null() {
                    eprintln!("\nError: can't open input file {}.\n", path);
                    return Err(5);
                }
                input
            }
            None => {
                let input = seq_reader_open_stdin();
                if input.is_null() {
                    eprintln!("\nError: can't open stdin.\n");
                    return Err(5);
                }
                input
            }
        };
        let stdout = stdout_handle();
        let open = |path: &Option<String>, what: &str, code: i32| match path {
            None => Ok(stdout),
            Some(path) => create_file(path).map_err(|_| {
                eprintln!("\nError: can't open {} {}.\n", what, path);
                code
            }),
        };
        Ok(Self {
            input,
            output: open(&config.output_file, "output file", 6)?,
            start: open(&config.start_file, "start file", 7)?,
            trans: open(&config.trans_file, "translation file", 8)?,
            nuc: open(&config.nuc_file, "gene nucleotide file", 16)?,
            stdout,
        })
    }

    pub unsafe fn close(&self) {
        seq_reader_close(self.input);
        for handle in [self.output, self.start, self.trans, self.nuc] {
            if handle != self.stdout {
                close_handle(handle);
            }
        }
    }
}

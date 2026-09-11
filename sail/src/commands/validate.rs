use std::path::PathBuf;

use anyhow::{Result, bail};
use clap::Args;
use libsail::collection::Indexable;
use libsail::format::Format;
use libsail::seq::fasta::IndexedFasta;
use libsail::seq::p7hmm::IndexedHmm;
use libsail::seq::stockholm::IndexedStockholm;

use crate::cli::FormatArg;
use crate::input::format_of_path;

#[derive(Args)]
pub struct ValidateArgs {
    /// files to check
    pub input: Vec<PathBuf>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,
}

impl ValidateArgs {
    pub fn run(self) -> Result<()> {
        if self.input.is_empty() {
            bail!("validate needs a path: it reports where in a file a record failed");
        }

        for path in &self.input {
            // a path rather than buffered bytes, and the disk
            // backends rather than the in-memory ones: only
            // those parse one record at a time, so only they
            // can say which record failed. Fasta::open parses
            // everything before it returns anything
            if path.as_os_str() == "-" {
                bail!("validate needs a path, not stdin: it reports a record's byte offset");
            }

            let count = match format_of_path(path, self.format)? {
                Format::Fasta => {
                    let c = IndexedFasta::open(path)?;
                    c.validate()?;
                    c.len()
                }
                Format::Stockholm => {
                    let c = IndexedStockholm::open(path)?;
                    c.validate()?;
                    c.len()
                }
                Format::Hmm => {
                    let c = IndexedHmm::open(path)?;
                    c.validate()?;
                    c.len()
                }
            };

            println!("{}\t{count} records\tok", path.display());
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::*;

    fn fixture(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../fixtures")
            .join(name)
    }

    fn validate(paths: &[PathBuf]) -> Result<()> {
        ValidateArgs {
            input: paths.to_vec(),
            format: None,
        }
        .run()
    }

    #[test]
    fn every_fixture_is_real_tool_output_and_validates() {
        for name in ["proteins.fa", "families.sto", "models.hmm"] {
            assert!(validate(&[fixture(name)]).is_ok(), "{name}");
        }
    }

    #[test]
    fn a_record_that_does_not_parse_is_reported_with_its_index() {
        // the second alignment is missing its header line,
        // which the index still frames as a record
        let path = std::env::temp_dir().join(format!("sail-val-{}.sto", std::process::id()));
        std::fs::write(&path, b"# STOCKHOLM 1.0\nseq1 AC\n//\nseq1 AC\n//\n").unwrap();

        let error = format!("{:#}", validate(std::slice::from_ref(&path)).unwrap_err());
        assert!(error.contains("record 1"), "{error}");

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn stdin_is_refused_because_a_byte_offset_needs_a_file() {
        assert!(validate(&[PathBuf::from("-")]).is_err());
    }
}

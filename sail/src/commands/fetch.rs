use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Args;
use libsail::collection::{Indexable, Iterable};

use crate::cli::FormatArg;
use crate::input::{Inputs, dispatch};
use crate::output::{emit, writer};

#[derive(Args)]
pub struct FetchArgs {
    /// file to read
    pub input: PathBuf,

    /// names to extract
    pub names: Vec<String>,

    /// read the names from this file, one per line
    #[arg(short = 'f', long)]
    pub names_file: Option<PathBuf>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// skip a name the input does not hold, rather than failing
    #[arg(long)]
    pub skip_missing: bool,
}

impl FetchArgs {
    pub fn run(self) -> Result<()> {
        let wanted = self.wanted()?;
        if wanted.is_empty() {
            bail!("fetch needs at least one name, on the command line or through -f");
        }

        let inputs = Inputs::read(std::slice::from_ref(&self.input), self.format)?;
        let entry = &inputs.entries()[0];
        let out = writer(self.output.as_deref())?;
        let skip_missing = self.skip_missing;

        dispatch!(inputs.format(), entry, |collection, _size, name| {
            let mut at: HashMap<&[u8], usize> = HashMap::with_capacity(collection.len());
            for (n, record) in collection.iter().enumerate() {
                if let Some(found) = name(record) {
                    // the first of a repeated name, because a
                    // duplicate identifier means a malformed
                    // file and the last is not the one the
                    // file leads with
                    at.entry(found).or_insert(n);
                }
            }

            let mut picked = Vec::with_capacity(wanted.len());
            for want in &wanted {
                match at.get(want.as_slice()) {
                    Some(&n) => picked.push(
                        collection
                            .cloned(n)
                            .expect("a position below len() is a record"),
                    ),
                    None if skip_missing => {}
                    None => bail!(
                        "{} holds no record named {:?}",
                        entry.name(),
                        String::from_utf8_lossy(want)
                    ),
                }
            }

            emit(picked, out)
        })
    }

    /// Every name asked for, in the order asked, from the command line and
    /// then from `-f`.
    //
    // bytes rather than String: a FASTA name need not be
    // utf-8, so a name file is read as bytes and split on
    // newlines rather than decoded
    fn wanted(&self) -> Result<Vec<Vec<u8>>> {
        let mut wanted: Vec<Vec<u8>> = self.names.iter().map(|n| n.as_bytes().to_vec()).collect();

        let Some(path) = &self.names_file else {
            return Ok(wanted);
        };

        let bytes = std::fs::read(path)
            .with_context(|| format!("reading the names in {}", path.display()))?;

        wanted.extend(
            bytes
                .split(|&b| b == b'\n')
                .map(|line| line.strip_suffix(b"\r").unwrap_or(line))
                .filter(|line| !line.is_empty())
                .map(<[u8]>::to_vec),
        );

        Ok(wanted)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args() -> FetchArgs {
        FetchArgs {
            input: PathBuf::new(),
            names: Vec::new(),
            names_file: None,
            format: None,
            output: None,
            skip_missing: false,
        }
    }

    #[test]
    fn names_come_from_the_command_line_and_then_from_the_file() {
        let path = std::env::temp_dir().join(format!("sail-fw-{}.txt", std::process::id()));
        std::fs::write(&path, b"PDZ\nSH3_1\n").unwrap();

        let fetch = FetchArgs {
            names: vec!["first".to_string()],
            names_file: Some(path.clone()),
            ..args()
        };

        assert_eq!(
            fetch.wanted().unwrap(),
            [&b"first"[..], b"PDZ", b"SH3_1"].map(<[u8]>::to_vec)
        );

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn a_name_file_drops_blank_lines_and_a_trailing_newline() {
        // a file written by another tool ends in \n, and a
        // trailing empty name would be looked up and missed
        let path = std::env::temp_dir().join(format!("sail-fb-{}.txt", std::process::id()));
        std::fs::write(&path, b"PDZ\n\r\n\nSH3_1\r\n").unwrap();

        let fetch = FetchArgs {
            names_file: Some(path.clone()),
            ..args()
        };

        assert_eq!(
            fetch.wanted().unwrap(),
            [&b"PDZ"[..], b"SH3_1"].map(<[u8]>::to_vec)
        );

        std::fs::remove_file(path).ok();
    }
}

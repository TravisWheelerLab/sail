use std::collections::HashSet;
use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::Iterable;

use crate::cli::FormatArg;
use crate::input::{Inputs, dispatch};
use crate::output::{emit, writer};

#[derive(Args)]
pub struct DedupArgs {
    /// files to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: Vec<PathBuf>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// drop every record with no name, rather than keeping them all
    #[arg(long)]
    pub drop_unnamed: bool,
}

impl DedupArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::read(&self.input, self.format)?;
        let mut out = writer(self.output.as_deref())?;
        let drop_unnamed = self.drop_unnamed;

        // across every input, not per file: two files
        // concatenated are one stream, and a name repeated
        // between them is still a repeat
        let mut seen: HashSet<Vec<u8>> = HashSet::new();

        for entry in inputs.entries() {
            dispatch!(inputs.format(), entry, |collection, _size, name| {
                let mut kept = Vec::new();

                for record in collection.iter() {
                    match name(record) {
                        Some(found) => {
                            // insert reports whether it was
                            // new, so the first of a repeat is
                            // the one kept
                            if seen.insert(found.to_vec()) {
                                kept.push(record);
                            }
                        }

                        // an unnamed record repeats nothing:
                        // there is no name for it to share
                        None if !drop_unnamed => kept.push(record),
                        None => {}
                    }
                }

                emit(kept, &mut out)?;
            });
        }

        Ok(())
    }
}

use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::Iterable;

use crate::cli::FormatArg;
use crate::input::{Inputs, dispatch};
use crate::output::{emit, writer};

#[derive(Args)]
pub struct CatArgs {
    /// files to concatenate, or - for stdin
    #[arg(default_value = "-")]
    pub input: Vec<PathBuf>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,
}

impl CatArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::read(&self.input, self.format)?;
        let mut out = writer(self.output.as_deref())?;

        // parsed and re-emitted rather than copied through:
        // concatenating the bytes of two Stockholm files
        // gives one stream with two headers in it, which is
        // not an alignment file. going through the parser is
        // what makes the result readable again
        for entry in inputs.entries() {
            dispatch!(inputs.format(), entry, |collection| {
                emit(collection.iter(), &mut out)?;
            });
        }

        Ok(())
    }
}

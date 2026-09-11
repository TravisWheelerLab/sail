use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::Iterable;

use crate::cli::FormatArg;
use crate::input::{Inputs, dispatch};
use crate::output::{emit, writer};

#[derive(Args)]
pub struct ReformatArgs {
    /// file to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: PathBuf,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,
}

impl ReformatArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::read(std::slice::from_ref(&self.input), self.format)?;
        let entry = &inputs.entries()[0];
        let out = writer(self.output.as_deref())?;

        // every record through the parser and back out through
        // the writer, so what comes out is what this crate
        // emits rather than whatever the producing tool did:
        // one line width, one field spacing, one block layout
        dispatch!(inputs.format(), entry, |collection| {
            emit(collection.iter(), out)
        })
    }
}

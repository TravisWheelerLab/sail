use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::{Indexable, Iterable};

use crate::cli::FormatArg;
use crate::input::{Inputs, dispatch};
use crate::output::{emit, writer};

#[derive(Args)]
pub struct TailArgs {
    /// file to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: PathBuf,

    /// how many records to keep
    #[arg(short = 'n', long, default_value_t = 10)]
    pub count: usize,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,
}

impl TailArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::read(std::slice::from_ref(&self.input), self.format)?;
        let entry = &inputs.entries()[0];
        let out = writer(self.output.as_deref())?;
        let count = self.count;

        dispatch!(inputs.format(), entry, |collection| {
            // saturating, so asking for more records than the
            // file holds skips none rather than wrapping to a
            // huge skip and emitting nothing
            let skip = collection.len().saturating_sub(count);

            emit(collection.skip(skip).iter(), out)
        })
    }
}

use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::{Indexable, Iterable};

use crate::cli::FormatArg;
use crate::input::{Inputs, dispatch};
use crate::output::{emit, writer};

#[derive(Args)]
pub struct ShuffleArgs {
    /// file to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: PathBuf,

    /// seed the shuffle, so the same input gives the same order
    #[arg(long)]
    pub seed: Option<u64>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,
}

impl ShuffleArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::read(std::slice::from_ref(&self.input), self.format)?;
        let entry = &inputs.entries()[0];
        let out = writer(self.output.as_deref())?;

        // an unseeded run picks a seed rather than taking a
        // second code path: libsail has one permutation and
        // it is the reproducible one
        let seed = self.seed.unwrap_or_else(rand::random);

        dispatch!(inputs.format(), entry, |collection, _size, _name, write| {
            emit(collection.permute(seed).iter(), write, out)
        })
    }
}

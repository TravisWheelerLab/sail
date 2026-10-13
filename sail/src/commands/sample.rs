use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::{Indexable, Iterable};

use crate::cli::FormatArg;
use crate::commands::{default_seed, draw_size};
use crate::input::{Inputs, dispatch};
use crate::output::{emit, writer};

#[derive(Args)]
pub struct SampleArgs {
    /// file to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: PathBuf,

    /// how many records to draw
    #[arg(short = 'n', long, conflicts_with = "fraction")]
    pub count: Option<usize>,

    /// what share of the records to draw, between 0 and 1
    #[arg(short = 'p', long)]
    pub fraction: Option<f64>,

    /// seed the draw, so the same input gives the same sample
    #[arg(long)]
    pub seed: Option<u64>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,
}

impl SampleArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::read(std::slice::from_ref(&self.input), self.format)?;
        let entry = &inputs.entries()[0];
        let out = writer(self.output.as_deref())?;
        let (count, fraction) = (self.count, self.fraction);

        // an unseeded run picks a seed rather than taking a
        // second code path: a draw that ends up in a paper
        // has to reproduce, and libsail has one draw
        let seed = self.seed.unwrap_or_else(default_seed);

        dispatch!(inputs.format(), entry, |collection, _size, _name, write| {
            let n = draw_size(count, fraction, collection.len())?;

            emit(collection.sample(n, seed).iter(), write, out)
        })
    }
}

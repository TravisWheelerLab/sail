use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::{Indexable, Iterable};
use rand::SeedableRng;
use rand::rngs::StdRng;

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
        let seed = self.seed;

        dispatch!(inputs.format(), entry, |collection| {
            let mut rng = match seed {
                Some(seed) => StdRng::seed_from_u64(seed),
                None => StdRng::from_rng(&mut rand::rng()),
            };

            emit(collection.permute_with(&mut rng).iter(), out)
        })
    }
}

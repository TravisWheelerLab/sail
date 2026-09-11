use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::Iterable;

use crate::cli::{Axis, FormatArg};
use crate::input::{Inputs, axis_for, dispatch};
use crate::output::{emit, writer};

#[derive(Args)]
pub struct FilterArgs {
    /// files to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: Vec<PathBuf>,

    /// keep records of at least this size
    #[arg(long)]
    pub min: Option<usize>,

    /// keep records of at most this size
    #[arg(long)]
    pub max: Option<usize>,

    /// for an alignment, measure rows or columns [default: depth]
    #[arg(long, value_enum)]
    pub by: Option<Axis>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// write the records that did not survive here
    #[arg(long)]
    pub rest: Option<PathBuf>,
}

impl FilterArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::read(&self.input, self.format)?;
        let axis = axis_for(inputs.format(), self.by)?;

        let (min, max) = (self.min.unwrap_or(0), self.max.unwrap_or(usize::MAX));
        let mut kept = writer(self.output.as_deref())?;
        let mut rest = match self.rest.as_deref() {
            Some(path) => Some(writer(Some(path))?),
            None => None,
        };

        for entry in inputs.entries() {
            dispatch!(inputs.format(), entry, axis, |collection, size, _name| {
                // one pass with a partition rather than two
                // filters, so --rest cannot disagree with the
                // kept side about where the boundary is
                let (keep, drop): (Vec<_>, Vec<_>) = collection
                    .iter()
                    .partition(|r| (min..=max).contains(&size(r)));

                emit(keep, &mut kept)?;
                if let Some(rest) = &mut rest {
                    emit(drop, rest)?;
                }
            });
        }

        Ok(())
    }
}

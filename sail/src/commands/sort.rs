use std::path::PathBuf;

use anyhow::Result;
use clap::{Args, ValueEnum};
use libsail::collection::{Indexable, Iterable};

use crate::cli::{Axis, FormatArg};
use crate::input::{Inputs, axis_for, dispatch};
use crate::output::{emit, writer};

/// What `sort` orders records by.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Key {
    Name,
    Size,
}

#[derive(Args)]
pub struct SortArgs {
    /// file to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: PathBuf,

    /// order by name or by size
    #[arg(long, value_enum, default_value = "name")]
    pub key: Key,

    /// largest or last first
    #[arg(short, long)]
    pub reverse: bool,

    /// for an alignment, measure rows or columns [default: depth]
    #[arg(long, value_enum)]
    pub by: Option<Axis>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,
}

impl SortArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::read(std::slice::from_ref(&self.input), self.format)?;
        let entry = &inputs.entries()[0];
        let axis = axis_for(inputs.format(), self.by)?;
        let out = writer(self.output.as_deref())?;
        let (key, reverse) = (self.key, self.reverse);

        dispatch!(
            inputs.format(),
            entry,
            axis,
            |collection, size, name, write| {
                let sorted = collection.sort_by(|a, b| {
                    let ordering = match key {
                        // a record with no name orders before
                        // every named one, because None < Some
                        // -- and the alternative, dropping it,
                        // would make sort lose records
                        Key::Name => name(a).cmp(&name(b)),
                        Key::Size => size(a).cmp(&size(b)),
                    };

                    if reverse {
                        ordering.reverse()
                    } else {
                        ordering
                    }
                });

                emit(sorted.iter(), write, out)
            }
        )
    }
}

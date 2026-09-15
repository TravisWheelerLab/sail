use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use std::io::Write;

use libsail::collection::Iterable;
use libsail::index::Reader;

use crate::cli::{Axis, FormatArg, ReadArgs};
use crate::input::{
    Backend, Inputs, Needs, axis_for, dispatch, indexed, indexed_path, size_framed,
};
use crate::output::{emit, put, writer};

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

    /// re-wrap the records written out, rather than copying their bytes through
    #[arg(long)]
    pub rewrap: bool,

    #[command(flatten)]
    pub read: ReadArgs,
}

impl FilterArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::plan(&self.input, self.format, self.read, Needs::Pass)?;
        let axis = axis_for(inputs.format(), self.by)?;

        let (min, max) = (self.min.unwrap_or(0), self.max.unwrap_or(usize::MAX));
        let mut kept = writer(self.output.as_deref())?;
        let mut rest = match self.rest.as_deref() {
            Some(path) => Some(writer(Some(path))?),
            None => None,
        };

        let format = inputs.format();
        let rewrap = self.rewrap;

        // one predicate and one write per record, so --rest
        // cannot disagree with the kept side about where the
        // boundary is
        let sort_one = |record: &[u8],
                        kept: &mut dyn Write,
                        rest: &mut Option<Box<dyn Write>>|
         -> Result<()> {
            let size = size_framed(format, record, axis)?;

            match ((min..=max).contains(&size), rest) {
                (true, _) => put(format, record, rewrap, kept)?,
                (false, Some(rest)) => put(format, record, rewrap, rest)?,
                (false, None) => {}
            }

            Ok(())
        };

        for entry in inputs.entries() {
            match inputs.backend() {
                Backend::Stream => {
                    let mut reader = Reader::new(entry.reader()?, inputs.format());

                    while reader.advance()? {
                        sort_one(reader.record(), &mut kept, &mut rest)?;
                    }
                }
                Backend::Indexed => {
                    let path = indexed_path(entry)?;

                    indexed!(inputs.format(), path, |collection| {
                        for n in 0..collection.index().len() {
                            let record = collection.record(n)?.expect("a counted record");

                            sort_one(&record, &mut kept, &mut rest)?;
                        }
                    });
                }
                Backend::Memory => {
                    dispatch!(
                        inputs.format(),
                        entry,
                        axis,
                        |collection, size, _name, write| {
                            // one pass with a partition rather
                            // than two filters, so --rest
                            // cannot disagree with the kept
                            // side about where the boundary is
                            let (keep, drop): (Vec<_>, Vec<_>) = collection
                                .iter()
                                .partition(|r| (min..=max).contains(&size(r)));

                            emit(keep, write, &mut kept)?;
                            if let Some(rest) = &mut rest {
                                emit(drop, write, rest)?;
                            }
                        }
                    );
                }
            }
        }

        Ok(())
    }
}

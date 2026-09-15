use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use std::io::Write;

use libsail::collection::{Indexable, Iterable};
use libsail::index::Reader;

use crate::cli::{FormatArg, ReadArgs};
use crate::input::{Backend, Inputs, Needs, dispatch, indexed, indexed_path};
use crate::output::{emit, write_framed, writer};

#[derive(Args)]
pub struct HeadArgs {
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

    #[command(flatten)]
    pub read: ReadArgs,
}

impl HeadArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::plan(
            std::slice::from_ref(&self.input),
            self.format,
            self.read,
            Needs::Pass,
        )?;
        let entry = &inputs.entries()[0];
        let out = writer(self.output.as_deref())?;
        let count = self.count;

        match inputs.backend() {
            // only the first `count` records are read: the
            // bytes past them are never touched
            Backend::Indexed => {
                let path = indexed_path(entry)?;
                let mut out = out;

                indexed!(inputs.format(), path, |collection| {
                    for n in 0..collection.index().len().min(count) {
                        let record = collection.record(n)?.expect("a counted record");

                        write_framed(inputs.format(), &record, &mut out)?;
                    }
                });
                out.flush()?;

                Ok(())
            }
            Backend::Stream => {
                let mut reader = Reader::new(entry.reader()?, inputs.format());
                let mut out = out;
                let mut n = 0;

                // clamping falls out of the loop: a file
                // shorter than `count` ends the pass, which
                // is the whole file rather than an error --
                // the same answer head(1) gives
                while n < count && reader.advance()? {
                    write_framed(inputs.format(), reader.record(), &mut out)?;
                    n += 1;
                }
                out.flush()?;

                Ok(())
            }
            Backend::Memory => {
                dispatch!(inputs.format(), entry, |collection, _size, _name, write| {
                    emit(collection.take(count).iter(), write, out)
                })
            }
        }
    }
}

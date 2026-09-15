use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use std::io::Write;

use libsail::collection::Iterable;
use libsail::index::Reader;

use crate::cli::{FormatArg, ReadArgs};
use crate::input::{Backend, Inputs, Needs, dispatch, indexed, indexed_path};
use crate::output::{emit, write_framed, writer};

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

    #[command(flatten)]
    pub read: ReadArgs,
}

impl ReformatArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::plan(
            std::slice::from_ref(&self.input),
            self.format,
            self.read,
            Needs::Pass,
        )?;
        let entry = &inputs.entries()[0];
        let out = writer(self.output.as_deref())?;

        // every record through the parser and back out through
        // the writer, so what comes out is what this crate
        // emits rather than whatever the producing tool did:
        // one line width, one field spacing, one block layout
        match inputs.backend() {
            Backend::Stream => {
                let mut reader = Reader::new(entry.reader()?, inputs.format());
                let mut out = out;

                while reader.advance()? {
                    write_framed(inputs.format(), reader.record(), &mut out)?;
                }
                out.flush()?;

                Ok(())
            }
            Backend::Indexed => {
                let path = indexed_path(entry)?;
                let mut out = out;

                indexed!(inputs.format(), path, |collection| {
                    for n in 0..collection.index().len() {
                        let record = collection.record(n)?.expect("a counted record");

                        write_framed(inputs.format(), &record, &mut out)?;
                    }
                });
                out.flush()?;

                Ok(())
            }
            Backend::Memory => {
                dispatch!(inputs.format(), entry, |collection, _size, _name, write| {
                    emit(collection.iter(), write, out)
                })
            }
        }
    }
}

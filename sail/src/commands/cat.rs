use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::Iterable;
use libsail::index::Reader;

use crate::cli::{FormatArg, ReadArgs};
use crate::input::{Backend, Inputs, Needs, dispatch, indexed, indexed_path};
use crate::output::{emit, write_framed, writer};

#[derive(Args)]
pub struct CatArgs {
    /// files to concatenate, or - for stdin
    #[arg(default_value = "-")]
    pub input: Vec<PathBuf>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    #[command(flatten)]
    pub read: ReadArgs,
}

impl CatArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::plan(&self.input, self.format, self.read, Needs::Pass)?;
        let mut out = writer(self.output.as_deref())?;

        // parsed and re-emitted rather than copied through:
        // concatenating the bytes of two Stockholm files
        // gives one stream with two headers in it, which is
        // not an alignment file
        for entry in inputs.entries() {
            match inputs.backend() {
                Backend::Stream => {
                    let mut reader = Reader::new(entry.reader()?, inputs.format());

                    while reader.advance()? {
                        write_framed(inputs.format(), reader.record(), &mut out)?;
                    }
                }
                Backend::Indexed => {
                    let path = indexed_path(entry)?;

                    indexed!(inputs.format(), path, |collection| {
                        for n in 0..collection.index().len() {
                            let record = collection.record(n)?.expect("a counted record");

                            write_framed(inputs.format(), &record, &mut out)?;
                        }
                    });
                }
                Backend::Memory => {
                    dispatch!(inputs.format(), entry, |collection, _size, _name, write| {
                        emit(collection.iter(), write, &mut out)?;
                    });
                }
            }
        }

        Ok(())
    }
}

use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use std::collections::VecDeque;
use std::io::Write;

use libsail::collection::{Indexable, Iterable};
use libsail::index::Reader;

use crate::cli::{FormatArg, ReadArgs};
use crate::input::{Backend, Inputs, Needs, dispatch};
use crate::output::{emit, write_framed, writer};

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

    #[command(flatten)]
    pub read: ReadArgs,
}

impl TailArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::plan(
            std::slice::from_ref(&self.input),
            self.format,
            self.read,
            Needs::Count,
        )?;
        let entry = &inputs.entries()[0];
        let out = writer(self.output.as_deref())?;
        let count = self.count;

        match inputs.backend() {
            // the end of a stream is only known by reaching
            // it, so the last n records are carried along in
            // a ring -- O(n asked for) rather than O(file)
            Backend::Stream | Backend::Indexed => {
                let mut reader = Reader::new(entry.reader()?, inputs.format());
                let mut last: VecDeque<Vec<u8>> = VecDeque::new();
                let mut out = out;

                while reader.advance()? {
                    if count == 0 {
                        continue;
                    }

                    if last.len() == count {
                        last.pop_front();
                    }
                    last.push_back(reader.record().to_vec());
                }

                for record in &last {
                    write_framed(inputs.format(), record, &mut out)?;
                }
                out.flush()?;

                Ok(())
            }
            Backend::Memory => {
                dispatch!(inputs.format(), entry, |collection, _size, _name, write| {
                    // saturating, so asking for more records
                    // than the file holds skips none rather
                    // than wrapping to a huge skip
                    let skip = collection.len().saturating_sub(count);

                    emit(collection.skip(skip).iter(), write, out)
                })
            }
        }
    }
}

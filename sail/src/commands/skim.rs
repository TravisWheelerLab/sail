use std::io::Write;
use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::Indexable;
use libsail::format::Format;
use libsail::index::Reader;
use libsail::select::Selection;

use crate::cli::{FormatArg, ReadArgs};
use crate::commands::draw_size;
use crate::input::{Backend, Input, Inputs, Needs, dispatch, indexed, indexed_path};
use crate::output::{emit, write_framed, writer};

#[derive(Args)]
pub struct SkimArgs {
    /// file to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: PathBuf,

    /// how many records to draw
    #[arg(short = 'n', long, conflicts_with = "fraction")]
    pub count: Option<usize>,

    /// what share of the records to draw, between 0 and 1
    #[arg(short = 'p', long)]
    pub fraction: Option<f64>,

    /// seed the draw, so the same input gives the same records
    #[arg(long)]
    pub seed: Option<u64>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    #[command(flatten)]
    pub read: ReadArgs,
}

impl SkimArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::plan(
            std::slice::from_ref(&self.input),
            self.format,
            self.read,
            Needs::Sweep,
        )?;
        let entry = &inputs.entries()[0];
        let out = writer(self.output.as_deref())?;
        let (count, fraction) = (self.count, self.fraction);

        // an unseeded run picks a seed rather than taking a
        // second code path: a draw that ends up in a paper
        // has to reproduce, and libsail has one draw
        let seed = self.seed.unwrap_or_else(rand::random);

        match inputs.backend() {
            Backend::Stream => {
                let held = total(entry, inputs.format())?;
                let n = draw_size(count, fraction, held)?;

                let mut drawn = Selection::new(held, n, seed).peekable();
                let mut reader = Reader::new(entry.reader()?, inputs.format());
                let mut at = 0;
                let mut out = out;

                while reader.advance()? {
                    // peek rather than a running skip count:
                    // the selection is already ascending, so
                    // the next index it names is the next one
                    // to keep
                    if drawn.peek() == Some(&at) {
                        drawn.next();
                        write_framed(inputs.format(), reader.record(), &mut out)?;
                    }

                    at += 1;
                }

                out.flush()?;

                Ok(())
            }

            Backend::Indexed => {
                let path = indexed_path(entry)?;
                let mut out = out;

                indexed!(inputs.format(), path, |collection| {
                    let n = draw_size(count, fraction, collection.len())?;

                    for at in Selection::new(collection.len(), n, seed) {
                        let record = collection.record(at)?.expect("a counted record");

                        write_framed(inputs.format(), &record, &mut out)?;
                    }
                });
                out.flush()?;

                Ok(())
            }

            Backend::Memory => {
                dispatch!(inputs.format(), entry, |collection, _size, _name, write| {
                    let n = draw_size(count, fraction, collection.len())?;

                    emit(collection.sample_in_order(n, seed), write, out)
                })
            }
        }
    }
}

// ---

/// How many records `entry` holds, off the index beside it when there is one
/// and by a pass over the records when there is not.
//
// the index file makes this free: count_in reads a 35-byte
// header where the pass reads the file. `sail index` first
// is what turns two passes into one
fn total(entry: &Input, format: Format) -> Result<usize> {
    if let Some(path) = entry.path() {
        let at = libsail::index::path_for(path);

        if libsail::index::is_current(&at, path)? {
            return Ok(libsail::index::count_in(&at)?);
        }
    }

    let mut reader = Reader::new(entry.reader()?, format);
    let mut held = 0;

    while reader.advance()? {
        held += 1;
    }

    Ok(held)
}

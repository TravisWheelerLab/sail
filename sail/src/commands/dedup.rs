use std::collections::HashSet;
use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::Iterable;
use libsail::index::Reader;

use crate::cli::{FormatArg, ReadArgs};
use crate::input::{Backend, Inputs, Needs, dispatch, indexed, indexed_path};
use crate::output::{emit, put, writer};

#[derive(Args)]
pub struct DedupArgs {
    /// files to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: Vec<PathBuf>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// drop every record with no name, rather than keeping them all
    #[arg(long)]
    pub drop_unnamed: bool,

    /// re-wrap the records written out, rather than copying their bytes through
    #[arg(long)]
    pub rewrap: bool,

    #[command(flatten)]
    pub read: ReadArgs,
}

impl DedupArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::plan(&self.input, self.format, self.read, Needs::Pass)?;
        let mut out = writer(self.output.as_deref())?;
        let drop_unnamed = self.drop_unnamed;

        // across every input, not per file: two files
        // concatenated are one stream, and a name repeated
        // between them is still a repeat
        let mut seen: HashSet<Vec<u8>> = HashSet::new();

        let format = inputs.format();

        for entry in inputs.entries() {
            match inputs.backend() {
                // only the names are held: each record is
                // written as the pass reaches it
                Backend::Stream => {
                    let mut reader = Reader::new(entry.reader()?, inputs.format());

                    while reader.advance()? {
                        let record = reader.record();

                        if first_of_its_name(&mut seen, format, record, drop_unnamed) {
                            put(format, record, self.rewrap, &mut out)?;
                        }
                    }
                }
                Backend::Indexed => {
                    let path = indexed_path(entry)?;

                    indexed!(format, path, |collection| {
                        for n in 0..collection.index().len() {
                            let record = collection.record(n)?.expect("a counted record");

                            if first_of_its_name(&mut seen, format, &record, drop_unnamed) {
                                put(format, &record, self.rewrap, &mut out)?;
                            }
                        }
                    });
                }
                Backend::Memory => {
                    dispatch!(inputs.format(), entry, |collection, _size, name, write| {
                        let mut kept = Vec::new();

                        for record in collection.iter() {
                            match name(record) {
                                Some(found) => {
                                    if seen.insert(found.to_vec()) {
                                        kept.push(record);
                                    }
                                }
                                None if !drop_unnamed => kept.push(record),
                                None => {}
                            }
                        }

                        emit(kept, write, &mut out)?;
                    });
                }
            }
        }

        Ok(())
    }
}

/// Whether this record is the first to carry its name.
fn first_of_its_name(
    seen: &mut HashSet<Vec<u8>>,
    format: libsail::format::Format,
    record: &[u8],
    drop_unnamed: bool,
) -> bool {
    match libsail::seq::name_of(format, record) {
        // insert reports whether the name was new, so the
        // first of a repeat is the one kept
        Some(found) => seen.insert(found.to_vec()),

        // an unnamed record repeats nothing: there is no
        // name for it to share
        None => !drop_unnamed,
    }
}

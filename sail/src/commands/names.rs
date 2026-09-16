use std::io::Write;
use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::Iterable;
use libsail::index::{Index, Reader};

use crate::cli::{FormatArg, ReadArgs};
use crate::input::{Backend, Inputs, Needs, dispatch};
use crate::output::writer;

#[derive(Args)]
pub struct NamesArgs {
    /// files to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: Vec<PathBuf>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// print a blank line for a record with no name, rather than skipping it
    #[arg(long)]
    pub keep_unnamed: bool,

    #[command(flatten)]
    pub read: ReadArgs,
}

impl NamesArgs {
    pub fn run(self) -> Result<()> {
        let mut out = writer(self.output.as_deref())?;
        self.write_names(&mut out)?;
        out.flush()?;

        Ok(())
    }

    /// One name per record, in file order.
    fn write_names(&self, out: &mut dyn Write) -> Result<()> {
        let inputs = Inputs::plan(&self.input, self.format, self.read, Needs::Pass)?;

        for entry in inputs.entries() {
            // written as they are found rather than
            // collected: a name is a handful of bytes, but
            // there is one per record, and holding them all
            // is holding something the size of the input
            let push = |found: Option<&[u8]>, out: &mut dyn Write| -> Result<()> {
                match found {
                    Some(found) => {
                        // bytes rather than String: a FASTA
                        // name is not required to be utf-8,
                        // and decoding here would make an
                        // ordinary latin-1 header unprintable
                        out.write_all(found)?;
                        out.write_all(b"\n")?;
                    }

                    // an alignment need not carry a #=GF ID,
                    // and an empty line for it would silently
                    // shift every name after it out of step
                    // with its record
                    None if self.keep_unnamed => out.write_all(b"\n")?,
                    None => {}
                }

                Ok(())
            };

            match inputs.backend() {
                // the name is a slice of the reader's buffer,
                // so nothing but the name itself is built --
                // no record, and none of its residues
                Backend::Stream => {
                    let mut reader = Reader::new(entry.reader()?, inputs.format());

                    while reader.advance()? {
                        push(libsail::seq::name_of(inputs.format(), reader.record()), out)?;
                    }
                }
                Backend::Indexed => {
                    let index = Index::build_named(entry.reader()?, inputs.format())?;

                    for n in 0..index.len() {
                        push(index.name(n)?.as_deref(), out)?;
                    }
                }
                Backend::Memory => {
                    dispatch!(inputs.format(), entry, |collection, _size, name, _write| {
                        for record in collection.iter() {
                            push(name(record), out)?;
                        }
                    });
                }
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::*;

    fn fixture(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../fixtures")
            .join(name)
    }

    fn names(paths: &[PathBuf], keep_unnamed: bool) -> Result<Vec<String>> {
        let mut out = Vec::new();

        NamesArgs {
            input: paths.to_vec(),
            format: None,
            output: None,
            keep_unnamed,
            read: ReadArgs::default(),
        }
        .write_names(&mut out)?;

        // every name ends in a newline, so the split leaves
        // one empty piece at the end that is not a name. an
        // empty piece anywhere else is --keep-unnamed's
        // placeholder and has to survive
        let mut lines: Vec<String> = out
            .split(|&b| b == b'\n')
            .map(|n| String::from_utf8(n.to_vec()).unwrap())
            .collect();
        lines.pop();

        Ok(lines)
    }

    #[test]
    fn a_fasta_name_stops_at_the_first_space_and_drops_the_description() {
        // every one of the five headers carries a UniProt
        // description after the identifier
        let found = names(&[fixture("proteins.fa")], false).unwrap();

        assert_eq!(found.len(), 5);
        assert!(found.iter().all(|n| !n.contains(' ')), "{found:?}");
    }

    #[test]
    fn an_alignment_is_named_by_its_gf_id_and_a_profile_by_its_name_line() {
        assert_eq!(
            names(&[fixture("families.sto")], false).unwrap(),
            names(&[fixture("models.hmm")], false).unwrap(),
            "hmmbuild carries the seed's ID into the profile it builds"
        );
    }

    #[test]
    fn an_unnamed_record_is_skipped_unless_it_is_asked_for() {
        // a blank line by default would put every later name
        // out of step with the record it belongs to, so the
        // caller has to ask for the placeholder
        let path = std::env::temp_dir().join(format!("sail-nm-{}.sto", std::process::id()));
        std::fs::write(
            &path,
            b"# STOCKHOLM 1.0\nseq1 AC\n//\n# STOCKHOLM 1.0\n#=GF ID FAM2\nseq1 AC\n//\n",
        )
        .unwrap();

        assert_eq!(names(std::slice::from_ref(&path), false).unwrap(), ["FAM2"]);
        assert_eq!(
            names(std::slice::from_ref(&path), true).unwrap(),
            ["", "FAM2"]
        );

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn several_files_report_their_names_in_the_order_they_were_given() {
        let found = names(&[fixture("models.hmm"), fixture("models.hmm")], false).unwrap();

        assert_eq!(found, ["SH3_1", "PDZ", "SH3_1", "PDZ"]);
    }
}

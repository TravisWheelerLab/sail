use std::io::Write;
use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::Iterable;

use crate::cli::FormatArg;
use crate::input::{Inputs, dispatch};
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
}

impl NamesArgs {
    pub fn run(self) -> Result<()> {
        let mut out = writer(self.output.as_deref())?;

        for line in self.names()? {
            out.write_all(&line)?;
            out.write_all(b"\n")?;
        }
        out.flush()?;

        Ok(())
    }

    /// One name per record, in file order.
    //
    // bytes rather than String: a FASTA name is not required
    // to be utf-8, and decoding here would make an ordinary
    // latin-1 header unprintable
    fn names(&self) -> Result<Vec<Vec<u8>>> {
        let inputs = Inputs::read(&self.input, self.format)?;
        let mut names = Vec::new();

        for entry in inputs.entries() {
            dispatch!(inputs.format(), entry, |collection, _size, name| {
                for record in collection.iter() {
                    match name(record) {
                        Some(found) => names.push(found.to_vec()),

                        // an alignment need not carry a #=GF
                        // ID, and an empty line for it would
                        // silently shift every name after it
                        // out of step with its record
                        None if self.keep_unnamed => names.push(Vec::new()),
                        None => {}
                    }
                }
            });
        }

        Ok(names)
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
        Ok(NamesArgs {
            input: paths.to_vec(),
            format: None,
            output: None,
            keep_unnamed,
        }
        .names()?
        .into_iter()
        .map(|n| String::from_utf8(n).unwrap())
        .collect())
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

use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::Indexable;
use libsail::index::Index;

use crate::cli::{FormatArg, ReadArgs};
use crate::input::{Backend, Inputs, Needs, dispatch};

#[derive(Args)]
pub struct CountArgs {
    /// files to count, or - for stdin
    #[arg(default_value = "-")]
    pub input: Vec<PathBuf>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    #[command(flatten)]
    pub read: ReadArgs,
}

impl CountArgs {
    pub fn run(self) -> Result<()> {
        for line in self.counts()? {
            println!("{line}");
        }

        Ok(())
    }

    /// One line per input, and a total when there is more than one.
    fn counts(&self) -> Result<Vec<String>> {
        let inputs = Inputs::plan(&self.input, self.format, self.read, Needs::Pass)?;
        let several = inputs.entries().len() > 1;

        let mut lines = Vec::new();
        let mut total = 0;

        for entry in inputs.entries() {
            let n = match inputs.backend() {
                // nothing is parsed: the frames are counted
                // and the bytes of each are dropped as the
                // next one overwrites them
                Backend::Stream => libsail::index::count(entry.reader()?, inputs.format())?,
                Backend::Indexed => Index::build(entry.reader()?, inputs.format())?.len(),
                Backend::Memory => dispatch!(inputs.format(), entry, |collection| collection.len()),
            };
            total += n;

            if several {
                lines.push(format!("{n}\t{}", entry.name()));
            }
        }

        lines.push(if several {
            format!("{total}\ttotal")
        } else {
            total.to_string()
        });

        Ok(lines)
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

    fn count(paths: &[PathBuf], format: Option<FormatArg>) -> Result<Vec<String>> {
        CountArgs {
            input: paths.to_vec(),
            format,
            read: ReadArgs::default(),
        }
        .counts()
    }

    #[test]
    fn each_format_is_counted_in_its_own_records_rather_than_in_lines() {
        // one FASTA sequence, one alignment and one model are
        // each one record, and none of the three is a line
        for (name, want) in [
            ("proteins.fa", "5"),
            ("families.sto", "2"),
            ("models.hmm", "2"),
        ] {
            assert_eq!(count(&[fixture(name)], None).unwrap(), [want], "{name}");
        }
    }

    #[test]
    fn several_files_are_reported_one_per_line_and_then_totalled() {
        let lines = count(&[fixture("proteins.fa"), fixture("proteins.fa")], None).unwrap();

        assert_eq!(lines.len(), 3);
        assert!(lines[0].starts_with("5\t"), "{lines:?}");
        assert!(lines[1].starts_with("5\t"), "{lines:?}");
        assert_eq!(lines[2], "10\ttotal");
    }

    #[test]
    fn a_single_file_is_reported_as_the_bare_number() {
        // so `sail count f.fa` reads straight into a shell
        // variable, with no name and no total to strip
        assert_eq!(count(&[fixture("proteins.fa")], None).unwrap(), ["5"]);
    }

    #[test]
    fn a_mismatched_format_assertion_is_refused_before_anything_is_counted() {
        assert!(count(&[fixture("proteins.fa")], Some(FormatArg::Hmm)).is_err());
    }
}

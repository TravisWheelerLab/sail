use std::path::PathBuf;

use anyhow::{Result, bail};
use clap::Args;
use libsail::format::Format;
use libsail::index::Index;
use libsail::seq::fasta::IndexedFasta;
use libsail::seq::p7hmm::IndexedHmm;
use libsail::seq::stockholm::IndexedStockholm;

use crate::cli::FormatArg;
use crate::input::format_of_path;

#[derive(Args)]
pub struct IndexArgs {
    /// file to index
    pub input: PathBuf,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write the index [default: <input>.saidx]
    #[arg(short, long)]
    pub output: Option<PathBuf>,
}

impl IndexArgs {
    pub fn run(self) -> Result<()> {
        // the one command that will not read stdin: an
        // index is byte offsets into a file that is still
        // there afterwards, and a pipe is neither
        if self.input.as_os_str() == "-" {
            bail!(
                "index needs a path, not stdin: an index is byte offsets into a file that stays put"
            );
        }

        let found = format_of_path(&self.input, self.format)?;

        let output = self
            .output
            .clone()
            .unwrap_or_else(|| Index::path_for(&self.input));

        // one arm per format because the collection type is
        // what carries the name scan, and each is its own
        // type. what the three arms do is the same
        let count = match found {
            Format::Fasta => {
                let mut collection = IndexedFasta::open(&self.input)?;
                collection.name_index()?;
                collection.index().save_as(&self.input, &output)?;

                collection.index().len()
            }
            Format::Stockholm => {
                let mut collection = IndexedStockholm::open(&self.input)?;
                collection.name_index()?;
                collection.index().save_as(&self.input, &output)?;

                collection.index().len()
            }
            Format::Hmm => {
                let mut collection = IndexedHmm::open(&self.input)?;
                collection.name_index()?;
                collection.index().save_as(&self.input, &output)?;

                collection.index().len()
            }
        };

        println!("{count} {found} records indexed to {}", output.display());

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use libsail::collection::Indexable;

    use super::*;

    /// Copy a fixture somewhere writable and hand back the copy.
    //
    // fixtures/ is checked in, and indexing a file writes an
    // index next to it
    fn scratch(tag: &str, fixture: &str) -> PathBuf {
        let from = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../fixtures")
            .join(fixture);
        let to = std::env::temp_dir().join(format!("sail-{tag}-{}-{fixture}", std::process::id()));

        std::fs::copy(&from, &to).unwrap();

        to
    }

    fn write(bytes: &[u8], tag: &str, extension: &str) -> PathBuf {
        let path =
            std::env::temp_dir().join(format!("sail-{tag}-{}.{extension}", std::process::id()));
        std::fs::write(&path, bytes).unwrap();

        path
    }

    fn index(input: &Path, format: Option<FormatArg>) -> Result<Index> {
        IndexArgs {
            input: input.to_path_buf(),
            format,
            output: None,
        }
        .run()?;

        Ok(Index::load(input).unwrap())
    }

    fn names(index: &Index) -> Vec<Option<String>> {
        (0..index.len())
            .map(|n| {
                index
                    .name(n)
                    .map(|n| String::from_utf8(n.to_vec()).unwrap())
            })
            .collect()
    }

    fn clean(input: &Path) {
        std::fs::remove_file(Index::path_for(input)).ok();
        std::fs::remove_file(input).ok();
    }

    #[test]
    fn indexing_a_fasta_file_names_every_record_from_its_header() {
        let input = scratch("fa-names", "proteins.fa");
        let index = index(&input, None).unwrap();

        assert_eq!(index.format(), Format::Fasta);
        assert_eq!(
            names(&index),
            [
                "DLG4_HUMAN",
                "SEM5_CAEEL",
                "NCF2_HUMAN",
                "NHRF1_HUMAN",
                "PPP5_HUMAN",
            ]
            .map(|n| Some(n.to_string()))
        );

        clean(&input);
    }

    #[test]
    fn indexing_an_alignment_names_every_record_from_its_gf_id() {
        let input = scratch("sto-names", "families.sto");
        let index = index(&input, None).unwrap();

        assert_eq!(index.format(), Format::Stockholm);
        assert_eq!(
            names(&index),
            [Some("SH3_1".to_string()), Some("PDZ".to_string())]
        );

        clean(&input);
    }

    #[test]
    fn an_alignment_with_no_gf_id_is_indexed_with_no_name_at_all() {
        // both records of families.sto carry an ID, so the
        // fixture cannot reach this path and a synthetic
        // alignment has to
        let input = write(
            b"# STOCKHOLM 1.0\nseq1 ACGT\nseq2 ACGA\n//\n",
            "sto-unnamed",
            "sto",
        );
        let index = index(&input, None).unwrap();

        assert_eq!(names(&index), [None]);

        clean(&input);
    }

    #[test]
    fn indexing_a_profile_database_names_every_model() {
        let input = scratch("hmm-names", "models.hmm");
        let index = index(&input, None).unwrap();

        assert_eq!(index.format(), Format::Hmm);
        assert_eq!(
            names(&index),
            [Some("SH3_1".to_string()), Some("PDZ".to_string())]
        );

        clean(&input);
    }

    #[test]
    fn the_written_index_reaches_the_records_a_fresh_scan_would() {
        // the end the whole command exists for: an index read
        // back off disk has to address what scanning the file
        // addresses
        let input = scratch("fa-reach", "proteins.fa");
        let loaded = IndexedFasta::with_index(&input, index(&input, None).unwrap()).unwrap();

        let scanned = IndexedFasta::open(&input).unwrap();
        assert_eq!(loaded.len(), scanned.len());
        for n in 0..scanned.len() {
            assert_eq!(loaded.cloned(n), scanned.cloned(n), "record {n}");
        }

        clean(&input);
    }

    #[test]
    fn an_index_is_fresh_against_the_file_it_was_just_built_from() {
        // loading is what refuses an index of a file that has
        // changed, so a load that succeeds is that check
        // passing
        let input = scratch("fa-fresh", "proteins.fa");
        index(&input, None).unwrap();

        assert!(Index::load(&input).is_ok());

        clean(&input);
    }

    #[test]
    fn a_mismatched_format_assertion_is_refused_before_any_index_is_written() {
        // the assertion has to land before the work, so a
        // rejected run leaves nothing behind to be mistaken
        // for an index of the file
        let input = scratch("fa-assert", "proteins.fa");

        assert!(index(&input, Some(FormatArg::Stockholm)).is_err());
        assert!(!Index::path_for(&input).exists());

        clean(&input);
    }

    #[test]
    fn a_matching_format_assertion_indexes_as_usual() {
        let input = scratch("fa-agree", "proteins.fa");

        assert!(index(&input, Some(FormatArg::Fasta)).is_ok());

        clean(&input);
    }

    #[test]
    fn stdin_is_refused_rather_than_read_as_a_file_named_dash() {
        let args = IndexArgs {
            input: PathBuf::from("-"),
            format: None,
            output: None,
        };

        assert!(args.run().is_err());
    }

    #[test]
    fn a_file_in_none_of_the_three_formats_is_refused() {
        let input = write(b"@read1\nACGT\n+\n!!!!\n", "fastq", "fq");

        assert!(index(&input, None).is_err());

        clean(&input);
    }
}

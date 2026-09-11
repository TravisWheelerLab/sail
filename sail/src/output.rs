//! Where an operation writes its records.

use std::fmt::Display;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

/// A stream for `path`, or stdout for `None` and for `-`.
pub fn writer(path: Option<&Path>) -> Result<Box<dyn Write>> {
    let Some(path) = path.filter(|p| p.as_os_str() != "-") else {
        return Ok(Box::new(BufWriter::new(std::io::stdout())));
    };

    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;

    Ok(Box::new(BufWriter::new(file)))
}

/// Write `records` through their own format's writer.
///
/// Takes an iterator rather than a collection, so a run of records selected by
/// position reaches the same writer a whole collection does.
pub fn emit<W: Write>(records: impl IntoIterator<Item = impl Display>, mut out: W) -> Result<()> {
    // Display rather than write_to, so this is one function
    // over all three formats. what it gives up is FASTA's
    // line width, which B6 will need and no table A
    // operation does
    for record in records {
        write!(out, "{record}")?;
    }

    out.flush()?;

    Ok(())
}

// ---

/// A record's NAME, which each format spells differently.
//
// free functions rather than closures in the dispatch macro:
// each returns a borrow of its argument, and a closure
// returning a reference derived from its parameter is a
// place closure inference gives up
pub mod name {
    use libsail::seq::fasta::FastaRecord;
    use libsail::seq::p7hmm::HmmRecord;
    use libsail::seq::stockholm::StockholmRecord;

    pub fn fasta(record: &FastaRecord) -> Option<&[u8]> {
        Some(&record.name)
    }

    /// `None` where the alignment carries no `#=GF ID`, which the spec does
    /// not require.
    pub fn stockholm(record: &StockholmRecord) -> Option<&[u8]> {
        record.id().map(str::as_bytes)
    }

    pub fn hmm(record: &HmmRecord) -> Option<&[u8]> {
        Some(record.header.name.as_bytes())
    }
}

#[cfg(test)]
mod tests {
    use libsail::collection::{Indexable, Iterable};
    use libsail::seq::fasta::Fasta;

    use super::*;

    #[test]
    fn a_collection_written_out_reads_back_as_the_same_records() {
        // the round trip is the test that matters: it pins
        // the writer against the parser with no hand-written
        // expected output in between
        let text = b">a desc\nACGT\n>b\nGGTT\n";
        let fasta = Fasta::new(&text[..]).unwrap();

        let mut out = Vec::new();
        emit(fasta.iter(), &mut out).unwrap();

        let back = Fasta::new(&out[..]).unwrap();
        assert_eq!(back.len(), 2);
        for n in 0..back.len() {
            assert_eq!(back.cloned(n), fasta.cloned(n), "record {n}");
        }
    }

    #[test]
    fn a_dash_is_stdout_and_not_a_file_of_that_name() {
        assert!(writer(Some(Path::new("-"))).is_ok());
        assert!(!Path::new("-").exists(), "a file named - was created");
    }
}

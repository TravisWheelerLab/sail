//! Where an operation writes its records.

use std::borrow::Borrow;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use libsail::format::Format;

/// A stream for `path`, or stdout for `None` and for `-`.
pub fn writer(path: Option<&Path>) -> Result<Box<dyn Write>> {
    let Some(path) = path.filter(|p| p.as_os_str() != "-") else {
        return Ok(Box::new(BufWriter::new(std::io::stdout())));
    };

    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;

    Ok(Box::new(BufWriter::new(file)))
}

/// How one format puts a record on a stream.
//
// a fn pointer rather than a trait: the three write_to
// methods take different arguments, FASTA's a line width
// and the other two none, so there is no one signature to
// put in a trait method
pub type Writer<R> = fn(&R, &mut dyn Write) -> libsail::Result<()>;

/// Write `records` through their own format's writer.
///
/// Takes an iterator rather than a collection, so a run of records selected by
/// position reaches the same writer a whole collection does.
pub fn emit<R, B, W>(
    records: impl IntoIterator<Item = B>,
    write: Writer<R>,
    mut out: W,
) -> Result<()>
where
    B: Borrow<R>,
    W: Write,
{
    for record in records {
        // write_to rather than Display: a Formatter takes
        // &str, so Display cannot carry a record whose name
        // is not utf-8, and a Display that fails panics out
        // of write! instead of returning the error
        write(record.borrow(), &mut out)?;
    }

    out.flush()?;

    Ok(())
}

// ---

/// Write one framed record the way its format's writer writes the parsed
/// record.
pub fn write_framed(format: Format, record: &[u8], w: &mut dyn Write) -> Result<()> {
    use libsail::parse::Parse;
    use libsail::seq::fasta::DEFAULT_LINE_WIDTH;
    use libsail::seq::p7hmm::HmmParser;
    use libsail::seq::stockholm::StockholmParser;

    match format {
        // rewrap writes the same bytes write_to does, without
        // building the record to write them from
        Format::Fasta => libsail::seq::fasta::rewrap(record, w, DEFAULT_LINE_WIDTH),
        Format::Stockholm => write::stockholm(&StockholmParser::parse(record)?, w),
        Format::Hmm => write::hmm(&HmmParser::parse(record)?, w),
    }?;

    Ok(())
}

/// Write one framed record, as it came or re-wrapped.
pub fn put(format: Format, record: &[u8], rewrap: bool, w: &mut dyn Write) -> Result<()> {
    if rewrap {
        write_framed(format, record, w)
    } else {
        // the input's own bytes, parsed by nothing: an
        // operation that only selected this record has
        // nothing to say about how it was wrapped
        w.write_all(record)?;

        Ok(())
    }
}

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

/// A record put on a stream, which each format spells differently.
pub mod write {
    use std::io::Write;

    use libsail::seq::fasta::{DEFAULT_LINE_WIDTH, FastaRecord};
    use libsail::seq::p7hmm::HmmRecord;
    use libsail::seq::stockholm::StockholmRecord;

    pub fn fasta(record: &FastaRecord, w: &mut dyn Write) -> libsail::Result<()> {
        record.write_to(w, DEFAULT_LINE_WIDTH)
    }

    pub fn stockholm(record: &StockholmRecord, w: &mut dyn Write) -> libsail::Result<()> {
        record.write_to(w)
    }

    pub fn hmm(record: &HmmRecord, w: &mut dyn Write) -> libsail::Result<()> {
        record.write_to(w)
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
        emit(fasta.iter(), super::write::fasta, &mut out).unwrap();

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

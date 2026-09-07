//! FASTA records.
//!
//! A record is a name, whatever followed it on the header line, and the
//! residues with all whitespace removed.

use std::fmt::Display;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use crate::Record;
use crate::collection::Indexable;
use crate::error::{Error, Result, bail_parse};
use crate::index::{Delimit, Delimiter, Offset, build_index};
use crate::parse::Parse;
#[cfg(feature = "mmap")]
use crate::source::MmapSource;
use crate::source::{FileSource, Source};

/// Residues per line written by [`FastaRecord::write_to`] and [`Display`].
pub const DEFAULT_LINE_WIDTH: usize = 60;

#[derive(Debug, Default, Clone, PartialEq, Eq, Hash)]
pub struct FastaRecord {
    /// The identifier, as bytes because FASTA specifies no encoding.
    /// [`name_str`](Self::name_str) converts when a caller wants text.
    pub name: Vec<u8>,
    /// Whatever followed the name on the header line, trimmed.
    pub extra: Vec<u8>,
    pub seq: Vec<u8>,
}

impl FastaRecord {
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Reverse the residues in place. Composition is preserved and homology is
    /// destroyed, which is what makes a reversed sequence usable as a decoy.
    pub fn reverse(&mut self) {
        self.seq.reverse();
    }

    /// The name as text, if it is valid UTF-8.
    pub fn name_str(&self) -> Result<&str> {
        as_str(&self.name, "a record name")
    }

    /// The description as text, if it is valid UTF-8.
    pub fn extra_str(&self) -> Result<&str> {
        as_str(&self.extra, "a record description")
    }

    pub fn write_to<W: Write>(&self, mut w: W, width: usize) -> Result<()> {
        w.write_all(b">")?;
        w.write_all(&self.name)?;
        if !self.extra.is_empty() {
            w.write_all(b" ")?;
            w.write_all(&self.extra)?;
        }
        writeln!(w)?;

        let width = width.max(1);
        for chunk in self.seq.chunks(width) {
            w.write_all(chunk)?;
            writeln!(w)?;
        }

        Ok(())
    }
}

impl Display for FastaRecord {
    /// Convenience wrapper over [`FastaRecord::write_to`] at
    /// [`DEFAULT_LINE_WIDTH`]. Use `write_to` where IO errors matter.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut buf = Vec::new();
        self.write_to(&mut buf, DEFAULT_LINE_WIDTH)
            .map_err(|_| std::fmt::Error)?;
        f.write_str(std::str::from_utf8(&buf).map_err(|_| std::fmt::Error)?)
    }
}

/// FASTA's record delimiter: a `>` at the start of a line, opening the record
/// it introduces.
pub struct FastaDelimiter {}

impl Delimit for FastaDelimiter {
    const DELIM: &'static [u8] = b">";
    const DELIM_TYPE: Delimiter = Delimiter::Initiating;
}

impl Record for FastaRecord {}

/// Reads a FASTA record out of the bytes of one record.
pub struct FastaParser {}

impl Parse for FastaParser {
    type Record = FastaRecord;

    /// The exact inverse of [`FastaRecord::write_to`]: a header line, then
    /// residues with all whitespace dropped however they were wrapped.
    fn parse(bytes: &[u8]) -> Result<FastaRecord> {
        let split = bytes
            .iter()
            .position(|&b| b == b'\n')
            .unwrap_or(bytes.len());
        let (header, residues) = bytes.split_at(split);

        // the only way this can fail, and build_index never
        // frames a range that does not begin with the
        // delimiter -- so an indexed range always parses,
        // which is what lets IndexedFasta keep Indexable's
        // promise
        let Some(header) = header.strip_prefix(b">") else {
            bail_parse!(
                "a record starts with '>', found {:?}",
                String::from_utf8_lossy(&bytes[..bytes.len().min(16)])
            )
        };

        // trimming ascii whitespace also drops a trailing
        // \r, so a CRLF file parses. the only concession to
        // CRLF in the crate -- the writer emits \n only
        let header = header.trim_ascii();
        let split = header
            .iter()
            .position(|b| b.is_ascii_whitespace())
            .unwrap_or(header.len());

        Ok(FastaRecord {
            name: header[..split].to_vec(),
            extra: header[split..].trim_ascii().to_vec(),
            seq: residues
                .iter()
                .copied()
                .filter(|b| !b.is_ascii_whitespace())
                .collect(),
        })
    }
}

/// Decode `bytes`, naming what they were so a failure says where it happened.
fn as_str<'a>(bytes: &'a [u8], what: &str) -> Result<&'a str> {
    match std::str::from_utf8(bytes) {
        Ok(s) => Ok(s),
        // an error rather than replacement characters: a
        // mangled name leaves the record unfindable by the
        // name it has
        Err(e) => bail_parse!("{what} is not valid utf-8: {e}"),
    }
}

// ---

/// Every record of a FASTA file, held in memory. Peak memory is roughly the
/// bytes plus the records; use [`IndexedFasta`] for a file too large for that.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Fasta {
    records: Vec<FastaRecord>,
}

impl Fasta {
    /// Read and parse every record from `bytes`.
    pub fn new<R: Read>(mut bytes: R) -> Result<Self> {
        let mut all = Vec::new();
        bytes.read_to_end(&mut all)?;

        let records = build_index::<_, FastaDelimiter>(&all[..])?
            .iter()
            .enumerate()
            .map(|(i, o)| {
                let range = &all[o.start as usize..][..o.n_bytes];

                // defensive: build_index only frames ranges
                // starting with the delimiter, so this cannot
                // currently fail. if it ever can, which record
                // failed matters as much as why
                FastaParser::parse(range)
                    .map_err(|e| Error::parse(format!("record {i} at byte {}: {e}", o.start)))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Fasta { records })
    }

    /// Read and parse every record in the file at `path`.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Fasta::new(File::open(path)?)
    }
}

impl Indexable for Fasta {
    type Record = FastaRecord;

    // the records are right here, so lend one
    type Yield<'a> = &'a FastaRecord;

    fn len(&self) -> usize {
        self.records.len()
    }

    fn get(&self, n: usize) -> Option<&FastaRecord> {
        self.records.get(n)
    }
}

// ---

/// Every record of a FASTA file, addressed on disk instead of held in memory.
/// Keeps 16 bytes per record and reads one on demand, so it suits sparse
/// access. Unix and Windows only; [`MmapFasta`] maps the file instead.
#[cfg(any(unix, windows))]
pub struct IndexedFasta<S: Source = FileSource> {
    source: S,
    offsets: Vec<Offset>,
}

#[cfg(any(unix, windows))]
impl IndexedFasta<FileSource> {
    /// Index the file at `path`, parsing none of it.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let source = FileSource::open(path)?;
        let offsets = build_index::<_, FastaDelimiter>(source.file())?;

        Ok(IndexedFasta { source, offsets })
    }
}

/// The same collection with the file mapped rather than read. See
/// [`MmapSource`] for what that costs.
#[cfg(feature = "mmap")]
pub type MmapFasta = IndexedFasta<MmapSource>;

#[cfg(feature = "mmap")]
impl IndexedFasta<MmapSource> {
    /// Map the file at `path` and index it, parsing none of it.
    ///
    /// # Safety
    ///
    /// The file must not be modified or truncated while the collection is
    /// alive; see [`MmapSource::open`].
    pub unsafe fn map(path: impl AsRef<Path>) -> Result<Self> {
        // SAFETY: forwarded to this function's own contract
        let source = unsafe { MmapSource::open(path)? };
        let offsets = build_index::<_, FastaDelimiter>(source.bytes())?;

        Ok(IndexedFasta { source, offsets })
    }
}

#[cfg(any(unix, windows))]
impl<S: Source> IndexedFasta<S> {
    /// Read every record, returning the first that fails and its index.
    pub fn validate(&self) -> Result<()> {
        // a full pass over the file, which is why open does
        // not do it -- that would be Fasta with extra steps
        for n in 0..self.len() {
            self.try_get(n)
                .map_err(|e| Error::parse(format!("record {n}: {e}")))?;
        }

        Ok(())
    }

    /// The `n`th record, or `None` past the end, keeping the failure reason
    /// that [`Indexable::get`] has to throw away.
    pub fn try_get(&self, n: usize) -> Result<Option<FastaRecord>> {
        let Some(o) = self.offsets.get(n) else {
            return Ok(None);
        };

        Ok(Some(FastaParser::parse(
            &self.source.range(o.start, o.n_bytes)?,
        )?))
    }
}

#[cfg(any(unix, windows))]
impl<S: Source> Indexable for IndexedFasta<S> {
    type Record = FastaRecord;

    // nothing to lend: the record is parsed inside the call
    // and did not exist before it
    type Yield<'a>
        = FastaRecord
    where
        Self: 'a;

    fn len(&self) -> usize {
        self.offsets.len()
    }

    /// Reads and parses the record's bytes on every call. A failure comes back
    /// as `None`; use [`try_get`](IndexedFasta::try_get) when the reason matters.
    fn get(&self, n: usize) -> Option<FastaRecord> {
        self.try_get(n).ok().flatten()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixture;

    fn rec(name: &str, extra: &str, seq: &str) -> FastaRecord {
        FastaRecord {
            name: name.as_bytes().to_vec(),
            extra: extra.as_bytes().to_vec(),
            seq: seq.as_bytes().to_vec(),
        }
    }

    #[test]
    fn length_is_residues_not_bytes_of_text() {
        let r = rec("a", "some description", "ACGTACGT");
        assert_eq!(r.len(), 8);
        assert!(!r.is_empty());
        assert!(rec("b", "", "").is_empty());
    }

    #[test]
    fn reversing_preserves_composition_and_destroys_order() {
        // what makes a reversed sequence usable as a decoy
        let mut r = rec("a", "", "ACGTTT");
        let before = r.len();
        r.reverse();
        assert_eq!(r.seq, b"TTTGCA");
        assert_eq!(r.len(), before);
        r.reverse();
        assert_eq!(r.seq, b"ACGTTT", "reversing twice is the identity");
    }

    #[test]
    fn writing_wraps_residues_at_the_requested_width() {
        let mut out = Vec::new();
        rec("a", "", "ACGTACGTAC").write_to(&mut out, 4).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">a\nACGT\nACGT\nAC\n");
    }

    #[test]
    fn a_width_of_zero_does_not_hang_or_panic() {
        let mut out = Vec::new();
        rec("a", "", "ACGT").write_to(&mut out, 0).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">a\nA\nC\nG\nT\n");
    }

    #[test]
    fn an_empty_description_leaves_no_trailing_space() {
        // a trailing space would come back as part of
        // the description
        let mut out = Vec::new();
        rec("a", "", "AC").write_to(&mut out, 60).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">a\nAC\n");

        let mut out = Vec::new();
        rec("a", "d", "AC").write_to(&mut out, 60).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">a d\nAC\n");
    }

    #[test]
    fn display_matches_write_to_at_the_default_width() {
        let r = rec("a", "d", &"A".repeat(DEFAULT_LINE_WIDTH + 5));
        let mut out = Vec::new();
        r.write_to(&mut out, DEFAULT_LINE_WIDTH).unwrap();
        assert_eq!(r.to_string(), String::from_utf8(out).unwrap());
        assert_eq!(
            r.to_string().lines().count(),
            3,
            "header plus two wrapped lines"
        );
    }

    #[test]
    fn an_empty_sequence_still_writes_its_header() {
        let mut out = Vec::new();
        rec("a", "", "").write_to(&mut out, 60).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">a\n");
    }

    #[test]
    fn the_delimiter_opens_a_record_at_a_line_leading_angle_bracket() {
        use crate::index::{Offset, build_index};

        let text = b">a d\nAC\n>b\nGT\n";
        let offsets = build_index::<_, FastaDelimiter>(&text[..]).unwrap();

        assert_eq!(
            offsets,
            vec![
                Offset {
                    start: 0,
                    n_bytes: 8
                },
                Offset {
                    start: 8,
                    n_bytes: 6
                }
            ],
            "each record runs from its own > to the next"
        );
    }

    #[test]
    fn parsing_a_written_record_returns_what_was_written() {
        // pinning parse against write_to checks both
        // directions with no hand-written expected value.
        // widths of 1 and wider-than-the-sequence are where
        // wrapping and parsing come apart
        for r in [
            rec("a", "", "ACGT"),
            rec("a", "some description here", "ACGTACGTAC"),
            rec("DEMO0001", "alpha subunit, synthetic", &"MKVLA".repeat(40)),
            rec("empty", "", ""),
        ] {
            for width in [1, 3, 7, 60, 1000] {
                let mut out = Vec::new();
                r.write_to(&mut out, width).unwrap();
                assert_eq!(
                    FastaParser::parse(&out).unwrap(),
                    r,
                    "{} at width {width}",
                    r.name_str().unwrap()
                );
            }
        }
    }

    #[test]
    fn residues_lose_their_wrapping_however_it_was_done() {
        let one_line = FastaParser::parse(b">a\nACGTACGT\n").unwrap();
        let wrapped = FastaParser::parse(b">a\nAC\nGT\nAC\nGT\n").unwrap();
        assert_eq!(one_line.seq, b"ACGTACGT");
        assert_eq!(one_line, wrapped, "wrapping is not part of the record");
    }

    #[test]
    fn a_header_with_no_description_leaves_extra_empty() {
        let r = FastaParser::parse(b">DEMO0004\nAC\n").unwrap();
        assert_eq!(r.name, b"DEMO0004");
        assert!(r.extra.is_empty());
    }

    #[test]
    fn a_record_need_not_end_in_a_newline() {
        // the last record of a file that does not end in one
        assert_eq!(
            FastaParser::parse(b">a d\nAC").unwrap(),
            rec("a", "d", "AC")
        );
        // a header and nothing else at all
        assert_eq!(FastaParser::parse(b">a").unwrap(), rec("a", "", ""));
    }

    #[test]
    fn a_space_after_the_angle_bracket_is_not_part_of_the_name() {
        // "> name desc" turns up in real files though write_to
        // never emits it; untrimmed, the name comes out empty
        assert_eq!(
            FastaParser::parse(b"> a desc\nAC\n").unwrap(),
            rec("a", "desc", "AC")
        );
        assert_eq!(
            FastaParser::parse(b">   a\nAC\n").unwrap(),
            rec("a", "", "AC")
        );
    }

    #[test]
    fn a_carriage_return_never_reaches_the_record() {
        // the one concession to CRLF: trimming the header and
        // dropping residue whitespace both absorb it
        let r = FastaParser::parse(b">a desc\r\nACGT\r\nAC\r\n").unwrap();
        assert_eq!(r, rec("a", "desc", "ACGTAC"));
    }

    #[test]
    fn bytes_that_do_not_start_a_record_say_so() {
        let err = FastaParser::parse(b"ACGT\nACGT\n").unwrap_err();
        assert!(matches!(err.kind(), crate::ErrorKind::Parse { .. }));
        assert!(
            err.to_string().contains('>'),
            "the message should name what was expected: {err}"
        );
    }

    #[test]
    fn a_name_that_is_not_utf8_survives_verbatim() {
        // FASTA specifies no encoding, so latin-1 is ordinary
        // data rather than an error. keeping the bytes is
        // what makes parsing an indexed range infallible
        let r = FastaParser::parse(b">a\xffb desc\nAC\n").unwrap();
        assert_eq!(r.name, b"a\xffb");
        assert_eq!(r.extra, b"desc");
        assert!(r.name_str().is_err(), "asking for text still reports it");
        assert!(r.name_str().unwrap_err().to_string().contains("utf-8"));

        // and it round-trips, which a lossy conversion would not
        let mut out = Vec::new();
        r.write_to(&mut out, 60).unwrap();
        assert_eq!(FastaParser::parse(&out).unwrap(), r);
    }

    #[test]
    fn every_record_in_the_real_fixture_parses() {
        use crate::index::build_index;

        let path = fixture!("proteins.fa");
        let bytes = std::fs::read(path).expect(path);

        let offsets = build_index::<_, FastaDelimiter>(&bytes[..]).unwrap();
        let records = offsets
            .iter()
            .map(|o| FastaParser::parse(&bytes[o.start as usize..][..o.n_bytes]).unwrap())
            .collect::<Vec<_>>();

        assert_eq!(records.len(), 5);
        assert_eq!(records.iter().map(|r| r.len()).sum::<usize>(), 2335);

        let bare = records.iter().find(|r| r.name == b"PPP5_HUMAN").unwrap();
        assert!(bare.extra.is_empty(), "PPP5_HUMAN has no description");

        for r in &records {
            assert_eq!(
                r.extra.trim_ascii(),
                r.extra,
                "no stray space in a description"
            );
            assert!(!r.seq.iter().any(|b| b.is_ascii_whitespace()));
        }
    }

    // ---

    use std::borrow::Borrow;

    use crate::collection::Iterable;

    /// Distinct lengths, one without a description, one with wrapped residues
    /// -- enough to tell order, framing and whitespace handling apart.
    const THREE: &[u8] = b">a alpha\nACGT\n>b\nAAAA\nCCCC\n>c gamma\nGG\n";

    fn three() -> Fasta {
        Fasta::new(THREE).unwrap()
    }

    fn names<C: Indexable<Record = FastaRecord>>(c: &C) -> Vec<String> {
        c.iter()
            .map(|r| r.borrow().name_str().unwrap().to_string())
            .collect()
    }

    #[test]
    fn a_file_becomes_its_records_in_file_order() {
        let fa = three();
        assert_eq!(fa.len(), 3);
        assert!(!fa.is_empty());

        // get lends here rather than cloning, so the
        // comparison is against a reference
        assert_eq!(fa.get(0).unwrap(), &rec("a", "alpha", "ACGT"));
        assert_eq!(fa.get(1).unwrap(), &rec("b", "", "AAAACCCC"));
        assert_eq!(fa.get(2).unwrap(), &rec("c", "gamma", "GG"));
        assert!(fa.get(3).is_none(), "past the end is None, not a panic");
    }

    #[test]
    fn an_input_with_no_records_is_an_empty_collection() {
        let fa = Fasta::new(&b""[..]).unwrap();
        assert_eq!(fa.len(), 0);
        assert!(fa.is_empty());
        assert!(fa.get(0).is_none());
        assert_eq!(fa.iter().count(), 0);
    }

    #[test]
    fn iterating_yields_every_record_in_file_order() {
        assert_eq!(names(&three()), ["a", "b", "c"]);
    }

    #[test]
    fn filter_narrows_and_keeps_the_order_it_found() {
        let long = three().filter(|r| r.len() >= 4);
        assert_eq!(long.len(), 2);
        assert_eq!(names(&long), ["a", "b"]);
    }

    #[test]
    fn sort_by_key_reorders_without_disturbing_the_records() {
        let by_len = three().sort_by_key(|r| r.len());
        assert_eq!(names(&by_len), ["c", "a", "b"], "2, 4, then 8 residues");
        assert_eq!(by_len.get(0).unwrap(), &rec("c", "gamma", "GG"));
    }

    #[test]
    fn permute_keeps_every_record_exactly_once() {
        let mut shuffled = names(&three().permute());
        shuffled.sort();
        assert_eq!(shuffled, ["a", "b", "c"]);
    }

    #[test]
    fn sample_returns_as_many_records_as_asked_for() {
        assert_eq!(three().sample(2).len(), 2);
        assert_eq!(three().sample(0).len(), 0);
        assert_eq!(three().sample(3).len(), 3, "all of them is still valid");
    }

    #[test]
    fn the_combinators_compose() {
        // Indexed<Fasta> is itself Indexable, so a
        // narrowed view sorts too
        let out = three().filter(|r| r.name != b"b").sort_by_key(|r| r.len());
        assert_eq!(names(&out), ["c", "a"]);
    }

    #[test]
    fn a_non_utf8_header_no_longer_stops_a_whole_file_loading() {
        // keeping the bytes is what makes an ordinary latin-1
        // header readable, and an indexed range always parse
        let fa = Fasta::new(&b">a\nAC\n>\xffb desc\nGG\n"[..]).unwrap();
        assert_eq!(fa.len(), 2);
        assert_eq!(fa.get(1).unwrap().name, b"\xffb");
        assert_eq!(fa.get(1).unwrap().extra, b"desc");
    }

    #[test]
    fn opening_the_real_fixture_gives_the_records_it_documents() {
        let fa = Fasta::open(fixture!("proteins.fa")).unwrap();

        assert_eq!(fa.len(), 5);
        assert_eq!(fa.iter().map(|r| r.len()).sum::<usize>(), 2335);
        assert_eq!(names(&fa)[0], "DLG4_HUMAN");

        let longest = fa.sort_by_key(|r| std::cmp::Reverse(r.len()));
        assert!(longest.get(0).unwrap().len() >= longest.get(4).unwrap().len());
    }

    // ---

    /// Write `bytes` somewhere `IndexedFasta` can open them, since it needs a
    /// real file descriptor to read at an offset.
    fn temp_fasta(tag: &str, bytes: &[u8]) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("libsail-{tag}-{}.fa", std::process::id()));
        std::fs::write(&p, bytes).unwrap();
        p
    }

    #[test]
    fn the_disk_backed_reader_agrees_with_the_in_memory_one() {
        // the claim of Indexable is that one API covers both
        // backends; this is what checks it beyond the types
        let path = temp_fasta("agree", THREE);
        let mem = Fasta::new(THREE).unwrap();
        let disk = IndexedFasta::open(&path).unwrap();

        assert_eq!(disk.len(), mem.len());
        for n in 0..mem.len() {
            assert_eq!(disk.cloned(n), mem.cloned(n), "record {n}");
        }

        assert!(disk.get(mem.len()).is_none(), "past the end is None");
        assert!(disk.try_get(mem.len()).unwrap().is_none());

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn two_fasta_files_address_as_one_collection() {
        // what proteins_extra.fa is for: a second file, so
        // the aggregate has something to aggregate. both
        // backends, since Aggregate forwards Yield and the
        // two differ in it
        use crate::collection::Aggregate;

        let paths = [fixture!("proteins.fa"), fixture!("proteins_extra.fa")];
        let disk = paths
            .iter()
            .map(IndexedFasta::open)
            .collect::<Result<Aggregate<_>>>()
            .unwrap();

        assert_eq!(disk.len(), 7);
        assert_eq!(names(&disk).len(), 7);
        assert_eq!(names(&disk)[0], "DLG4_HUMAN");
        assert_eq!(names(&disk)[5], "SRC_HUMAN", "the second file starts here");

        // the boundary, from both sides
        assert_eq!(disk.locate(4), Some((0, 4)));
        assert_eq!(disk.locate(5), Some((1, 0)));
        assert_eq!(disk.locate(7), None);

        // and the same files held in memory agree record for
        // record, which is the aggregate keeping Indexable's
        // promise rather than only its shape
        let mem: Aggregate<Fasta> = paths
            .iter()
            .map(Fasta::open)
            .collect::<Result<_>>()
            .unwrap();

        assert_eq!(mem.len(), disk.len());
        for n in 0..mem.len() {
            assert_eq!(mem.cloned(n), disk.cloned(n), "record {n}");
        }

        // the combinators reach across the file boundary
        let longest = mem.sort_by_key(|r| std::cmp::Reverse(r.len()));
        assert!(longest.get(0).unwrap().len() >= longest.get(6).unwrap().len());
    }

    #[test]
    fn the_two_backends_agree_on_the_real_fixture() {
        let path = fixture!("proteins.fa");
        let (mem, disk) = (
            Fasta::open(path).expect(path),
            IndexedFasta::open(path).expect(path),
        );

        assert_eq!(disk.len(), mem.len());
        for n in 0..mem.len() {
            assert_eq!(disk.cloned(n), mem.cloned(n), "record {n}");
        }
    }

    #[test]
    fn a_file_not_ending_in_a_newline_keeps_its_last_residue() {
        // every other record ends in \n and the parser drops
        // trailing whitespace, so a read stopping one byte
        // early would be invisible. here the last byte is a
        // residue and has to survive
        const NO_EOL: &[u8] = b">a\nAC\n>b\nGTGT";

        let path = temp_fasta("noeol", NO_EOL);
        let disk = IndexedFasta::open(&path).unwrap();
        let mem = Fasta::new(NO_EOL).unwrap();

        assert_eq!(disk.get(1).unwrap(), rec("b", "", "GTGT"));
        assert_eq!(disk.cloned(1), mem.cloned(1));

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn try_get_keeps_the_reason_that_get_has_to_drop() {
        let path = temp_fasta("truncated", b">a\nAC\n>b\nGTGT\n");
        let disk = IndexedFasta::open(&path).unwrap();
        assert_eq!(disk.len(), 2);

        // the file shrinks under the open index -- with
        // parsing infallible, the only way a record still fails
        //
        // **note: FileSource only, always. under MmapSource a
        //         truncation raises SIGBUS, so this test
        //         would abort the process rather than fail
        std::fs::write(&path, b">a\nAC\n").unwrap();

        assert!(disk.get(1).is_none(), "get has nowhere to put the reason");
        let err = disk.try_get(1).unwrap_err();
        assert!(err.to_string().contains("fill whole buffer"), "{err}");

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn the_combinators_work_on_the_disk_backed_reader_too() {
        // a backend satisfying len/get but not the rest makes
        // the unified API a fiction
        let path = temp_fasta("combi", THREE);
        let disk = IndexedFasta::open(&path).unwrap();

        assert_eq!(names(&disk), ["a", "b", "c"]);
        assert_eq!((&disk).sample(2).len(), 2);

        // through &disk, so one open file feeds every adapter
        // rather than each one needing its own
        let long = (&disk).filter(|r| r.len() >= 4);
        assert_eq!(names(&long), ["a", "b"]);

        let by_len = (&disk).sort_by_key(|r| r.len());
        assert_eq!(names(&by_len), ["c", "a", "b"]);

        // and by name, which is what an index on disk is for
        let by_name = disk.index_by(|r| r.name.clone());
        assert_eq!(disk.get(by_name[&b"b".to_vec()]).unwrap().seq, b"AAAACCCC");

        std::fs::remove_file(path).ok();
    }

    #[cfg(feature = "mmap")]
    #[test]
    fn the_mapped_backend_agrees_with_the_read_one() {
        // a third backend under one API. Source is only
        // worth having if these agree record for record
        let path = temp_fasta("mapped", THREE);
        let mem = Fasta::new(THREE).unwrap();
        let read = IndexedFasta::open(&path).unwrap();
        // SAFETY: this test owns the file and does not touch
        // it while the map is alive
        let mapped = unsafe { MmapFasta::map(&path) }.unwrap();

        assert_eq!(mapped.len(), mem.len());
        for n in 0..mem.len() {
            assert_eq!(mapped.cloned(n), mem.cloned(n), "record {n} against memory");
            assert_eq!(mapped.cloned(n), read.cloned(n), "record {n} against pread");
        }

        assert!(mapped.get(mem.len()).is_none(), "past the end is None");
        assert!(mapped.try_get(mem.len()).unwrap().is_none());
        mapped.validate().unwrap();
        assert_eq!(names(&mapped), ["a", "b", "c"]);

        // SAFETY: as above
        let sorted = unsafe { MmapFasta::map(&path) }.unwrap().sample(2);
        assert_eq!(sorted.len(), 2, "the combinators reach it too");

        std::fs::remove_file(path).ok();
    }

    #[cfg(feature = "mmap")]
    #[test]
    fn mapping_an_empty_file_is_an_empty_collection_not_a_failure() {
        // mmap of a zero-length file is EINVAL on Linux;
        // memmap2 special-cases it, and this is what says so
        let path = temp_fasta("mapped-empty", b"");
        // SAFETY: this test owns the file and leaves it alone
        let mapped = unsafe { MmapFasta::map(&path) }.unwrap();

        assert!(mapped.is_empty());
        assert!(mapped.get(0).is_none());
        assert_eq!(mapped.iter().count(), 0);

        std::fs::remove_file(path).ok();
    }
}

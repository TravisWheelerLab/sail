//! Stockholm alignment records.
//!
//! Rows are interleaved in the file — one sequence name may appear in several
//! blocks and its pieces concatenate — so a record holds each row whole. Column
//! annotations (`#=GC`) and per-residue annotations (`#=GR`) are the same width
//! as the rows and stay aligned to them.

use std::borrow::Cow;
use std::fmt::Display;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use indexmap::IndexMap;

use crate::Record;
use crate::collection::Indexable;
use crate::error::{Error, Result, bail_parse};
use crate::format::Format;
use crate::index::{Delimit, Delimiter, Index, build_index};
use crate::parse::Parse;
#[cfg(feature = "mmap")]
use crate::source::MmapSource;
use crate::source::{FileSource, Source};

pub const HEADER: &str = "# STOCKHOLM 1.0";

/// One `//`-terminated alignment.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct StockholmRecord {
    /// `#=GF` features as (feature, value), in file order. A list rather than
    /// a map because a feature may repeat.
    pub gf: Vec<(String, String)>,
    /// `#=GS` annotations as (sequence name, feature, value), in file order.
    pub gs: Vec<(String, String, String)>,
    /// `#=GC` column annotations, concatenated across interleaved blocks.
    pub gc: IndexMap<String, String>,
    /// `#=GR` per-residue annotations keyed by (sequence name, feature),
    /// concatenated across interleaved blocks.
    pub gr: IndexMap<(String, String), String>,
    /// Row names, in the order they first appear.
    pub names: Vec<String>,
    /// Aligned rows, all of equal width.
    pub seqs: Vec<Vec<u8>>,
    index: IndexMap<String, usize>,
}

impl StockholmRecord {
    /// Alignment width in columns.
    pub fn width(&self) -> usize {
        self.seqs.first().map(|s| s.len()).unwrap_or(0)
    }

    /// Number of aligned rows.
    pub fn depth(&self) -> usize {
        self.seqs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seqs.is_empty()
    }

    /// The row named `name`.
    pub fn get(&self, name: &str) -> Option<&[u8]> {
        self.index.get(name).map(|&i| self.seqs[i].as_slice())
    }

    /// Residues down one alignment column.
    pub fn column(&self, i: usize) -> impl Iterator<Item = u8> + '_ {
        self.seqs.iter().filter_map(move |s| s.get(i).copied())
    }

    /// The `#=GF ID` value, if the alignment carries one — the spec does not
    /// require it, and `hmmalign` output omits it.
    pub fn id(&self) -> Option<&str> {
        // read out of gf rather than stored beside it:
        // write_to emits gf, so a stored id could be written
        // into a record and then silently dropped
        self.gf_first("ID")
    }

    /// The first value for a `#=GF` feature.
    pub fn gf_first(&self, feature: &str) -> Option<&str> {
        self.gf
            .iter()
            .find(|(f, _)| f == feature)
            .map(|(_, v)| v.as_str())
    }

    /// Every value for a `#=GF` feature, in file order.
    pub fn gf_all<'a>(&'a self, feature: &'a str) -> impl Iterator<Item = &'a str> + 'a {
        self.gf
            .iter()
            .filter(move |(f, _)| f == feature)
            .map(|(_, v)| v.as_str())
    }

    /// The consensus secondary structure, if the alignment carries one.
    pub fn ss_cons(&self) -> Option<&str> {
        self.gc.get("SS_cons").map(|s| s.as_str())
    }

    /// The index of `name`'s row, appending an empty one if it is new. The only
    /// way to add a row: appending to `names`/`seqs` desynchronizes the index.
    pub fn row(&mut self, name: &str) -> usize {
        if let Some(&i) = self.index.get(name) {
            return i;
        }
        let i = self.seqs.len();
        self.index.insert(name.to_string(), i);
        self.names.push(name.to_string());
        self.seqs.push(Vec::new());
        i
    }

    /// Check that every row and annotation spans the same number of columns.
    /// Not enforced by construction: a record is assembled block by block and
    /// is legitimately ragged until its last block lands.
    pub fn validate(&self) -> Result<()> {
        let width = self.width();

        for (name, seq) in self.names.iter().zip(&self.seqs) {
            if seq.len() != width {
                bail_parse!(
                    "sequence {name:?} spans {} columns but the alignment is {width} wide",
                    seq.len()
                );
            }
        }

        for (feature, value) in &self.gc {
            if value.len() != width {
                bail_parse!(
                    "#=GC {feature} spans {} columns but the alignment is {width} wide",
                    value.len()
                );
            }
        }

        for ((name, feature), value) in &self.gr {
            if value.len() != width {
                bail_parse!(
                    "#=GR {name} {feature} spans {} columns but the alignment is {width} wide",
                    value.len()
                );
            }
        }

        Ok(())
    }

    pub fn write_to<W: Write>(&self, mut w: W) -> Result<()> {
        writeln!(w, "{HEADER}")?;

        for (feature, value) in &self.gf {
            writeln!(w, "#=GF {feature} {value}")?;
        }
        for (name, feature, value) in &self.gs {
            writeln!(w, "#=GS {name} {feature} {value}")?;
        }

        for (name, seq) in self.names.iter().zip(&self.seqs) {
            write!(w, "{name} ")?;
            w.write_all(seq)?;
            writeln!(w)?;

            for ((gr_name, feature), value) in &self.gr {
                if gr_name == name {
                    writeln!(w, "#=GR {gr_name} {feature} {value}")?;
                }
            }
        }

        for (feature, value) in &self.gc {
            writeln!(w, "#=GC {feature} {value}")?;
        }

        writeln!(w, "//")?;
        Ok(())
    }
}

impl Display for StockholmRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut buf = Vec::new();
        self.write_to(&mut buf).map_err(|_| std::fmt::Error)?;
        f.write_str(std::str::from_utf8(&buf).map_err(|_| std::fmt::Error)?)
    }
}

// ---

/// Stockholm's record delimiter: `//` at the start of a line, closing the
/// alignment it ends.
pub struct StockholmDelimiter {}

impl Delimit for StockholmDelimiter {
    const DELIM: &'static [u8] = b"//";
    const DELIM_TYPE: Delimiter = Delimiter::Terminating;
}

impl Record for StockholmRecord {}

/// Reads a Stockholm alignment out of the bytes of one record.
pub struct StockholmParser {}

impl Parse for StockholmParser {
    type Record = StockholmRecord;

    /// The inverse of [`StockholmRecord::write_to`], accepting what a real file
    /// adds: any number of blocks, `#` comments, and column padding.
    fn parse(bytes: &[u8]) -> Result<StockholmRecord> {
        // **note: unlike FASTA, this and the bails below are
        //         reachable from a range build_index framed,
        //         because an alignment is a grammar and not
        //         a header plus residues. they are the None
        //         IndexedStockholm::get can return and
        //         IndexedFasta::get cannot
        let Ok(text) = std::str::from_utf8(bytes) else {
            bail_parse!("an alignment is utf-8 text, unlike FASTA's byte names")
        };

        let mut rec = StockholmRecord::default();
        let mut seen_header = false;

        for line in text.lines() {
            // trimming the end drops a trailing \r, so a CRLF
            // file parses. the writer emits \n only
            let line = line.trim_end();

            let Some((tag, rest)) = split_first(line) else {
                continue;
            };

            // a terminator need not end its line, so this
            // tests the prefix; see index::Delimit::DELIM
            if tag.starts_with("//") {
                break;
            }

            if !seen_header {
                if !line.trim_start().starts_with("# STOCKHOLM") {
                    bail_parse!(
                        "an alignment starts with {HEADER:?}, found {:?}",
                        line.chars().take(32).collect::<String>()
                    )
                }
                seen_header = true;
                continue;
            }

            match tag {
                "#=GF" => {
                    let (feature, value) = split_two(rest, "#=GF", line)?;
                    rec.gf.push((feature.to_string(), value.to_string()));
                }
                "#=GS" => {
                    let (name, rest) = split_two(rest, "#=GS", line)?;
                    let (feature, value) = split_two(rest, "#=GS", line)?;
                    rec.gs
                        .push((name.to_string(), feature.to_string(), value.to_string()));
                }
                "#=GC" => {
                    let (feature, value) = split_two(rest, "#=GC", line)?;
                    // concatenated, not replaced: one feature
                    // reappears in every block and the pieces
                    // are consecutive stretches of one row
                    rec.gc
                        .entry(feature.to_string())
                        .or_default()
                        .push_str(aligned(value, line)?);
                }
                "#=GR" => {
                    let (name, rest) = split_two(rest, "#=GR", line)?;
                    let (feature, value) = split_two(rest, "#=GR", line)?;
                    rec.gr
                        .entry((name.to_string(), feature.to_string()))
                        .or_default()
                        .push_str(aligned(value, line)?);
                }
                _ if tag.starts_with('#') => {}
                name => {
                    let residues = aligned(rest, line)?;

                    // through row(), which is what keeps the
                    // private index in step with names and seqs
                    let i = rec.row(name);
                    rec.seqs[i].extend_from_slice(residues.as_bytes());
                }
            }
        }

        if !seen_header {
            bail_parse!("an alignment starts with {HEADER:?}, found none")
        }

        Ok(rec)
    }
}

/// The `#=GF ID` value in the bytes of one framed alignment, read without
/// parsing it. `None` where the alignment carries none.
pub fn name_of(record: &[u8]) -> Option<&[u8]> {
    for line in record.split(|&b| b == b'\n') {
        let Some((tag, rest)) = split_first_bytes(line) else {
            continue;
        };

        // a terminator need not end its line, and the parser
        // stops at one; see index::Delimit::DELIM
        if tag.starts_with(b"//") {
            return None;
        }

        if tag == b"#=GF"
            && let Some((feature, value)) = split_first_bytes(rest)
            && feature == b"ID"
        {
            return Some(value);
        }
    }

    None
}

/// [`split_first`] over bytes, so a name can be found without decoding the
/// alignment.
fn split_first_bytes(line: &[u8]) -> Option<(&[u8], &[u8])> {
    let line = line.trim_ascii_start();
    if line.is_empty() {
        return None;
    }

    let end = line
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(line.len());

    Some((&line[..end], line[end..].trim_ascii()))
}

/// The first whitespace-separated token and the rest of the line, trimmed.
fn split_first(line: &str) -> Option<(&str, &str)> {
    let line = line.trim_start();
    if line.is_empty() {
        return None;
    }

    let end = line.find(char::is_whitespace).unwrap_or(line.len());

    Some((&line[..end], line[end..].trim()))
}

/// [`split_first`], naming the line that was short of a field.
fn split_two<'a>(rest: &'a str, what: &str, line: &str) -> Result<(&'a str, &'a str)> {
    match split_first(rest) {
        Some(pair) => Ok(pair),
        None => bail_parse!("a {what} line is missing a field: {line:?}"),
    }
}

/// One column-aligned field: whitespace-free, since a space inside one would
/// shift every column after it.
fn aligned<'a>(value: &'a str, line: &str) -> Result<&'a str> {
    if value.is_empty() || value.contains(char::is_whitespace) {
        bail_parse!("expected one run of aligned columns, found {line:?}")
    }

    Ok(value)
}

// ---

/// Every alignment of a Stockholm file, held in memory. Peak memory is roughly
/// the bytes plus the records; [`IndexedStockholm`] suits a larger file.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Stockholm {
    records: Vec<StockholmRecord>,
}

impl Stockholm {
    /// Read and parse every alignment from `bytes`.
    pub fn new<R: Read>(mut bytes: R) -> Result<Self> {
        let mut all = Vec::new();
        bytes.read_to_end(&mut all)?;

        let records = build_index::<_, StockholmDelimiter>(&all[..])?
            .iter()
            .enumerate()
            .map(|(i, o)| {
                let range = &all[o.start as usize..][..o.n_bytes];

                // which alignment failed matters as much as
                // why: the grammar makes this reachable, so
                // the message has to place the record
                StockholmParser::parse(range)
                    .map_err(|e| Error::parse(format!("record {i} at byte {}: {e}", o.start)))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Stockholm { records })
    }

    /// Read and parse every alignment in the file at `path`.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Stockholm::new(File::open(path)?)
    }
}

impl Indexable for Stockholm {
    type Record = StockholmRecord;

    // the alignments are right here, so lend one
    type Yield<'a> = &'a StockholmRecord;

    fn len(&self) -> usize {
        self.records.len()
    }

    fn get(&self, n: usize) -> Option<&StockholmRecord> {
        self.records.get(n)
    }
}

// ---

/// Every alignment of a Stockholm file, addressed on disk instead of held in
/// memory. Keeps 16 bytes per record and reads one on demand, so it suits
/// sparse access. Unix and Windows only; [`MmapStockholm`] maps the file.
///
/// # Panics
///
/// Iterating panics on a record that does not parse, since [`Indexable`] reads
/// a `None` from `get` as a broken contract. [`validate`](Self::validate)
/// clears a file of unknown provenance up front.
#[cfg(any(unix, windows))]
pub struct IndexedStockholm<S: Source = FileSource> {
    source: S,
    index: Index,
}

#[cfg(any(unix, windows))]
impl IndexedStockholm<FileSource> {
    /// Index the file at `path`, parsing none of it.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let source = FileSource::open(path)?;
        let offsets = build_index::<_, StockholmDelimiter>(source.file())?;

        Ok(IndexedStockholm {
            source,
            index: Index::new(Format::Stockholm, offsets),
        })
    }

    /// Address the file at `path` with `index` instead of scanning for it.
    ///
    /// # Errors
    ///
    /// An index of another format. Nothing here checks that it describes this
    /// file; [`Index::load_from`] does that when one comes off disk.
    pub fn with_index(path: impl AsRef<Path>, index: Index) -> Result<Self> {
        if index.format() != Format::Stockholm {
            bail_parse!(
                "an index of {} cannot address {}, which is stockholm",
                index.format(),
                path.as_ref().display()
            )
        }

        let source = FileSource::open(path)?;

        Ok(IndexedStockholm { source, index })
    }
}

/// The same collection with the file mapped rather than read. See
/// [`MmapSource`] for what that costs.
#[cfg(feature = "mmap")]
pub type MmapStockholm = IndexedStockholm<MmapSource>;

#[cfg(feature = "mmap")]
impl IndexedStockholm<MmapSource> {
    /// Map the file at `path` and index it, parsing none of it.
    ///
    /// # Safety
    ///
    /// The file must not be modified or truncated while the collection is
    /// alive; see [`MmapSource::open`].
    pub unsafe fn map(path: impl AsRef<Path>) -> Result<Self> {
        // SAFETY: forwarded to this function's own contract
        let source = unsafe { MmapSource::open(path)? };
        let offsets = build_index::<_, StockholmDelimiter>(source.bytes())?;

        Ok(IndexedStockholm {
            source,
            index: Index::new(Format::Stockholm, offsets),
        })
    }
}

#[cfg(any(unix, windows))]
impl<S: Source> IndexedStockholm<S> {
    /// Where each alignment sits, and what each is called once
    /// [`name_index`](Self::name_index) has run.
    pub fn index(&self) -> &Index {
        &self.index
    }

    /// Name every alignment, by a byte scan of the range the index frames rather
    /// than a parse of it.
    pub fn name_index(&mut self) -> Result<()> {
        // the fields are borrowed apart because name_with
        // takes &mut self.index while the scan reads
        // self.source
        let source = &self.source;

        self.index.name_with(|o| {
            let bytes = source.range(o.start, o.n_bytes)?;

            Ok(name_of(&bytes).map(|name| Cow::Owned(name.to_vec())))
        })
    }

    /// Read every alignment, returning the first that fails and its index.
    pub fn validate(&self) -> Result<()> {
        // a full pass over the file, which is why open does
        // not do it -- that would be Stockholm with extra
        // steps
        for n in 0..self.len() {
            self.try_get(n)
                .map_err(|e| Error::parse(format!("record {n}: {e}")))?;
        }

        Ok(())
    }

    /// The `n`th alignment, or `None` past the end, keeping the failure reason
    /// that [`Indexable::get`] has to throw away.
    pub fn try_get(&self, n: usize) -> Result<Option<StockholmRecord>> {
        let Some(o) = self.index.offset(n) else {
            return Ok(None);
        };

        Ok(Some(StockholmParser::parse(
            &self.source.range(o.start, o.n_bytes)?,
        )?))
    }
}

#[cfg(any(unix, windows))]
impl<S: Source> Indexable for IndexedStockholm<S> {
    type Record = StockholmRecord;

    // nothing to lend: the alignment is parsed inside the call
    // and did not exist before it
    type Yield<'a>
        = StockholmRecord
    where
        Self: 'a;

    fn len(&self) -> usize {
        self.index.len()
    }

    /// Reads and parses the record's bytes on every call. A failure comes back
    /// as `None`; [`try_get`](Self::try_get) keeps the reason.
    fn get(&self, n: usize) -> Option<StockholmRecord> {
        self.try_get(n).ok().flatten()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_name_scan_agrees_with_the_parsed_record_on_every_fixture_record() {
        // two independent readings of one file: the scan
        // reads bytes and the parser reads an alignment, and
        // over real tool output agreeing is what catches
        // either of them drifting
        let path = fixture!("families.sto");
        let bytes = std::fs::read(path).expect(path);
        let offsets = build_index::<_, StockholmDelimiter>(&bytes[..]).unwrap();
        assert_eq!(offsets.len(), 2);

        for o in &offsets {
            let range = &bytes[o.start as usize..][..o.n_bytes];
            let record = StockholmParser::parse(range).unwrap();

            assert_eq!(
                name_of(range),
                record.id().map(str::as_bytes),
                "the alignment at byte {}",
                o.start
            );
        }
    }

    #[test]
    fn the_name_scan_reports_no_name_for_an_alignment_with_no_gf_id() {
        // both fixture records carry an ID, so the file
        // cannot reach this and a synthetic alignment has to
        let bytes = b"# STOCKHOLM 1.0\nseq1 ACGT\nseq2 ACGA\n//\n";

        assert_eq!(name_of(bytes), None);
        assert_eq!(StockholmParser::parse(bytes).unwrap().id(), None);
    }

    #[test]
    fn a_gf_id_after_the_terminator_belongs_to_no_alignment() {
        // the framed range ends at the terminator's line, but
        // a caller handing over more must not be answered
        // with the next record's name
        let bytes = b"# STOCKHOLM 1.0\nseq1 AC\n//\n# STOCKHOLM 1.0\n#=GF ID FAM2\n//\n";

        assert_eq!(name_of(bytes), None);
    }

    use crate::fixture;

    /// Two rows, one column annotation, one per-residue annotation.
    fn aln() -> StockholmRecord {
        let mut rec = StockholmRecord::default();
        rec.gf.push(("ID".into(), "FAM1".into()));
        let a = rec.row("seq1");
        rec.seqs[a].extend_from_slice(b"ACGT");
        let b = rec.row("seq2");
        rec.seqs[b].extend_from_slice(b"AC-T");
        rec.gc.insert("SS_cons".into(), "<..>".into());
        rec.gr.insert(("seq1".into(), "PP".into()), "9876".into());
        rec.gf.push(("DE".into(), "a family".into()));
        rec
    }

    #[test]
    fn row_appends_once_and_finds_the_same_row_again() {
        // the side index has to stay in step with names and seqs
        let mut rec = StockholmRecord::default();
        let first = rec.row("seq1");
        let again = rec.row("seq1");
        assert_eq!(first, again, "asking twice must not append twice");
        assert_eq!(rec.names.len(), 1);
        assert_eq!(rec.seqs.len(), 1);

        let other = rec.row("seq2");
        assert_ne!(first, other);
        assert_eq!(rec.names, ["seq1", "seq2"]);
    }

    #[test]
    fn get_finds_a_row_added_through_row() {
        let rec = aln();
        assert_eq!(rec.get("seq1"), Some(&b"ACGT"[..]));
        assert_eq!(rec.get("seq2"), Some(&b"AC-T"[..]));
        assert_eq!(rec.get("absent"), None);
    }

    #[test]
    fn width_and_depth_are_columns_and_rows() {
        let rec = aln();
        assert_eq!(rec.width(), 4);
        assert_eq!(rec.depth(), 2);
        assert!(!rec.is_empty());
        assert!(StockholmRecord::default().is_empty());
    }

    #[test]
    fn a_column_reads_down_the_alignment() {
        let rec = aln();
        assert_eq!(rec.column(0).collect::<Vec<_>>(), b"AA");
        assert_eq!(rec.column(2).collect::<Vec<_>>(), b"G-");
        assert!(rec.column(99).next().is_none(), "past the width is empty");
    }

    #[test]
    fn gf_features_may_repeat_and_keep_their_order() {
        // Pfam uses several RN, RM and DR lines per
        // family; a map would collapse them to one
        let mut rec = StockholmRecord::default();
        rec.gf.push(("DR".into(), "first".into()));
        rec.gf.push(("DR".into(), "second".into()));
        rec.gf.push(("AC".into(), "PF00001".into()));

        assert_eq!(rec.gf_first("DR"), Some("first"));
        assert_eq!(
            rec.gf_all("DR").collect::<Vec<_>>(),
            ["first", "second"],
            "both survive, in file order"
        );
        assert_eq!(rec.gf_first("nope"), None);
    }

    #[test]
    fn the_consensus_structure_is_a_column_annotation() {
        // Rfam's SS_cons lives in #=GC
        assert_eq!(aln().ss_cons(), Some("<..>"));
        assert_eq!(StockholmRecord::default().ss_cons(), None);
    }

    #[test]
    fn validate_accepts_an_alignment_whose_rows_all_line_up() {
        aln().validate().unwrap();
        StockholmRecord::default().validate().unwrap();
    }

    #[test]
    fn validate_rejects_a_ragged_row() {
        let mut rec = aln();
        let i = rec.row("seq2");
        rec.seqs[i].push(b'A');
        let err = rec.validate().unwrap_err();
        assert!(err.to_string().contains("seq2"), "{err}");
        assert!(err.to_string().contains("5 columns"), "{err}");
    }

    #[test]
    fn validate_rejects_an_annotation_that_does_not_span_the_alignment() {
        let mut rec = aln();
        rec.gc.insert("SS_cons".into(), "<>".into());
        let err = rec.validate().unwrap_err();
        assert!(err.to_string().contains("SS_cons"), "{err}");

        let mut rec = aln();
        rec.gr.insert(("seq1".into(), "PP".into()), "9".into());
        let err = rec.validate().unwrap_err();
        assert!(err.to_string().contains("PP"), "{err}");
    }

    #[test]
    fn writing_emits_a_header_the_rows_and_a_terminator() {
        let mut out = Vec::new();
        aln().write_to(&mut out).unwrap();
        let text = String::from_utf8(out).unwrap();

        assert!(text.starts_with(HEADER), "{text}");
        assert!(text.trim_end().ends_with("//"), "{text}");
        assert_eq!(aln().id(), Some("FAM1"));
        for needle in [
            "#=GF ID FAM1",
            "#=GF DE a family",
            "seq1",
            "AC-T",
            "SS_cons",
            "PP",
        ] {
            assert!(text.contains(needle), "missing {needle:?} in:\n{text}");
        }
    }

    #[test]
    fn an_alignment_with_no_id_still_writes() {
        // the spec does not require #=GF ID, and
        // hmmalign output omits it
        let mut rec = aln();
        rec.gf.retain(|(f, _)| f != "ID");
        assert_eq!(rec.id(), None);
        let mut out = Vec::new();
        rec.write_to(&mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(!text.contains("#=GF ID"), "{text}");
        assert!(text.contains("seq1"), "{text}");
    }

    // ---

    use std::borrow::Borrow;

    use crate::collection::Iterable;

    /// Two alignments: the first interleaved across two blocks with every
    /// annotation kind, the second a single block with none.
    const TWO: &[u8] = b"\
# STOCKHOLM 1.0
#=GF ID   FAM1
#=GF DE   two blocks

seq1     ACGT
seq2     AC-T
#=GR seq1 PP  9876
#=GC SS_cons  <..>

seq1     TTTT
seq2     TT-T
#=GR seq1 PP  5432
#=GC SS_cons  >..>
//
# STOCKHOLM 1.0
#=GF ID   FAM2
seq9     GGGG
//
";

    fn two() -> Stockholm {
        Stockholm::new(TWO).unwrap()
    }

    fn ids<C: Indexable<Record = StockholmRecord>>(c: &C) -> Vec<String> {
        c.iter()
            .map(|r| r.borrow().id().unwrap_or("?").to_string())
            .collect()
    }

    #[test]
    fn parsing_a_written_alignment_returns_what_was_written() {
        // pinning parse against write_to checks both
        // directions with no hand-written expected value
        let mut out = Vec::new();
        aln().write_to(&mut out).unwrap();

        assert_eq!(StockholmParser::parse(&out).unwrap(), aln());
    }

    #[test]
    fn an_alignment_with_no_rows_at_all_still_round_trips() {
        // #=GF only, which is what a stub family looks like
        let mut rec = StockholmRecord::default();
        rec.gf.push(("ID".into(), "STUB".into()));

        let mut out = Vec::new();
        rec.write_to(&mut out).unwrap();

        assert_eq!(StockholmParser::parse(&out).unwrap(), rec);
    }

    #[test]
    fn interleaved_blocks_concatenate_into_one_row_each() {
        // the reason a record holds each row whole: a name
        // reappears in every block and its pieces are
        // consecutive stretches of the same sequence
        // bound, because get lends out of the collection
        // rather than handing back a copy of a record
        let sto = two();
        let rec = sto.get(0).unwrap();

        assert_eq!(rec.depth(), 2, "two rows, not four");
        assert_eq!(rec.get("seq1"), Some(&b"ACGTTTTT"[..]));
        assert_eq!(rec.get("seq2"), Some(&b"AC-TTT-T"[..]));
        assert_eq!(rec.ss_cons(), Some("<..>>..>"));
        assert_eq!(
            rec.gr.get(&("seq1".to_string(), "PP".to_string())),
            Some(&"98765432".to_string())
        );

        rec.validate()
            .expect("the rows and both annotations span 8 columns");
        assert_eq!(rec.width(), 8);
    }

    #[test]
    fn gf_and_gs_values_keep_the_spaces_inside_them() {
        // a #=GF DE is free text; only the feature name is a
        // token, so splitting on every space would truncate it
        let rec = StockholmParser::parse(
            b"# STOCKHOLM 1.0\n\
              #=GF DE   a family of things\n\
              #=GS seq1 DE alpha subunit\n\
              seq1 AC\n//\n",
        )
        .unwrap();

        assert_eq!(rec.gf_first("DE"), Some("a family of things"));
        assert_eq!(
            rec.gs,
            [(
                "seq1".to_string(),
                "DE".to_string(),
                "alpha subunit".to_string()
            )]
        );
    }

    #[test]
    fn a_plain_comment_is_not_a_row() {
        // "# some note" is a legal Stockholm comment, and
        // reading it as a row would invent a sequence named #
        let rec = StockholmParser::parse(b"# STOCKHOLM 1.0\n# a note\nseq1 AC\n//\n").unwrap();

        assert_eq!(rec.names, ["seq1"]);
        assert_eq!(rec.depth(), 1);
    }

    #[test]
    fn a_terminator_ends_the_alignment_and_nothing_after_it_is_read() {
        // parse is handed one framed record, but a caller may
        // hand it a whole file; the second family must not
        // bleed into the first
        let rec = StockholmParser::parse(TWO).unwrap();

        assert_eq!(rec.id(), Some("FAM1"));
        assert!(rec.get("seq9").is_none(), "that row is the next family's");
    }

    #[test]
    fn bytes_that_do_not_start_an_alignment_say_so() {
        for bad in [&b">a\nACGT\n"[..], &b""[..], &b"//\n"[..]] {
            let err = StockholmParser::parse(bad).unwrap_err();
            assert!(matches!(err.kind(), crate::ErrorKind::Parse { .. }));
            assert!(
                err.to_string().contains("STOCKHOLM"),
                "the message should name what was expected: {err}"
            );
        }
    }

    #[test]
    fn a_row_whose_residues_hold_a_space_is_rejected() {
        // aligned text is one run of columns; a space in it
        // shifts every column after it, and silently absorbing
        // the tail would corrupt the width invariant unseen
        let err = StockholmParser::parse(b"# STOCKHOLM 1.0\nseq1 AC GT\n//\n").unwrap_err();
        assert!(err.to_string().contains("aligned columns"), "{err}");

        let err = StockholmParser::parse(b"# STOCKHOLM 1.0\n#=GC SS_cons <. .>\n//\n").unwrap_err();
        assert!(err.to_string().contains("aligned columns"), "{err}");
    }

    #[test]
    fn a_carriage_return_never_reaches_a_row() {
        let rec = StockholmParser::parse(b"# STOCKHOLM 1.0\r\nseq1 ACGT\r\n//\r\n").unwrap();
        assert_eq!(rec.get("seq1"), Some(&b"ACGT"[..]));
    }

    #[test]
    fn the_delimiter_closes_the_alignment_it_ends() {
        use crate::index::build_index;

        let offsets = build_index::<_, StockholmDelimiter>(TWO).unwrap();
        assert_eq!(offsets.len(), 2, "one record per terminator");
        assert_eq!(offsets[0].start, 0);
        assert_eq!(
            &TWO[offsets[0].start as usize..][..offsets[0].n_bytes]
                .strip_suffix(b"//\n")
                .is_some(),
            &true,
            "a record ends with its own terminator line"
        );
    }

    #[test]
    fn a_file_becomes_its_alignments_in_file_order() {
        let sto = two();
        assert_eq!(sto.len(), 2);
        assert!(!sto.is_empty());
        assert_eq!(ids(&sto), ["FAM1", "FAM2"]);
        assert!(sto.get(2).is_none(), "past the end is None, not a panic");
    }

    #[test]
    fn an_input_with_no_alignments_is_an_empty_collection() {
        // no terminator means no framed record, so the bytes
        // belong to nothing rather than to a ragged last one
        let sto = Stockholm::new(&b""[..]).unwrap();
        assert_eq!(sto.len(), 0);
        assert!(sto.is_empty());
        assert_eq!(sto.iter().count(), 0);

        assert_eq!(
            Stockholm::new(&b"# STOCKHOLM 1.0\nseq1 AC\n"[..])
                .unwrap()
                .len(),
            0
        );
    }

    #[test]
    fn a_malformed_alignment_names_the_record_that_failed() {
        // the in-memory reader parses everything up front, so
        // it can refuse the file rather than hand back a
        // collection with a hole in it
        let err =
            Stockholm::new(&b"# STOCKHOLM 1.0\nseq1 AC\n//\nnot stockholm\n//\n"[..]).unwrap_err();

        assert!(err.to_string().contains("record 1"), "{err}");
        assert!(err.to_string().contains("byte"), "{err}");
    }

    #[test]
    fn the_combinators_work_on_alignments_too() {
        let wide = two().filter(|r| r.width() > 4);
        assert_eq!(ids(&wide), ["FAM1"]);

        let by_depth = two().sort_by_key(|r| r.depth());
        assert_eq!(ids(&by_depth), ["FAM2", "FAM1"], "one row, then two");
        assert_eq!(two().sample(2).len(), 2);
        assert_eq!(two().permute().len(), 2);
    }

    #[test]
    fn every_alignment_in_the_real_fixture_parses() {
        let path = fixture!("families.sto");
        let sto = Stockholm::open(path).expect(path);

        // the Pfam SEEDs for SH3_1 and PDZ, six sequences
        // kept from each by fixtures/make.sh
        assert_eq!(sto.len(), 2);
        assert_eq!(ids(&sto), ["SH3_1", "PDZ"]);

        let sh3 = sto.get(0).unwrap();
        assert_eq!(sh3.depth(), 6);
        assert_eq!(sh3.width(), 56);
        assert_eq!(
            sh3.get("VAV_HUMAN/788-834"),
            Some(&b"KARYDFCARD..RSELSLKEGDIIKILNKKGQ..QGWWRGEIY.....GRVGWFPA"[..])
        );
        assert_eq!(
            sh3.gf_all("DR").collect::<Vec<_>>(),
            [
                "INTERPRO; IPR001452;",
                "PROSITE; PDOC50002;",
                "PRINTS; PR00452;",
                "SCOP; 1shf; sf;",
                "HOMSTRAD; sh3;",
                "SO; 0000417; polypeptide_domain;",
            ],
            "a repeated feature keeps every value, in file order"
        );
        assert_eq!(sh3.gs.len(), 6, "one #=GS AC per kept sequence");

        for n in 0..sto.len() {
            sto.get(n).unwrap().validate().unwrap();
        }
    }

    // ---

    /// Write `bytes` somewhere [`IndexedStockholm`] can open them, since it
    /// needs a real file descriptor to read at an offset.
    fn temp_sto(tag: &str, bytes: &[u8]) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("libsail-{tag}-{}.sto", std::process::id()));
        std::fs::write(&p, bytes).unwrap();
        p
    }

    #[test]
    fn the_disk_backed_reader_agrees_with_the_in_memory_one() {
        // the claim of Indexable is that one API covers both
        // backends; this is what checks it beyond the types
        let path = temp_sto("agree", TWO);
        let mem = two();
        let disk = IndexedStockholm::open(&path).unwrap();

        assert_eq!(disk.len(), mem.len());
        for n in 0..mem.len() {
            assert_eq!(disk.cloned(n), mem.cloned(n), "record {n}");
        }

        assert!(disk.get(mem.len()).is_none(), "past the end is None");
        assert!(disk.try_get(mem.len()).unwrap().is_none());
        disk.validate().unwrap();

        assert_eq!(ids(&disk), ["FAM1", "FAM2"]);
        assert_eq!(
            IndexedStockholm::open(&path)
                .unwrap()
                .sort_by_key(|r| r.depth())
                .get(0)
                .unwrap()
                .id(),
            Some("FAM2")
        );

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn validate_finds_the_record_get_can_only_call_none() {
        // the gap FASTA does not have: parsing an indexed
        // range can fail here, so get() has a None that
        // Indexable treats as a broken contract
        let path = temp_sto("bad", b"# STOCKHOLM 1.0\nseq1 AC\n//\nnot stockholm\n//\n");
        let disk = IndexedStockholm::open(&path).unwrap();

        assert_eq!(disk.len(), 2);
        assert!(disk.get(0).is_some());
        assert!(disk.get(1).is_none(), "get has nowhere to put the reason");

        let err = disk.validate().unwrap_err();
        assert!(err.to_string().contains("record 1"), "{err}");
        assert!(err.to_string().contains("STOCKHOLM"), "{err}");

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn a_file_that_shrinks_under_the_index_is_an_io_failure_not_a_parse_one() {
        // the other way a record fails here: the offsets
        // outlive the bytes they were taken from, and the
        // positioned read is what reports it
        //
        // **note: FileSource only, always. under MmapSource a
        //         truncation raises SIGBUS, so this test
        //         would abort the process rather than fail
        let path = temp_sto("truncated", TWO);
        let disk = IndexedStockholm::open(&path).unwrap();
        assert_eq!(disk.len(), 2);

        std::fs::write(&path, b"# STOCKHOLM 1.0\nseq1 AC\n//\n").unwrap();

        assert!(disk.get(1).is_none(), "get has nowhere to put the reason");
        let err = disk.try_get(1).unwrap_err();
        assert!(matches!(err.kind(), crate::ErrorKind::Io { .. }), "{err}");
        assert!(err.to_string().contains("fill whole buffer"), "{err}");

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn the_two_backends_agree_on_the_real_fixture() {
        let path = fixture!("families.sto");
        let (mem, disk) = (
            Stockholm::open(path).expect(path),
            IndexedStockholm::open(path).expect(path),
        );

        assert_eq!(disk.len(), mem.len());
        for n in 0..mem.len() {
            assert_eq!(disk.cloned(n), mem.cloned(n), "record {n}");
        }
    }

    #[test]
    fn an_index_handed_in_addresses_the_file_the_same_as_one_scanned_for() {
        // what writing an index down is worth: with_index
        // skips the scan, and nothing downstream may be able
        // to tell that it did
        let path = fixture!("families.sto");
        let scanned = IndexedStockholm::open(path).expect(path);

        let offsets = build_index::<_, StockholmDelimiter>(File::open(path).unwrap()).unwrap();
        let handed =
            IndexedStockholm::with_index(path, Index::new(Format::Stockholm, offsets)).expect(path);

        assert_eq!(handed.len(), scanned.len());
        for n in 0..scanned.len() {
            assert_eq!(handed.cloned(n), scanned.cloned(n), "record {n}");
        }
    }

    #[test]
    fn an_index_of_another_format_is_refused() {
        // offsets alone cannot say what they frame, so the
        // format the index carries is the only thing between
        // a hmm index and a file it does not describe
        let path = fixture!("families.sto");
        let offsets = build_index::<_, StockholmDelimiter>(File::open(path).unwrap()).unwrap();

        assert!(IndexedStockholm::with_index(path, Index::new(Format::Hmm, offsets)).is_err());
    }

    #[test]
    fn naming_the_index_finds_every_record_by_the_name_its_parser_reports() {
        let path = fixture!("families.sto");
        let mut collection = IndexedStockholm::open(path).expect(path);
        collection.name_index().unwrap();

        for n in 0..collection.len() {
            let Some(name) = collection.index().name(n) else {
                continue;
            };

            assert_eq!(collection.index().position(name), Some(n), "record {n}");
        }
    }

    #[cfg(feature = "mmap")]
    #[test]
    fn the_mapped_backend_agrees_with_the_read_one() {
        // a third backend under one API. Source is only
        // worth having if these agree alignment for alignment
        let path = temp_sto("mapped", TWO);
        let mem = two();
        let read = IndexedStockholm::open(&path).unwrap();
        // SAFETY: this test owns the file and does not touch
        // it while the map is alive
        let mapped = unsafe { MmapStockholm::map(&path) }.unwrap();

        assert_eq!(mapped.len(), mem.len());
        for n in 0..mem.len() {
            assert_eq!(
                mapped.cloned(n),
                mem.cloned(n),
                "alignment {n} against memory"
            );
            assert_eq!(
                mapped.cloned(n),
                read.cloned(n),
                "alignment {n} against pread"
            );
        }

        assert!(mapped.get(mem.len()).is_none(), "past the end is None");
        assert!(mapped.try_get(mem.len()).unwrap().is_none());
        mapped.validate().unwrap();
        assert_eq!(ids(&mapped), ["FAM1", "FAM2"]);

        // SAFETY: as above
        let sorted = unsafe { MmapStockholm::map(&path) }.unwrap().sample(2);
        assert_eq!(sorted.len(), 2, "the combinators reach it too");

        std::fs::remove_file(path).ok();
    }
}

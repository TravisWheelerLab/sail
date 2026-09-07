//! Pairwise search result tables.
//!
//! Every tool here reports the same thing — one query matched one target, with
//! a score, an E-value, and optionally where the match sat in each sequence.
//! What differs is which column holds what, so a [`Hit`] is the shape and a
//! [`HitColumns`] impl is one tool\'s layout of it.

use std::fmt::Display;
use std::fs::File;
use std::io::{Read, Write};
use std::marker::PhantomData;
use std::path::Path;

use crate::Record;
use crate::collection::Indexable;
use crate::error::{Error, Result, bail_parse};
use crate::index::Offset;
use crate::lines::line_index;
use crate::parse::Parse;
#[cfg(feature = "mmap")]
use crate::source::MmapSource;
use crate::source::{FileSource, Source};

pub mod blast;
pub mod hmmer;
pub mod nail;

/// A coordinate span as these formats write them: 1-based and inclusive on
/// both ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Span {
    /// The first position covered.
    pub start: usize,

    /// The last position covered, which may be **below** [`start`](Self::start).
    //
    // why this is not a Range: ranges are half-open and
    // imply start <= end, but BLAST and nhmmer report a
    // minus-strand match with start > end, which a Range
    // renders empty
    pub end: usize,
}

#[allow(
    clippy::len_without_is_empty,
    reason = "a span always covers at least one position, so `is_empty` would \
              be a method that only ever returns false"
)]
impl Span {
    pub fn new(start: usize, end: usize) -> Self {
        Span { start, end }
    }

    /// Number of positions covered, orientation ignored.
    pub fn len(&self) -> usize {
        self.start.abs_diff(self.end) + 1
    }

    /// Whether the span runs backwards, which these formats use to mean the
    /// match is on the minus strand.
    pub fn is_reverse(&self) -> bool {
        self.start > self.end
    }
}

impl Display for Span {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}-{}", self.start, self.end)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Hit {
    pub query: String,
    pub target: String,
    /// Present when the layout reports query coordinates.
    pub query_span: Option<Span>,
    /// Present when the layout reports target coordinates.
    pub target_span: Option<Span>,
    pub score: f32,
    pub e_value: f64,
}

impl Display for Hit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {} {:.2e}", self.query, self.target, self.e_value)
    }
}

/// Where each value sits in one row of a particular table format.
pub trait HitColumns {
    const QUERY: usize;
    const TARGET: usize;
    const Q_START: Option<usize>;
    const Q_END: Option<usize>;
    const T_START: Option<usize>;
    const T_END: Option<usize>;
    const SCORE: usize;
    const E_VALUE: usize;
    /// Mandatory fields in a row.
    //
    // a row with fewer is an error rather than a panic, which
    // catches a table in the wrong format at its first row
    const N_FIELDS: usize;
}

// ---

/// A `#` opens a comment in every layout the crate models.
fn is_comment(b: u8) -> bool {
    b == b'#'
}

impl Record for Hit {}

/// Reads a [`Hit`] out of one row, in `C`'s column order.
//
// one impl covers every layout, including one a caller
// writes. a blanket impl<C: HitColumns> Parse for C is
// not available: coherence does no negative reasoning, so
// it would collide with impl Parse for FastaParser
pub struct HitParser<C>(PhantomData<fn() -> C>);

impl<C: HitColumns> Parse for HitParser<C> {
    type Record = Hit;

    /// The inverse of [`Hit::write_to`] for the same layout.
    fn parse(bytes: &[u8]) -> Result<Hit> {
        let Ok(text) = std::str::from_utf8(bytes) else {
            bail_parse!("a row is utf-8 text")
        };

        // whitespace and not tabs: BLAST writes tabs, nail
        // and HMMER write aligned runs of spaces, and none
        // of them permit a space inside a modelled column
        let fields: Vec<&str> = text.split_whitespace().collect();

        // >= and not ==: HMMER's description is free text past
        // the last column, and `-outfmt "6 std qlen slen"`
        // appends columns after the standard twelve, which
        // this reads correctly and == would reject
        if fields.len() < C::N_FIELDS {
            bail_parse!(
                "a row of this layout has {} fields, found {}: {:?}",
                C::N_FIELDS,
                fields.len(),
                text.trim_end()
            )
        }

        Ok(Hit {
            query: field(&fields, C::QUERY, "query")?.to_string(),
            target: field(&fields, C::TARGET, "target")?.to_string(),
            query_span: span(&fields, C::Q_START, C::Q_END, "query")?,
            target_span: span(&fields, C::T_START, C::T_END, "target")?,
            score: number(field(&fields, C::SCORE, "score")?, "score")?,
            e_value: number(field(&fields, C::E_VALUE, "e-value")?, "e-value")?,
        })
    }
}

/// Column `i`, or a failure naming the column rather than the index.
fn field<'a>(fields: &[&'a str], i: usize, what: &str) -> Result<&'a str> {
    match fields.get(i) {
        Some(f) => Ok(*f),
        // reachable only from a layout whose index is past
        // its own N_FIELDS.
        // every_layout_indexes_within_its_own_row rules that
        // out for the four here, not for a caller's
        None => bail_parse!("this layout puts {what} at column {i}, past the row"),
    }
}

/// Both ends of a span, or `None` where the layout reports no coordinates.
fn span(
    fields: &[&str],
    start: Option<usize>,
    end: Option<usize>,
    what: &str,
) -> Result<Option<Span>> {
    match (start, end) {
        (Some(s), Some(e)) => Ok(Some(Span::new(
            number(field(fields, s, what)?, what)?,
            number(field(fields, e, what)?, what)?,
        ))),
        (None, None) => Ok(None),

        // a layout naming one end and not the other would
        // otherwise drop the coordinate it does report
        _ => bail_parse!("this layout reports half a {what} span"),
    }
}

fn number<T>(value: &str, what: &str) -> Result<T>
where
    T: std::str::FromStr,
    T::Err: Display,
{
    match value.parse() {
        Ok(v) => Ok(v),
        Err(e) => bail_parse!("{what}: {value:?} is not a number: {e}"),
    }
}

impl Hit {
    /// Write this hit as one row in `C`'s column order. Round-trips a [`Hit`],
    /// not a file: columns the layout has that a `Hit` lacks are written `-`.
    pub fn write_to<C: HitColumns, W: Write>(&self, mut w: W) -> Result<()> {
        let mut fields = vec!["-".to_string(); C::N_FIELDS];

        let mut set = |i: usize, value: String| -> Result<()> {
            match fields.get_mut(i) {
                Some(f) => {
                    *f = value;
                    Ok(())
                }
                None => bail_parse!(
                    "this layout puts a column at {i}, past its own {} fields",
                    C::N_FIELDS
                ),
            }
        };

        set(C::QUERY, self.query.clone())?;
        set(C::TARGET, self.target.clone())?;

        // Display for a float is the shortest text that reads
        // back as the same value, which is what the round trip
        // needs. {:e} because an e-value like 3.2e-88 spells
        // out to ninety digits otherwise
        set(C::SCORE, self.score.to_string())?;
        set(C::E_VALUE, format!("{:e}", self.e_value))?;

        for (start, end, span) in [
            (C::Q_START, C::Q_END, self.query_span),
            (C::T_START, C::T_END, self.target_span),
        ] {
            // a span is dropped where the layout has
            // nowhere to put it -- HmmerTable reports no
            // coordinates at all -- so a round trip through
            // one needs a Hit that layout can hold
            if let (Some(i), Some(j), Some(span)) = (start, end, span) {
                set(i, span.start.to_string())?;
                set(j, span.end.to_string())?;
            }
        }

        writeln!(w, "{}", fields.join("\t"))?;

        Ok(())
    }
}

// ---

/// Every row of a hit table, held in memory.
pub struct Table<C: HitColumns> {
    hits: Vec<Hit>,

    // fn() -> C rather than C: the layout is a parsing
    // decision and not something the collection owns, and
    // this keeps the auto traits without asking C for them
    layout: PhantomData<fn() -> C>,
}

/// BLAST `-outfmt 6`, in memory.
pub type BlastHits = Table<crate::tbl::blast::BlastTable>;
/// nail, in memory.
pub type NailHits = Table<crate::tbl::nail::NailTable>;
/// HMMER `--tblout`, in memory.
pub type HmmerHits = Table<crate::tbl::hmmer::HmmerTable>;
/// HMMER `--domtblout`, in memory.
pub type HmmerDomHits = Table<crate::tbl::hmmer::HmmerDomTable>;

impl<C: HitColumns> Table<C> {
    /// Read and parse every row from `bytes`.
    pub fn new<R: Read>(mut bytes: R) -> Result<Self> {
        let mut all = Vec::new();
        bytes.read_to_end(&mut all)?;

        let hits = line_index(&all[..], is_comment)?
            .iter()
            .enumerate()
            .map(|(i, o)| {
                let range = &all[o.start as usize..][..o.n_bytes];

                // which row failed matters as much as why: a
                // table in the wrong layout fails at its first
                // one, and the byte offset is how to find it
                HitParser::<C>::parse(range)
                    .map_err(|e| Error::parse(format!("row {i} at byte {}: {e}", o.start)))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Table {
            hits,
            layout: PhantomData,
        })
    }

    /// Read and parse every row in the file at `path`.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Table::new(File::open(path)?)
    }

    /// Write every row in `C`'s column order. See [`Hit::write_to`] for what a
    /// row can and cannot carry back.
    pub fn write_to<W: Write>(&self, mut w: W) -> Result<()> {
        for hit in &self.hits {
            hit.write_to::<C, _>(&mut w)?;
        }

        Ok(())
    }
}

// written out rather than derived: a derive bounds every
// type parameter, so it would demand a HitColumns impl be
// Debug and Clone to describe a table never holding one
impl<C: HitColumns> std::fmt::Debug for Table<C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Table").field("hits", &self.hits).finish()
    }
}

impl<C: HitColumns> Clone for Table<C> {
    fn clone(&self) -> Self {
        Table {
            hits: self.hits.clone(),
            layout: PhantomData,
        }
    }
}

impl<C: HitColumns> PartialEq for Table<C> {
    fn eq(&self, other: &Self) -> bool {
        self.hits == other.hits
    }
}

impl<C: HitColumns> Default for Table<C> {
    fn default() -> Self {
        Table {
            hits: Vec::new(),
            layout: PhantomData,
        }
    }
}

impl<C: HitColumns> Indexable for Table<C> {
    type Record = Hit;

    // the hits are right here, so lend one
    type Yield<'a>
        = &'a Hit
    where
        Self: 'a;

    fn len(&self) -> usize {
        self.hits.len()
    }

    fn get(&self, n: usize) -> Option<&Hit> {
        self.hits.get(n)
    }
}

// ---

/// Every row of a hit table, addressed on disk instead of held in memory.
/// Keeps 16 bytes per row and reads one on demand. Unix and Windows only.
///
/// # Panics
///
/// Iterating panics on a row that does not parse, since [`Indexable`] reads a
/// `None` from `get` as a broken contract. [`validate`](Self::validate) clears
/// a file of unknown provenance up front.
#[cfg(any(unix, windows))]
pub struct IndexedTable<C: HitColumns, S: Source = FileSource> {
    source: S,
    offsets: Vec<Offset>,
    layout: PhantomData<fn() -> C>,
}

#[cfg(any(unix, windows))]
impl<C: HitColumns> IndexedTable<C, FileSource> {
    /// Index the file at `path`, parsing none of it.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let source = FileSource::open(path)?;
        let offsets = line_index(source.file(), is_comment)?;

        Ok(IndexedTable {
            source,
            offsets,
            layout: PhantomData,
        })
    }
}

#[cfg(feature = "mmap")]
impl<C: HitColumns> IndexedTable<C, MmapSource> {
    /// Map the file at `path` and index it, parsing none of it.
    ///
    /// # Safety
    ///
    /// The file must not be modified or truncated while the collection is
    /// alive; see [`MmapSource::open`].
    pub unsafe fn map(path: impl AsRef<Path>) -> Result<Self> {
        // SAFETY: forwarded to this function's own contract
        let source = unsafe { MmapSource::open(path)? };
        let offsets = line_index(source.bytes(), is_comment)?;

        Ok(IndexedTable {
            source,
            offsets,
            layout: PhantomData,
        })
    }
}

#[cfg(any(unix, windows))]
impl<C: HitColumns, S: Source> IndexedTable<C, S> {
    /// Read every row, returning the first that fails and its index.
    pub fn validate(&self) -> Result<()> {
        for n in 0..self.len() {
            self.try_get(n)
                .map_err(|e| Error::parse(format!("row {n}: {e}")))?;
        }

        Ok(())
    }

    /// The `n`th row, or `None` past the end, keeping the failure reason that
    /// [`Indexable::get`] has to throw away.
    pub fn try_get(&self, n: usize) -> Result<Option<Hit>> {
        let Some(o) = self.offsets.get(n) else {
            return Ok(None);
        };

        Ok(Some(HitParser::<C>::parse(
            &self.source.range(o.start, o.n_bytes)?,
        )?))
    }
}

#[cfg(any(unix, windows))]
impl<C: HitColumns, S: Source> Indexable for IndexedTable<C, S> {
    type Record = Hit;

    // nothing to lend: the hit is parsed inside the call and
    // did not exist before it
    type Yield<'a>
        = Hit
    where
        Self: 'a;

    fn len(&self) -> usize {
        self.offsets.len()
    }

    /// Reads and parses the row's bytes on every call. A failure comes back as
    /// `None`; [`try_get`](Self::try_get) keeps the reason.
    fn get(&self, n: usize) -> Option<Hit> {
        self.try_get(n).ok().flatten()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixture;
    use crate::tbl::blast::BlastTable;
    use crate::tbl::hmmer::{HmmerDomTable, HmmerTable};
    use crate::tbl::nail::NailTable;

    #[test]
    fn a_span_counts_both_endpoints() {
        // 1-based and inclusive at both ends
        assert_eq!(Span::new(1, 1).len(), 1);
        assert_eq!(Span::new(10, 20).len(), 11);
    }

    #[test]
    fn a_backwards_span_is_a_minus_strand_match_not_an_empty_one() {
        // why Span is not a Range: BLAST and nhmmer report a
        // minus-strand match with start > end
        let rev = Span::new(900, 400);
        assert!(rev.is_reverse());
        assert_eq!(rev.len(), 501, "orientation does not change the length");
        assert_eq!(rev.len(), Span::new(400, 900).len());
        assert!(!Span::new(400, 900).is_reverse());
    }

    #[test]
    fn hit_display_is_a_summary_and_says_so() {
        // pinned: Display is a summary line, not a rendering
        // of the record. unlike the seq types it drops the
        // score and both spans, and cannot round-trip
        let hit = Hit {
            query: "q1".into(),
            target: "t1".into(),
            query_span: Some(Span::new(1, 10)),
            target_span: Some(Span::new(5, 14)),
            score: 12.5,
            e_value: 1e-9,
        };
        let shown = hit.to_string();
        assert_eq!(shown, "q1 t1 1.00e-9");
        assert!(!shown.contains("12.5"), "the score is not shown");
    }

    /// Every layout must index within the row width it declares, or a parser
    /// built on it reads the wrong column — or panics.
    fn layout_is_self_consistent<C: HitColumns>(what: &str) {
        let n = C::N_FIELDS;
        assert!(C::QUERY < n, "{what}: QUERY past the row");
        assert!(C::TARGET < n, "{what}: TARGET past the row");
        assert!(C::SCORE < n, "{what}: SCORE past the row");
        assert!(C::E_VALUE < n, "{what}: E_VALUE past the row");
        for (name, col) in [
            ("Q_START", C::Q_START),
            ("Q_END", C::Q_END),
            ("T_START", C::T_START),
            ("T_END", C::T_END),
        ] {
            if let Some(i) = col {
                assert!(i < n, "{what}: {name} past the row");
            }
        }
        // coordinates come in pairs or not at all
        assert_eq!(
            C::Q_START.is_some(),
            C::Q_END.is_some(),
            "{what}: half a query span"
        );
        assert_eq!(
            C::T_START.is_some(),
            C::T_END.is_some(),
            "{what}: half a target span"
        );
    }

    #[test]
    fn every_layout_indexes_within_its_own_row() {
        layout_is_self_consistent::<BlastTable>("blast");
        layout_is_self_consistent::<NailTable>("nail");
        layout_is_self_consistent::<HmmerTable>("hmmer tblout");
        layout_is_self_consistent::<HmmerDomTable>("hmmer domtblout");
    }

    #[test]
    fn the_layouts_are_actually_different() {
        // if two of these agreed on every column,
        // one of them would be wrong
        assert_ne!(BlastTable::N_FIELDS, HmmerTable::N_FIELDS);
        assert_ne!(HmmerTable::N_FIELDS, HmmerDomTable::N_FIELDS);
        assert_ne!(BlastTable::E_VALUE, HmmerTable::E_VALUE);
    }

    // ---

    use std::borrow::Borrow;

    use crate::collection::Iterable;

    /// A hit with both spans, which every layout but [`HmmerTable`] can hold.
    fn hit(query: &str, target: &str) -> Hit {
        Hit {
            query: query.into(),
            target: target.into(),
            query_span: Some(Span::new(12, 101)),
            target_span: Some(Span::new(88, 21)),
            score: 68.4,
            e_value: 1.1e-12,
        }
    }

    fn row_of<C: HitColumns>(h: &Hit) -> Vec<u8> {
        let mut out = Vec::new();
        h.write_to::<C, _>(&mut out).unwrap();
        out
    }

    #[test]
    fn parsing_a_written_row_returns_what_was_written() {
        // pinning parse against write_to checks both
        // directions with no hand-written expected value.
        // the minus-strand span is in here on purpose: a Range
        // would render it empty and it must survive
        for h in [
            hit("DEMO0001", "DEMO0102"),
            Hit {
                score: -0.5,
                e_value: 3.2e-88,
                ..hit("q", "t")
            },
            Hit {
                query_span: Some(Span::new(1, 1)),
                target_span: Some(Span::new(9999, 1)),
                ..hit("long-name_with.punct", "t|2")
            },
        ] {
            assert_eq!(
                HitParser::<BlastTable>::parse(&row_of::<BlastTable>(&h)).unwrap(),
                h
            );
            assert_eq!(
                HitParser::<NailTable>::parse(&row_of::<NailTable>(&h)).unwrap(),
                h
            );
            assert_eq!(
                HitParser::<HmmerDomTable>::parse(&row_of::<HmmerDomTable>(&h)).unwrap(),
                h
            );
        }
    }

    #[test]
    fn a_layout_with_no_coordinates_round_trips_a_hit_without_them() {
        // HmmerTable reports no spans at all, so a hit written
        // through it comes back with None -- which is the
        // round trip for the hits that layout can hold, and a
        // documented loss for the ones it cannot
        let spanned = hit("q", "t");
        let flat = Hit {
            query_span: None,
            target_span: None,
            ..spanned.clone()
        };

        assert_eq!(
            HitParser::<HmmerTable>::parse(&row_of::<HmmerTable>(&flat)).unwrap(),
            flat
        );
        assert_eq!(
            HitParser::<HmmerTable>::parse(&row_of::<HmmerTable>(&spanned)).unwrap(),
            flat,
            "the spans are dropped, not mangled"
        );
    }

    #[test]
    fn a_written_row_says_nothing_about_the_columns_a_hit_does_not_model() {
        // the writer emits a well-formed row of the layout, not
        // a reproduction of a tool's output: BLAST's pident and
        // the rest are '-' because a Hit never held them
        let text = String::from_utf8(row_of::<BlastTable>(&hit("q", "t"))).unwrap();
        let fields: Vec<&str> = text.trim_end().split('\t').collect();

        assert_eq!(fields.len(), BlastTable::N_FIELDS);
        assert_eq!(fields[BlastTable::QUERY], "q");
        assert_eq!(fields[2], "-", "pident is not modelled");
        assert_eq!(fields[3], "-", "length is not modelled");
    }

    #[test]
    fn a_score_and_an_e_value_survive_the_text_they_are_written_as() {
        // float Display is the shortest text that reads back as
        // the same value; if it were not, the round trip above
        // would be checking rounding rather than the layout
        for (score, e_value) in [
            (256.1f32, 3.2e-88f64),
            (0.0, 1.0),
            (-9.75, 4.7e-5),
            (f32::MAX, f64::MIN_POSITIVE),
        ] {
            let h = Hit {
                score,
                e_value,
                ..hit("q", "t")
            };
            let back = HitParser::<BlastTable>::parse(&row_of::<BlastTable>(&h)).unwrap();

            assert_eq!(back.score.to_bits(), score.to_bits(), "{score}");
            assert_eq!(back.e_value.to_bits(), e_value.to_bits(), "{e_value}");
        }
    }

    #[test]
    fn a_row_too_narrow_for_its_layout_says_so() {
        // N_FIELDS catching a table in the wrong layout: a
        // 10-field nail row read as a 12-field BLAST one
        let nail = row_of::<NailTable>(&hit("q", "t"));
        let err = HitParser::<BlastTable>::parse(&nail).unwrap_err();

        assert!(matches!(err.kind(), crate::ErrorKind::Parse { .. }));
        assert!(err.to_string().contains("12 fields"), "{err}");
        assert!(err.to_string().contains("found 10"), "{err}");
    }

    #[test]
    fn a_wider_row_is_read_by_its_first_columns() {
        // why the check is >= and not ==: `-outfmt "6 std qlen
        // slen"` appends columns after the standard twelve, and
        // the twelve in front are still the standard set
        let mut row = row_of::<BlastTable>(&hit("q", "t"));
        row.pop();
        row.extend_from_slice(b"\t128\t203\n");

        assert_eq!(
            HitParser::<BlastTable>::parse(&row).unwrap(),
            hit("q", "t"),
            "the appended columns are ignored, not an error"
        );
    }

    #[test]
    fn a_column_that_is_not_a_number_names_itself() {
        let mut row = row_of::<BlastTable>(&hit("q", "t"));
        let text = String::from_utf8(row.clone()).unwrap();
        row = text.replace("1.1e-12", "not-a-number").into_bytes();

        let err = HitParser::<BlastTable>::parse(&row).unwrap_err();
        assert!(err.to_string().contains("e-value"), "{err}");
    }

    // ---

    const ROWS: &[u8] = b"# a header\nq1\tt1\t-\t-\t-\t-\t1\t10\t5\t14\t1e-9\t12.5\n\
                          q2\tt2\t-\t-\t-\t-\t1\t20\t30\t11\t1e-3\t8\n\
                          #\n# Program: made up\n";

    fn two() -> BlastHits {
        Table::new(ROWS).unwrap()
    }

    fn queries<C: Indexable<Record = Hit>>(c: &C) -> Vec<String> {
        c.iter().map(|h| h.borrow().query.clone()).collect()
    }

    #[test]
    fn a_table_becomes_its_rows_and_no_comments() {
        let t = two();

        assert_eq!(t.len(), 2, "four lines, two of them comments");
        assert!(!t.is_empty());
        assert_eq!(queries(&t), ["q1", "q2"]);
        assert!(t.get(2).is_none(), "past the end is None, not a panic");
        assert_eq!(t.get(1).unwrap().target_span, Some(Span::new(30, 11)));
    }

    #[test]
    fn a_table_of_only_comments_is_an_empty_collection() {
        let t = BlastHits::new(&b"# nothing here\n#\n"[..]).unwrap();
        assert!(t.is_empty());
        assert_eq!(t.iter().count(), 0);
        assert_eq!(BlastHits::new(&b""[..]).unwrap().len(), 0);
    }

    #[test]
    fn a_malformed_row_names_the_row_that_failed() {
        let mut bytes = ROWS.to_vec();
        bytes.extend_from_slice(b"q3\tt3\n");

        let err = BlastHits::new(&bytes[..]).unwrap_err();
        assert!(err.to_string().contains("row 2"), "{err}");
        assert!(err.to_string().contains("byte"), "{err}");
    }

    #[test]
    fn a_whole_table_round_trips_through_its_own_writer() {
        let t = two();
        let mut out = Vec::new();
        t.write_to(&mut out).unwrap();

        let back = BlastHits::new(&out[..]).unwrap();
        assert_eq!(back.len(), t.len());
        for n in 0..t.len() {
            assert_eq!(back.get(n), t.get(n), "row {n}");
        }
    }

    #[test]
    fn the_combinators_work_on_hits() {
        // sorting by e-value is the operation this exists for
        let best = two().sort_by_key(|h| h.e_value.to_bits());
        assert_eq!(queries(&best), ["q1", "q2"], "1e-9 before 1e-3");

        let strong = two().filter(|h| h.e_value < 1e-6);
        assert_eq!(queries(&strong), ["q1"]);
        assert_eq!(two().sample(2).len(), 2);
    }

    // ---

    fn temp_tbl(tag: &str, bytes: &[u8]) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("libsail-{tag}-{}.tsv", std::process::id()));
        std::fs::write(&p, bytes).unwrap();
        p
    }

    #[test]
    fn the_disk_backed_reader_agrees_with_the_in_memory_one() {
        // the claim of Indexable is that one API covers both
        // backends; this is what checks it beyond the types
        let path = temp_tbl("agree", ROWS);
        let mem = two();
        let disk: IndexedTable<BlastTable> = IndexedTable::open(&path).unwrap();

        assert_eq!(disk.len(), mem.len());
        for n in 0..mem.len() {
            assert_eq!(disk.cloned(n), mem.cloned(n), "row {n}");
        }

        assert!(disk.get(mem.len()).is_none(), "past the end is None");
        assert!(disk.try_get(mem.len()).unwrap().is_none());
        disk.validate().unwrap();
        assert_eq!(queries(&disk), ["q1", "q2"]);

        std::fs::remove_file(path).ok();
    }

    #[cfg(feature = "mmap")]
    #[test]
    fn the_mapped_backend_agrees_with_the_read_one() {
        let path = temp_tbl("mapped", ROWS);
        let mem = two();
        // SAFETY: this test owns the file and does not touch
        // it while the map is alive
        let mapped: IndexedTable<BlastTable, MmapSource> =
            unsafe { IndexedTable::map(&path) }.unwrap();

        assert_eq!(mapped.len(), mem.len());
        for n in 0..mem.len() {
            assert_eq!(mapped.cloned(n), mem.cloned(n), "row {n}");
        }
        mapped.validate().unwrap();

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn validate_finds_the_row_get_can_only_call_none() {
        let mut bytes = ROWS.to_vec();
        bytes.extend_from_slice(b"q3\tt3\n");

        let path = temp_tbl("bad", &bytes);
        let disk: IndexedTable<BlastTable> = IndexedTable::open(&path).unwrap();

        assert_eq!(disk.len(), 3);
        assert!(disk.get(0).is_some());
        assert!(disk.get(2).is_none(), "get has nowhere to put the reason");

        let err = disk.validate().unwrap_err();
        assert!(err.to_string().contains("row 2"), "{err}");

        std::fs::remove_file(path).ok();
    }

    // ---

    #[test]
    fn domain_rows_group_by_the_pair_they_describe() {
        // the domain table holds 7 rows over 5 query-target
        // pairs, so grouping has something to group
        let dom = HmmerDomHits::open(fixture!("hits.hmmer.domtblout")).unwrap();

        let groups = dom.group_by(|h| (h.query.clone(), h.target.clone()));

        assert_eq!(groups.len(), 5, "7 rows over 5 pairs");
        assert_eq!(
            groups.iter().map(|(_, g)| g.len()).sum::<usize>(),
            7,
            "every row lands in exactly one group"
        );

        // the tandem-domain proteins are the ones worth finding:
        // Grb2's two SH3 domains and NHERF-1's two PDZ domains
        let tandem: Vec<&str> = groups
            .iter()
            .filter(|(_, g)| g.len() == 2)
            .map(|(k, _)| k.1.as_str())
            .collect();
        assert_eq!(tandem, ["SEM5_CAEEL", "NHRF1_HUMAN"]);

        // and the groups arrive in the order the table lists
        // them rather than however they hash
        assert_eq!(groups[0].0.1, "SEM5_CAEEL", "the first row's target leads");

        // the best hit per pair, which is what a caller
        // actually wants out of a grouping
        let best = groups[0].1.min_by_key(|h| h.e_value.to_bits()).unwrap();
        assert!(best.e_value <= groups[0].1.get(0).unwrap().e_value);
    }

    #[test]
    fn the_real_fixtures_parse_to_the_counts_they_hold() {
        let blast = BlastHits::open(fixture!("hits.blast.tbl")).unwrap();
        let nail = NailHits::open(fixture!("hits.nail.tbl")).unwrap();
        let tblout = HmmerHits::open(fixture!("hits.hmmer.tblout")).unwrap();
        let dom = HmmerDomHits::open(fixture!("hits.hmmer.domtblout")).unwrap();

        // nail reports fewer rows than hmmsearch: its mmseqs
        // seeding does not reach NHERF-1's weaker PDZ, so the
        // two tables are separate searches and not one search
        // counted twice
        assert_eq!(blast.len(), 4);
        assert_eq!(nail.len(), 4);
        assert_eq!(tblout.len(), 5);
        assert_eq!(dom.len(), 7);

        // every hit runs forward: these are protein searches,
        // so no span is reversed the way a nucleotide one can be
        assert!(blast.iter().all(|h| !h.target_span.unwrap().is_reverse()));

        // HmmerTable reports no coordinates, and the free-text
        // description past the last column is not a field
        assert!(tblout.iter().all(|h| h.query_span.is_none()));
        assert_eq!(tblout.get(0).unwrap().query, "SH3_1");
        assert_eq!(tblout.get(0).unwrap().target, "SEM5_CAEEL");
    }

    #[test]
    fn blast_and_nail_describe_the_same_search() {
        // the strongest check the fixtures allow: two layouts
        // over one set of hits must agree once parsed, which no
        // single-format test can catch a transposed column with
        let (blast, nail) = (
            BlastHits::open(fixture!("hits.blast.tbl")).unwrap(),
            NailHits::open(fixture!("hits.nail.tbl")).unwrap(),
        );

        let key = |h: &Hit| (h.query.clone(), h.target.clone());
        let mut a: Vec<Hit> = blast.iter().cloned().collect();
        let mut b: Vec<Hit> = nail.iter().cloned().collect();
        a.sort_by_key(key);
        b.sort_by_key(key);

        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(&b) {
            assert_eq!((&x.query, &x.target), (&y.query, &y.target));
            assert_eq!(x.query_span, y.query_span, "{}/{}", x.query, x.target);
            assert_eq!(x.target_span, y.target_span, "{}/{}", x.query, x.target);
            assert_eq!(x.score, y.score, "{}/{}", x.query, x.target);
        }
    }

    #[test]
    fn every_fixture_row_survives_its_own_writer() {
        let blast = BlastHits::open(fixture!("hits.blast.tbl")).unwrap();

        for (name, bytes) in [
            ("blast", {
                let mut o = Vec::new();
                blast.write_to(&mut o).unwrap();
                o
            }),
            ("dom", {
                let mut o = Vec::new();
                HmmerDomHits::open(fixture!("hits.hmmer.domtblout"))
                    .unwrap()
                    .write_to(&mut o)
                    .unwrap();
                o
            }),
        ] {
            let reread = match name {
                "blast" => {
                    BlastHits::new(&bytes[..]).map(|t| t.iter().cloned().collect::<Vec<_>>())
                }
                _ => HmmerDomHits::new(&bytes[..]).map(|t| t.iter().cloned().collect::<Vec<_>>()),
            }
            .unwrap();

            let want: Vec<Hit> = match name {
                "blast" => blast.iter().cloned().collect(),
                _ => HmmerDomHits::open(fixture!("hits.hmmer.domtblout"))
                    .unwrap()
                    .iter()
                    .cloned()
                    .collect(),
            };

            assert_eq!(reread, want, "{name} did not survive a write and a read");
        }
    }
}

//! Pairwise search result tables.
//!
//! Every tool here reports the same thing — one query matched one target, with
//! a score, an E-value, and optionally where the match sat in each sequence.
//! What differs is which column holds what, so a [`Hit`] is the shape and a
//! [`HitColumns`] impl is one tool's layout of it.
//!
//! That shape is the intersection of the formats. Everything a single tool
//! reports past it — BLAST's `pident`, nail's `cellfrac`, HMMER's accessions,
//! sequence lengths and description — lives in that format's own row type
//! beside its layout: [`blast::BlastRow`], [`nail::NailRow`],
//! [`hmmer::HmmerRow`], [`hmmer::HmmerDomRow`]. Each converts into a [`Hit`].

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

/// One tool's report that a query matched a target, holding what every layout
/// modelled here has in common.
#[derive(Debug, Clone, PartialEq)]
pub struct Hit {
    pub query: String,
    pub target: String,

    /// The first query position covered, 1-based and inclusive, where the
    /// layout reports coordinates.
    //
    // above the end coordinate on a minus-strand match, which
    // is how BLAST and nhmmer spell one -- and the reason
    // these are bare coordinates: a Range would render such a
    // match empty
    pub query_start: Option<usize>,

    /// The last query position covered, which may be below
    /// [`query_start`](Self::query_start).
    pub query_end: Option<usize>,

    /// The first target position covered, 1-based and inclusive, where the
    /// layout reports coordinates.
    pub target_start: Option<usize>,

    /// The last target position covered, which may be below
    /// [`target_start`](Self::target_start).
    pub target_end: Option<usize>,

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
    /// Mandatory columns in a row.
    //
    // a row with fewer is an error rather than a panic, which
    // catches a table in the wrong format at its first row
    const N_COLUMNS: usize;
}

// ---

/// Whether `b` opens a comment, which every layout modelled here spells `#`.
fn is_comment(b: u8) -> bool {
    b == b'#'
}

/// The columns at `want`, borrowed from `row`. `None` where the row holds too
/// few of them.
pub fn fields<const K: usize>(row: &[u8], want: [usize; K]) -> Option<[&[u8]; K]> {
    let mut out = [&row[..0]; K];
    let mut found = [false; K];

    // nothing asked for is nothing to walk the row for, and
    // max over no columns has no answer to compare against
    let Some(&last) = want.iter().max() else {
        return Some(out);
    };

    let mut column = 0;
    let mut at = 0;

    loop {
        while at < row.len() && row[at].is_ascii_whitespace() {
            at += 1;
        }

        if at == row.len() {
            break;
        }

        let start = at;
        while at < row.len() && !row[at].is_ascii_whitespace() {
            at += 1;
        }

        // want may name one column more than once, so this
        // fills every position that asked for it rather than
        // the first
        for (i, &w) in want.iter().enumerate() {
            if w == column {
                out[i] = &row[start..at];
                found[i] = true;
            }
        }

        // the row past the last wanted column is never read,
        // which is most of it where a caller wants three of
        // a domain table's twenty-two
        if column == last {
            break;
        }

        column += 1;
    }

    found.iter().all(|&f| f).then_some(out)
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
        let (cols, _) = split_columns(bytes, C::N_COLUMNS, "a row of this layout")?;

        let (query_start, query_end) = pair(&cols, C::Q_START, C::Q_END, "query")?;
        let (target_start, target_end) = pair(&cols, C::T_START, C::T_END, "target")?;

        Ok(Hit {
            query: column(&cols, C::QUERY, "query")?.to_string(),
            target: column(&cols, C::TARGET, "target")?.to_string(),
            query_start,
            query_end,
            target_start,
            target_end,
            score: number(column(&cols, C::SCORE, "score")?, "score")?,
            e_value: number(column(&cols, C::E_VALUE, "e-value")?, "e-value")?,
        })
    }
}

/// The first `n` columns of a row, plus the rest of the line unsplit — which
/// is where HMMER writes its free-text description.
///
/// # Errors
///
/// A row of fewer than `n` columns, reported as `what`.
pub(crate) fn split_columns<'a>(
    bytes: &'a [u8],
    n: usize,
    what: &str,
) -> Result<(Vec<&'a str>, &'a str)> {
    let Ok(text) = std::str::from_utf8(bytes) else {
        bail_parse!("a row is utf-8 text")
    };

    let mut cols = Vec::with_capacity(n);
    let mut rest = text.trim_end();

    for _ in 0..n {
        rest = rest.trim_start();

        // any whitespace and not tabs: BLAST writes tabs, nail
        // and HMMER write aligned runs of spaces, and none of
        // them permit a space inside a modelled column
        let end = rest.find(char::is_whitespace).unwrap_or(rest.len());
        if end == 0 {
            break;
        }

        cols.push(&rest[..end]);
        rest = &rest[end..];
    }

    // narrower than the format is an error and wider is not:
    // HMMER's description is free text past the last column,
    // and `-outfmt "6 std qlen slen"` appends columns after
    // the standard twelve. stopping the loop at n is what
    // reads those by their first n columns
    if cols.len() < n {
        bail_parse!(
            "{what} has {n} columns, found {}: {:?}",
            cols.len(),
            text.trim_end()
        )
    }

    Ok((cols, rest.trim()))
}

/// Column `i`, or a failure naming the column rather than the index.
fn column<'a>(cols: &[&'a str], i: usize, what: &str) -> Result<&'a str> {
    match cols.get(i) {
        Some(c) => Ok(*c),
        // reachable only from a layout whose index is past
        // its own N_COLUMNS, which
        // every_layout_indexes_within_its_own_row rules out
        // for the four here but not for a caller's
        None => bail_parse!("this layout puts {what} at column {i}, past the row"),
    }
}

/// Both coordinates of a span, or `(None, None)` where the layout reports
/// neither.
fn pair(
    cols: &[&str],
    start: Option<usize>,
    end: Option<usize>,
    what: &str,
) -> Result<(Option<usize>, Option<usize>)> {
    match (start, end) {
        (Some(s), Some(e)) => Ok((
            Some(number(column(cols, s, what)?, what)?),
            Some(number(column(cols, e, what)?, what)?),
        )),
        (None, None) => Ok((None, None)),

        // a layout naming one end and not the other would
        // otherwise drop the coordinate it does report
        _ => bail_parse!("this layout reports half a {what} span"),
    }
}

pub(crate) fn number<T>(value: &str, what: &str) -> Result<T>
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
        let mut cols = vec!["-".to_string(); C::N_COLUMNS];

        let mut set = |i: usize, value: String| -> Result<()> {
            match cols.get_mut(i) {
                Some(c) => {
                    *c = value;
                    Ok(())
                }
                None => bail_parse!(
                    "this layout puts a column at {i}, past its own {} columns",
                    C::N_COLUMNS
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

        for (column, coordinate) in [
            (C::Q_START, self.query_start),
            (C::Q_END, self.query_end),
            (C::T_START, self.target_start),
            (C::T_END, self.target_end),
        ] {
            // a coordinate is dropped where the layout has
            // nowhere to put it -- HmmerTable reports none at
            // all -- so a round trip through one needs a Hit
            // that layout can hold
            if let (Some(i), Some(c)) = (column, coordinate) {
                set(i, c.to_string())?;
            }
        }

        writeln!(w, "{}", cols.join("\t"))?;

        Ok(())
    }
}

// ---

/// Rendering one row, the reverse of the format's [`Parse`] impl.
pub trait WriteRow: Parse {
    /// Write `record` as the row a parser of this format reads back.
    //
    // no receiver, and on the parser rather than on the
    // record, because writing a Hit needs a layout and only
    // the parser has one. that is also what keeps Table to a
    // single write_to: an inherent block over
    // Table<HitParser<C>> beside a generic one over Table<P>
    // collides on the method name, where clauses or not
    fn write_row<W: Write>(record: &Self::Record, w: W) -> Result<()>;
}

impl<C: HitColumns> WriteRow for HitParser<C> {
    fn write_row<W: Write>(record: &Hit, w: W) -> Result<()> {
        record.write_to::<C, _>(w)
    }
}

/// Every row of a hit table, held in memory.
pub struct Table<P: Parse> {
    /// The parsed rows, in file order.
    //
    // P::Record rather than Hit, so that one collection
    // covers both readings of a table: the cross-format hit
    // and the format's own row type. every format here is
    // framed the same way -- one row per non-comment line --
    // so the parser is all that varies
    rows: Vec<P::Record>,

    // fn() -> P rather than P: the format is a parsing
    // decision and not something the collection owns, and
    // this keeps the auto traits without asking P for them
    format: PhantomData<fn() -> P>,
}

/// BLAST `-outfmt 6` rows as cross-format [`Hit`]s, in memory.
pub type BlastHits = Table<HitParser<crate::tbl::blast::BlastTable>>;
/// nail rows as cross-format [`Hit`]s, in memory.
pub type NailHits = Table<HitParser<crate::tbl::nail::NailTable>>;
/// HMMER `--tblout` rows as cross-format [`Hit`]s, in memory.
pub type HmmerHits = Table<HitParser<crate::tbl::hmmer::HmmerTable>>;
/// HMMER `--domtblout` rows as cross-format [`Hit`]s, in memory.
pub type HmmerDomHits = Table<HitParser<crate::tbl::hmmer::HmmerDomTable>>;

/// BLAST `-outfmt 6` rows with every column, in memory.
pub type BlastRows = Table<crate::tbl::blast::BlastRow>;
/// nail rows with every column, in memory.
pub type NailRows = Table<crate::tbl::nail::NailRow>;
/// HMMER `--tblout` rows with every column, in memory.
pub type HmmerRows = Table<crate::tbl::hmmer::HmmerRow>;
/// HMMER `--domtblout` rows with every column, in memory.
pub type HmmerDomRows = Table<crate::tbl::hmmer::HmmerDomRow>;

impl<P: Parse> Table<P> {
    /// Read and parse every row from `bytes`.
    pub fn new<R: Read>(mut bytes: R) -> Result<Self> {
        let mut all = Vec::new();
        bytes.read_to_end(&mut all)?;

        let rows = line_index(&all[..], is_comment)?
            .iter()
            .enumerate()
            .map(|(i, o)| {
                let range = &all[o.start as usize..][..o.n_bytes];

                // which row failed matters as much as why: a
                // table in the wrong layout fails at its first
                // one, and the byte offset is how to find it
                P::parse(range)
                    .map_err(|e| Error::parse(format!("row {i} at byte {}: {e}", o.start)))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Table {
            rows,
            format: PhantomData,
        })
    }

    /// Read and parse every row in the file at `path`.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Table::new(File::open(path)?)
    }
}

impl<P: WriteRow> Table<P> {
    /// Write every row in this format's column order. See [`Hit::write_to`]
    /// for what a row of cross-format hits can and cannot carry back.
    pub fn write_to<W: Write>(&self, mut w: W) -> Result<()> {
        for row in &self.rows {
            P::write_row(row, &mut w)?;
        }

        Ok(())
    }
}

// written out rather than derived: a derive bounds every
// type parameter, so it would demand a parser be Debug and
// Clone to describe a table never holding one
impl<P: Parse> std::fmt::Debug for Table<P>
where
    P::Record: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Table").field("rows", &self.rows).finish()
    }
}

impl<P: Parse> Clone for Table<P> {
    fn clone(&self) -> Self {
        Table {
            rows: self.rows.clone(),
            format: PhantomData,
        }
    }
}

impl<P: Parse> PartialEq for Table<P>
where
    P::Record: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.rows == other.rows
    }
}

impl<P: Parse> Default for Table<P> {
    fn default() -> Self {
        Table {
            rows: Vec::new(),
            format: PhantomData,
        }
    }
}

impl<P: Parse> Indexable for Table<P> {
    type Record = P::Record;

    // the rows are right here, so lend one
    type Yield<'a>
        = &'a P::Record
    where
        Self: 'a;

    fn len(&self) -> usize {
        self.rows.len()
    }

    fn get(&self, n: usize) -> Option<&P::Record> {
        self.rows.get(n)
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
pub struct IndexedTable<P: Parse, S: Source = FileSource> {
    source: S,
    offsets: Vec<Offset>,
    format: PhantomData<fn() -> P>,
}

#[cfg(any(unix, windows))]
impl<P: Parse> IndexedTable<P, FileSource> {
    /// Index the file at `path`, parsing none of it.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let source = FileSource::open(path)?;
        let offsets = line_index(source.file(), is_comment)?;

        Ok(IndexedTable {
            source,
            offsets,
            format: PhantomData,
        })
    }
}

#[cfg(feature = "mmap")]
impl<P: Parse> IndexedTable<P, MmapSource> {
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
            format: PhantomData,
        })
    }
}

#[cfg(any(unix, windows))]
impl<P: Parse, S: Source> IndexedTable<P, S> {
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
    pub fn try_get(&self, n: usize) -> Result<Option<P::Record>> {
        let Some(o) = self.offsets.get(n) else {
            return Ok(None);
        };

        Ok(Some(P::parse(&self.source.range(o.start, o.n_bytes)?)?))
    }
}

#[cfg(any(unix, windows))]
impl<P: Parse, S: Source> Indexable for IndexedTable<P, S> {
    type Record = P::Record;

    // nothing to lend: the row is parsed inside the call and
    // did not exist before it
    type Yield<'a>
        = P::Record
    where
        Self: 'a;

    fn len(&self) -> usize {
        self.offsets.len()
    }

    /// Reads and parses the row's bytes on every call. A failure comes back as
    /// `None`; [`try_get`](Self::try_get) keeps the reason.
    fn get(&self, n: usize) -> Option<P::Record> {
        self.try_get(n).ok().flatten()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixture;
    use crate::tbl::blast::{BlastRow, BlastTable};
    use crate::tbl::hmmer::{HmmerDomRow, HmmerDomTable, HmmerRow, HmmerTable};
    use crate::tbl::nail::{NailRow, NailTable};

    #[test]
    fn hit_display_is_a_summary_and_says_so() {
        // Display is a summary line, not a rendering of the
        // record: unlike the seq types it drops the score and
        // both spans, and cannot round-trip
        let hit = Hit {
            query: "q1".into(),
            target: "t1".into(),
            query_start: Some(1),
            query_end: Some(10),
            target_start: Some(5),
            target_end: Some(14),
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
        let n = C::N_COLUMNS;
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
        assert_ne!(BlastTable::N_COLUMNS, HmmerTable::N_COLUMNS);
        assert_ne!(HmmerTable::N_COLUMNS, HmmerDomTable::N_COLUMNS);
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
            query_start: Some(12),
            query_end: Some(101),
            target_start: Some(88),
            target_end: Some(21),
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
    fn a_backwards_pair_of_coordinates_is_a_minus_strand_match() {
        // why the coordinates are bare and not a Range: BLAST
        // and nhmmer report a minus-strand match with start
        // above end, which a Range renders empty. it must
        // survive a round trip in the order it arrived
        let h = Hit {
            target_start: Some(900),
            target_end: Some(400),
            ..hit("q", "t")
        };
        let back = HitParser::<BlastTable>::parse(&row_of::<BlastTable>(&h)).unwrap();

        assert_eq!(back, h);
        assert!(back.target_start > back.target_end, "still reversed");
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
                query_start: Some(1),
                query_end: Some(1),
                target_start: Some(9999),
                target_end: Some(1),
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
            query_start: None,
            query_end: None,
            target_start: None,
            target_end: None,
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
        let cols: Vec<&str> = text.trim_end().split('\t').collect();

        assert_eq!(cols.len(), BlastTable::N_COLUMNS);
        assert_eq!(cols[BlastTable::QUERY], "q");
        assert_eq!(cols[2], "-", "pident is not modelled");
        assert_eq!(cols[3], "-", "length is not modelled");
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
        // N_COLUMNS catching a table in the wrong layout: a
        // 10-column nail row read as a 12-column BLAST one
        let nail = row_of::<NailTable>(&hit("q", "t"));
        let err = HitParser::<BlastTable>::parse(&nail).unwrap_err();

        assert!(matches!(err.kind(), crate::ErrorKind::Parse { .. }));
        assert!(err.to_string().contains("12 columns"), "{err}");
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
        assert_eq!(t.get(1).unwrap().target_start, Some(30));
        assert_eq!(t.get(1).unwrap().target_end, Some(11));
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
        let disk: IndexedTable<HitParser<BlastTable>> = IndexedTable::open(&path).unwrap();

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
        let mapped: IndexedTable<HitParser<BlastTable>, MmapSource> =
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
        let disk: IndexedTable<HitParser<BlastTable>> = IndexedTable::open(&path).unwrap();

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

        // the best hit per pair, which is what a
        // grouping is for
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
        assert!(
            blast
                .iter()
                .all(|h| h.target_start.unwrap() < h.target_end.unwrap())
        );

        // HmmerTable reports no coordinates, and the free-text
        // description past the last column is not a field
        assert!(tblout.iter().all(|h| h.query_start.is_none()));
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
            assert_eq!(
                (x.query_start, x.query_end),
                (y.query_start, y.query_end),
                "{}/{}",
                x.query,
                x.target
            );
            assert_eq!(
                (x.target_start, x.target_end),
                (y.target_start, y.target_end),
                "{}/{}",
                x.query,
                x.target
            );
            assert_eq!(x.score, y.score, "{}/{}", x.query, x.target);
        }
    }

    /// Every row of `path`, read through the format's own row type and through
    /// its [`HitColumns`] layout, must describe the same hit.
    fn rows_agree<R, C>(path: &str)
    where
        R: Parse<Record = R>,
        C: HitColumns,
        for<'a> &'a R: Into<Hit>,
    {
        let rows: Table<R> = Table::open(path).unwrap();
        let hits: Table<HitParser<C>> = Table::open(path).unwrap();

        assert!(!rows.is_empty(), "{path} has rows to compare");
        assert_eq!(rows.len(), hits.len(), "{path}");

        for n in 0..rows.len() {
            let converted: Hit = rows.get(n).unwrap().into();
            assert_eq!(&converted, hits.get(n).unwrap(), "row {n} of {path}");
        }
    }

    #[test]
    fn the_full_row_and_the_generic_hit_agree_on_every_fixture_row() {
        // a format's row type and its HitColumns layout are two
        // independent readings of the same file. over real tool
        // output, agreement is what catches a column either one
        // has transposed, and nothing else here does
        rows_agree::<BlastRow, BlastTable>(fixture!("hits.blast.tbl"));
        rows_agree::<NailRow, NailTable>(fixture!("hits.nail.tbl"));
        rows_agree::<HmmerRow, HmmerTable>(fixture!("hits.hmmer.tblout"));
        rows_agree::<HmmerDomRow, HmmerDomTable>(fixture!("hits.hmmer.domtblout"));
    }

    #[test]
    fn the_collections_read_a_full_row_type_on_every_backend() {
        // Table and IndexedTable are generic over the parser
        // and not over a hit layout, so the same pair has to
        // carry a format's own row type
        let path = fixture!("hits.hmmer.domtblout");
        let mem = HmmerDomRows::open(path).unwrap();
        let disk: IndexedTable<HmmerDomRow> = IndexedTable::open(path).unwrap();

        assert_eq!(disk.len(), mem.len());
        for n in 0..mem.len() {
            assert_eq!(disk.cloned(n), mem.cloned(n), "row {n}");
        }

        disk.validate().unwrap();

        #[cfg(feature = "mmap")]
        {
            // SAFETY: a fixture is read-only test data, and
            // nothing here writes to it while the map is alive
            let mapped: IndexedTable<HmmerDomRow, MmapSource> =
                unsafe { IndexedTable::map(path) }.unwrap();

            assert_eq!(mapped.len(), mem.len());
            for n in 0..mem.len() {
                assert_eq!(mapped.cloned(n), mem.cloned(n), "row {n}");
            }
        }

        // and the combinators come with them
        let tandem = mem.filter(|r| r.domain_count == 2);
        assert_eq!(tandem.len(), 4, "two proteins with two domains each");
    }

    /// Every row of `path` written back and read again must be unchanged.
    fn rows_survive_a_write<R>(path: &str)
    where
        R: Parse<Record = R> + WriteRow + PartialEq + std::fmt::Debug,
    {
        let rows: Table<R> = Table::open(path).unwrap();
        assert!(!rows.is_empty(), "{path} has rows to write");

        let mut out = Vec::new();
        rows.write_to(&mut out).unwrap();

        assert_eq!(Table::<R>::new(&out[..]).unwrap(), rows, "{path}");
    }

    #[test]
    fn a_whole_table_of_full_rows_survives_its_own_writer() {
        // every_fixture_row_survives_its_own_writer can only
        // check the six values a Hit models, since the writer
        // fills the rest with '-'. a row type holds every
        // column, so here the whole file comes back
        rows_survive_a_write::<BlastRow>(fixture!("hits.blast.tbl"));
        rows_survive_a_write::<NailRow>(fixture!("hits.nail.tbl"));
        rows_survive_a_write::<HmmerRow>(fixture!("hits.hmmer.tblout"));
        rows_survive_a_write::<HmmerDomRow>(fixture!("hits.hmmer.domtblout"));
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

    // ---

    #[test]
    fn fields_borrows_the_columns_asked_for() {
        let row = b"query1  target1  42.5  1e-10";

        let [q, e] = fields(row, [0, 3]).expect("four columns are there");
        assert_eq!(q, b"query1");
        assert_eq!(e, b"1e-10");

        let [only] = fields(row, [2]).expect("column 2 is there");
        assert_eq!(only, b"42.5");
    }

    #[test]
    fn a_run_of_whitespace_separates_one_column_from_the_next() {
        // hmmer pads its columns to align them, so a row
        // carries runs of spaces that are one separator
        let row = b"  query1 \t\t target1    42.5  ";

        let [a, b, c] = fields(row, [0, 1, 2]).expect("three columns are there");
        assert_eq!((a, b, c), (&b"query1"[..], &b"target1"[..], &b"42.5"[..]));
    }

    #[test]
    fn a_row_with_too_few_columns_is_refused_rather_than_padded() {
        // a table of the wrong tool reads as a short row, and
        // answering with an empty field would score it
        let row = b"query1 target1";

        assert!(fields(row, [0, 1]).is_some());
        assert!(fields(row, [0, 5]).is_none());
        assert!(fields(b"", [0]).is_none());
    }

    #[test]
    fn a_column_asked_for_twice_is_answered_twice() {
        let row = b"a b c";

        let [x, y, z] = fields(row, [1, 1, 0]).expect("both columns are there");
        assert_eq!((x, y, z), (&b"b"[..], &b"b"[..], &b"a"[..]));
    }
}

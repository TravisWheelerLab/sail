//! Byte offsets of the records in a file.
//!
//! Records in these formats are separated by a fixed pattern at the start of a
//! line — FASTA's `>` opens the record it introduces, Stockholm's `//` closes
//! the one it ends.
//!
//! `scan` reads bytes and reports delimiter positions, which is where the
//! chunk-boundary care lives. `frame` turns those into [`Offset`]s, which is
//! pure arithmetic and the only part that consults [`Delimiter`].

use std::borrow::Cow;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::num::NonZeroU64;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use memchr::memmem;

use crate::error::{Error, Result, bail_parse};
use crate::format::Format;
#[cfg(feature = "mmap")]
use crate::source::MmapSource;
#[cfg(any(unix, windows))]
use crate::source::{FileSource, Source};

/// Where one record sits in the bytes it was indexed from.
#[derive(Default, Clone, Copy, Debug, PartialEq)]
pub struct Offset {
    /// Absolute offset of the record's first byte.
    pub start: u64,
    /// Length of the record in bytes.
    pub n_bytes: usize,
}

/// Whether a delimiter opens the record after it or closes the record before
/// it.
pub enum Delimiter {
    /// FASTA's `>`: the delimiter is the first byte of its own record.
    Initiating,
    /// Stockholm's and p7hmm's `//`: the delimiter is the last line of its own
    /// record.
    Terminating,
}

/// A format whose records are separated by a fixed pattern at the start of a
/// line. An impl supplies data only; everything else here is derived from it.
pub trait Delimit {
    /// The record delimiter.
    const DELIM: &'static [u8];

    /// Whether [`DELIM`](Self::DELIM) opens or closes the record it marks.
    const DELIM_TYPE: Delimiter;

    /// Width of the anchored pattern actually searched for: a newline plus
    /// [`DELIM`](Self::DELIM).
    //
    // derived, with CARRY below, so the two cannot drift
    // apart and so an impl can neither forget the anchor nor
    // apply it twice
    const WIDTH: usize = 1 + Self::DELIM.len();

    /// Bytes carried across a chunk boundary.
    //
    // forced rather than tuned: one less than the pattern is
    // the longest run that cannot contain a whole pattern
    const CARRY: usize = Self::WIDTH - 1;

    /// Fails to compile for a delimiter `scan` cannot handle.
    const VALID: () = {
        // CARRY of 0 carries nothing, so buf[0] -- the
        // seeded newline standing in for the line ending
        // before byte 0 -- is lost on the first refill
        assert!(
            !Self::DELIM.is_empty(),
            "a delimiter must be at least one byte"
        );

        // a newline in DELIM gives the anchored needle a
        // proper border, since a border of length k needs
        // needle[WIDTH - k] == needle[0] == '\n'. memmem
        // reports non-overlapping matches, which agree with
        // testing every position only when no border exists
        //
        // slice methods are not const, hence the bare while
        let mut i = 0;
        while i < Self::DELIM.len() {
            assert!(
                Self::DELIM[i] != b'\n',
                "a delimiter cannot contain a newline"
            );
            i += 1;
        }
    };
}

// ---

/// Fresh bytes read per round of [`scan`].
const CHUNK: usize = 64 * 1024;

/// Fill `buf` as far as repeated reads will take it, returning 0 only at EOF.
pub(crate) fn read_some<R: Read>(bytes: &mut R, buf: &mut [u8]) -> Result<usize> {
    let mut filled = 0;

    // Read::read may stop short without being at the end, so
    // one call cannot tell a short read from exhaustion
    while filled < buf.len() {
        match bytes.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e.into()),
        }
    }

    Ok(filled)
}

/// Absolute positions of every anchored delimiter in `bytes`, and the total
/// number of bytes read.
fn scan<R, D>(mut bytes: R, chunk: usize) -> Result<(Vec<u64>, u64)>
where
    R: Read,
    D: Delimit,
{
    // fresh bytes per round: chunk changes how often the
    // buffer refills and nothing else, per the carry
    // argument below
    let chunk = chunk.max(1);

    // naming VALID is what evaluates it; see Delimit::DELIM
    let () = D::VALID;

    // the anchor is added here, never by the impl: without
    // it the // inside a Pfam #=GF WK http://... line ends
    // the record halfway through
    let mut needle = Vec::with_capacity(D::WIDTH);
    needle.push(b'\n');
    needle.extend_from_slice(D::DELIM);
    let finder = memmem::Finder::new(&needle);

    // buf holds the carry from the previous round followed
    // by the bytes just read:
    //
    //         0        CARRY             filled
    //         |          |                 |
    //   buf   [ c c c c c  f f f f f f f f ]
    //           carried     freshly read
    //
    // invariant, held across every round:
    //     base  = file offset of buf[1]
    //     buf[i] is at file offset base + i - 1
    //
    // **note: buf[0] is a newline that is not in the file. it
    //         stands in for the line ending before byte 0, so
    //         a file opening with a delimiter needs no first
    //         round special case
    let mut buf = vec![0u8; D::CARRY + chunk];
    buf[0] = b'\n';

    let mut filled = 1;
    let mut base = 0u64;
    let mut total = 0u64;
    let mut positions = Vec::new();

    loop {
        let n = read_some(&mut bytes, &mut buf[filled..])?;
        if n == 0 {
            break;
        }

        filled += n;
        total += n as u64;

        // note: a delimiter need not *end* its line -- this
        //       matches "\n//" wherever it sits, so
        //       "//short ACGT" terminates a record. Easel
        //       ends the alignment there too, and a stricter
        //       rule would need lookahead CARRY cannot have
        //
        // no bound of its own: find_iter reports only starts
        // where the whole needle fits in buf[..filled], which
        // is the bound the carry argument below assumes. a
        // needle running off the end is undecided until the
        // next round, which the carry guarantees it reaches
        for i in finder.find_iter(&buf[..filled]) {
            // the match starts at the anchor, so the delimiter
            // is at buf[i + 1], at file offset base + i
            positions.push(base + i as u64);
        }

        // why CARRY = WIDTH - 1 reports every delimiter
        // exactly once, writing p for a match's file offset:
        //
        //   scannable this round:
        //     p - base <= filled - WIDTH
        //   next round's origin:
        //     base' = base + filled - CARRY
        //
        // so a match already reported lands at:
        //   p - base' <= (filled - WIDTH) - (filled - CARRY)
        //              = CARRY - WIDTH
        //              = -1
        //     ** behind the buffer, so never reported twice
        //
        // and one deferred this round lands at:
        //   p - base' > CARRY - WIDTH = -1, so >= 0
        //     ** still in the buffer, so never dropped
        let keep = filled.min(D::CARRY);
        buf.copy_within(filled - keep..filled, 0);
        base += (filled - keep) as u64;
        filled = keep;
    }

    Ok((positions, total))
}

// ---

/// Turn delimiter positions into one [`Offset`] per record.
fn frame<D: Delimit>(positions: &[u64], total: u64) -> Vec<Offset> {
    let span = |start: u64, end: u64| Offset {
        start,
        n_bytes: (end - start) as usize,
    };

    match D::DELIM_TYPE {
        // the delimiter opens its record, so a record runs to
        // the next delimiter, or to the end for the last one.
        // bytes before the first belong to no record
        Delimiter::Initiating => positions
            .iter()
            .zip(positions.iter().skip(1).chain(std::iter::once(&total)))
            .map(|(&start, &end)| span(start, end))
            .collect(),

        // the delimiter closes its record, so a record runs
        // from the end of the previous one through the
        // delimiter's own line. bytes after the last are a
        // truncated record and belong to none
        Delimiter::Terminating => {
            let mut start = 0u64;
            positions
                .iter()
                .map(|&pos| {
                    // + 1 for the newline ending the delimiter's
                    // line; a final delimiter with no trailing
                    // newline clamps to the end
                    let end = (pos + D::DELIM.len() as u64 + 1).min(total);
                    let offset = span(start, end);
                    start = end;
                    offset
                })
                .collect()
        }
    }
}

// ---

/// Index every record in `bytes`, in file order. `R: Read` is enough — the
/// position is tracked by hand, so a pipe indexes the same as a file.
pub fn build_index<R, D>(bytes: R) -> Result<Vec<Offset>>
where
    R: Read,
    D: Delimit,
{
    let (positions, total) = scan::<R, D>(bytes, CHUNK)?;

    Ok(frame::<D>(&positions, total))
}

// ---

/// The extension [`path_for`] appends.
pub const EXTENSION: &str = "saidx";

const MAGIC: [u8; 4] = *b"SLIX";

/// The on-disk layout written here, and the only one read back.
//
// there is no cross-version promise: read_from rejects
// anything else rather than guessing at a layout, which is
// what makes bumping this a safe way to change what follows
// the header
const VERSION: u8 = 1;

/// Longest record name stored.
//
// bounds the allocation a read makes from a length nothing
// has yet corroborated -- a corrupt span would otherwise ask
// for 4 GiB before the read that fails
const MAX_NAME: usize = 64 * 1024;

/// Where every record of one file sits, and what each is called.
///
/// [`new`](Index::new) holds only offsets; [`name_with`](Index::name_with)
/// adds the names, and [`write`](Index::write) saves the pair, stamped with the
/// file they describe.
#[derive(Debug)]
pub struct Index {
    format: Format,
    body: Body,
}

/// Bytes before the first offset: the header, and nothing else.
const HEADER: u64 = 35;

/// One offset on disk: a start and a length, both little-endian.
const OFFSET: u64 = 16;

/// One name span on disk: a present flag, a start, and a length.
const SPAN: u64 = 13;

/// Where an index keeps what it addresses.
//
// an enum rather than a type parameter: a second axis on
// IndexedFasta would show in every signature and could not
// be taken back, and this one is private, so a third
// backing later breaks no caller
#[derive(Debug)]
enum Body {
    /// Read whole, by [`Index::read`].
    Memory {
        offsets: Vec<Offset>,

        /// Every record's name, absent until
        /// [`name_with`](Index::name_with) runs.
        //
        // a layer rather than a field per record: an index
        // is built to address a file, and naming it costs a
        // pass over every record's bytes that most callers
        // never want
        names: Option<Names>,
    },

    /// Left in the file, by [`Index::open`], and read an entry at a time.
    File { file: File, layout: Layout },
}

/// Where each part of an index sits in the file holding it.
#[derive(Debug)]
struct Layout {
    count: usize,
    names: Option<NameLayout>,
}

/// The name layer's three runs, once the blob's length is known.
#[derive(Debug)]
struct NameLayout {
    blob_at: u64,
    spans_at: u64,
    by_name_at: u64,
    by_name_len: usize,
}

impl Body {
    fn len(&self) -> usize {
        match self {
            Body::Memory { offsets, .. } => offsets.len(),
            Body::File { layout, .. } => layout.count,
        }
    }

    fn offset(&self, n: usize) -> Result<Option<Offset>> {
        match self {
            Body::Memory { offsets, .. } => Ok(offsets.get(n).copied()),

            Body::File { file, layout } => {
                if n >= layout.count {
                    return Ok(None);
                }

                let mut buf = [0u8; OFFSET as usize];
                crate::pread::read_exact_at(file, &mut buf, HEADER + n as u64 * OFFSET)?;

                Ok(Some(Offset {
                    start: u64::from_le_bytes(buf[..8].try_into().expect("eight bytes")),
                    n_bytes: narrow(u64::from_le_bytes(
                        buf[8..].try_into().expect("eight bytes"),
                    ))?,
                }))
            }
        }
    }

    /// The record the `i`th entry of the name ordering points at.
    fn ordered(&self, i: usize) -> Result<usize> {
        let Body::File { file, layout } = self else {
            return Ok(0);
        };
        let names = layout
            .names
            .as_ref()
            .expect("an ordering means a name layer");

        let mut buf = [0u8; 4];
        crate::pread::read_exact_at(file, &mut buf, names.by_name_at + i as u64 * 4)?;

        Ok(u32::from_le_bytes(buf) as usize)
    }

    fn has_names(&self) -> bool {
        match self {
            Body::Memory { names, .. } => names.is_some(),
            Body::File { layout, .. } => layout.names.is_some(),
        }
    }

    fn name(&self, n: usize) -> Result<Option<Cow<'_, [u8]>>> {
        match self {
            Body::Memory { names, .. } => {
                Ok(names.as_ref().and_then(|ns| ns.name(n)).map(Cow::Borrowed))
            }

            Body::File { file, layout } => {
                let Some(names) = &layout.names else {
                    return Ok(None);
                };
                if n >= layout.count {
                    return Ok(None);
                }

                let mut buf = [0u8; SPAN as usize];
                crate::pread::read_exact_at(file, &mut buf, names.spans_at + n as u64 * SPAN)?;

                // a record with no name of its own, kept
                // distinct from one named the empty string
                if buf[0] == 0 {
                    return Ok(None);
                }

                let start = u64::from_le_bytes(buf[1..9].try_into().expect("eight bytes"));
                let len = u32::from_le_bytes(buf[9..].try_into().expect("four bytes"));

                let mut name = vec![0u8; len as usize];
                crate::pread::read_exact_at(file, &mut name, names.blob_at + start)?;

                Ok(Some(Cow::Owned(name)))
            }
        }
    }
}

/// The name layer: the names themselves, and an ordering of them for lookup.
#[derive(Clone, Debug, PartialEq)]
struct Names {
    /// Every name, concatenated.
    blob: Vec<u8>,

    /// Where each record's name sits in `blob`, one entry per record.
    //
    // Option, so absence stays distinct from an empty name:
    // a Stockholm alignment with no #=GF ID has no name at
    // all, and must not answer to ""
    spans: Vec<Option<Span>>,

    /// The named records' positions, ordered by their names.
    //
    // a sorted permutation rather than a map: position is a
    // partition_point over it, and it writes out as a flat
    // run of u32
    by_name: Vec<u32>,
}

/// Where one name sits in [`Names::blob`].
//
// one word per record, where a u64 and a u32 under an
// Option take 24 bytes:
//
//     bit  63       set, so the word is never zero and
//                   Option<Span> stays 8 bytes
//     bits 17..63   start, below MAX_BLOB
//     bits  0..17   len, at most MAX_NAME
//
#[derive(Clone, Copy, Debug, PartialEq)]
struct Span(NonZeroU64);

const LEN_BITS: u32 = 17;

/// The most bytes of names one index holds.
const MAX_BLOB: u64 = 1 << 46;

const _: () = assert!(MAX_NAME < 1 << LEN_BITS);

impl Span {
    fn new(start: u64, len: u32) -> Span {
        debug_assert!(start + len as u64 <= MAX_BLOB && len as usize <= MAX_NAME);

        Span(NonZeroU64::new(1 << 63 | start << LEN_BITS | len as u64).expect("bit 63 is set"))
    }

    fn start(self) -> u64 {
        (self.0.get() & !(1 << 63)) >> LEN_BITS
    }

    fn len(self) -> u32 {
        (self.0.get() & ((1 << LEN_BITS) - 1)) as u32
    }
}

/// The name `span` frames.
fn name_at(blob: &[u8], span: Span) -> &[u8] {
    &blob[span.start() as usize..][..span.len() as usize]
}

/// The name `span` frames, or `None` where a record carries none.
fn name_in(blob: &[u8], span: Option<Span>) -> Option<&[u8]> {
    Some(name_at(blob, span?))
}

/// A name layer built one record at a time, in file order.
#[derive(Default)]
struct NameBuilder {
    blob: Vec<u8>,
    spans: Vec<Option<Span>>,
}

impl NameBuilder {
    fn with_capacity(records: usize) -> NameBuilder {
        NameBuilder {
            blob: Vec::new(),
            spans: Vec::with_capacity(records),
        }
    }

    /// Add the next record's name, `None` for a record that carries none.
    fn push(&mut self, name: Option<&[u8]>) -> Result<()> {
        let Some(name) = name else {
            self.spans.push(None);

            return Ok(());
        };

        if name.len() > MAX_NAME {
            bail_parse!(
                "a record name of {} bytes exceeds the {MAX_NAME} an index holds",
                name.len()
            )
        }

        let start = self.blob.len() as u64;
        if start + name.len() as u64 > MAX_BLOB {
            bail_parse!("record names past {MAX_BLOB} bytes exceed what an index holds")
        }

        self.spans.push(Some(Span::new(start, name.len() as u32)));
        self.blob.extend_from_slice(name);

        Ok(())
    }

    fn finish(self) -> Names {
        let NameBuilder { blob, spans } = self;

        // each entry carries its span, so a comparison
        // reads the blob straight away. looking the span up
        // through spans[n] cost a second cache miss per
        // side: 0.70 s against 0.23 s on 2.4M names in
        // shuffled order
        let mut keyed: Vec<(Span, u32)> = spans
            .iter()
            .enumerate()
            .filter_map(|(n, span)| Some(((*span)?, n as u32)))
            .collect();

        // unstable, with the position breaking ties, which is
        // the stable order: equal names stay in file order
        // and position answers with the first of a repeat --
        // the rule collection::Indexable::index_by follows
        keyed.sort_unstable_by(|&(a, m), &(b, n)| {
            name_at(&blob, a).cmp(name_at(&blob, b)).then(m.cmp(&n))
        });

        // copied out rather than collected in place, which
        // would keep keyed's 16-byte slots as by_name's
        // allocation
        let mut by_name = Vec::with_capacity(keyed.len());
        by_name.extend(keyed.iter().map(|&(_, n)| n));

        Names {
            blob,
            spans,
            by_name,
        }
    }
}

impl Names {
    fn name(&self, n: usize) -> Option<&[u8]> {
        name_in(&self.blob, *self.spans.get(n)?)
    }
}

impl Index {
    pub fn new(format: Format, offsets: Vec<Offset>) -> Self {
        Index {
            format,
            body: Body::Memory {
                offsets,
                names: None,
            },
        }
    }

    pub fn format(&self) -> Format {
        self.format
    }

    pub fn len(&self) -> usize {
        self.body.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Where the `n`th record sits, `None` past the end.
    ///
    /// # Errors
    ///
    /// A read failure, for an index left in its file by [`open`](Index::open).
    pub fn offset(&self, n: usize) -> Result<Option<Offset>> {
        self.body.offset(n)
    }

    /// Every offset at once, `None` for an index left in its file.
    //
    // the branch a walk over every record would otherwise
    // pay per record, paid once
    pub fn offsets(&self) -> Option<&[Offset]> {
        match &self.body {
            Body::Memory { offsets, .. } => Some(offsets),
            Body::File { .. } => None,
        }
    }

    pub fn has_names(&self) -> bool {
        self.body.has_names()
    }

    /// The `n`th record's name, `None` when the index is unnamed or that
    /// record carries no name.
    ///
    /// Borrowed from an index read whole and owned from one left in its file.
    pub fn name(&self, n: usize) -> Result<Option<Cow<'_, [u8]>>> {
        self.body.name(n)
    }

    /// The position of the first record called `name`.
    pub fn position(&self, name: &[u8]) -> Result<Option<usize>> {
        // by_name is ordered by name bytes, so the first
        // entry not ordered before `name` is the only one
        // that can equal it. every entry names a record that
        // carries one -- name_with builds it that way and
        // read_from refuses anything else -- so the fallback
        // is unreachable rather than a default
        match &self.body {
            Body::Memory { names, .. } => {
                let Some(names) = names else {
                    return Ok(None);
                };

                let at = names
                    .by_name
                    .partition_point(|&n| names.name(n as usize).unwrap_or_default() < name);

                let Some(&n) = names.by_name.get(at) else {
                    return Ok(None);
                };
                let n = n as usize;

                Ok((names.name(n) == Some(name)).then_some(n))
            }

            Body::File { layout, .. } => {
                let Some(names) = &layout.names else {
                    return Ok(None);
                };

                // the same partition point, one read a probe
                // rather than an indexed lookup
                let (mut lo, mut hi) = (0usize, names.by_name_len);
                while lo < hi {
                    let mid = lo + (hi - lo) / 2;
                    let at = self.body.ordered(mid)?;

                    if self.body.name(at)?.unwrap_or_default().as_ref() < name {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }

                if lo == names.by_name_len {
                    return Ok(None);
                }

                let n = self.body.ordered(lo)?;

                Ok((self.body.name(n)?.as_deref() == Some(name)).then_some(n))
            }
        }
    }

    /// Name every record from the bytes its [`Offset`] frames.
    ///
    /// # Errors
    ///
    /// Whatever `name_of` reports, and a name past the length an index holds.
    pub fn name_with<F>(&mut self, mut name_of: F) -> Result<()>
    where
        F: for<'a> FnMut(&'a Offset) -> Result<Option<Cow<'a, [u8]>>>,
    {
        let offsets = self.nameable()?;
        let mut names = NameBuilder::with_capacity(offsets.len());

        for offset in offsets {
            names.push(name_of(offset)?.as_deref())?;
        }

        self.attach(names.finish());

        Ok(())
    }

    /// The offsets of an index that can take a name layer.
    fn nameable(&self) -> Result<&[Offset]> {
        // naming walks every record, and the layer it
        // builds has nowhere to go in a file already
        // written
        let Body::Memory { offsets, .. } = &self.body else {
            bail_parse!("an index opened in place cannot be named, only one read whole")
        };

        if u32::try_from(offsets.len()).is_err() {
            bail_parse!(
                "an index of {} records is past the u32 a name lookup addresses",
                offsets.len()
            )
        }

        Ok(offsets)
    }

    /// Attach a name layer, once [`nameable`](Self::nameable) has passed.
    fn attach(&mut self, built: Names) {
        let Body::Memory { names, .. } = &mut self.body else {
            unreachable!("nameable refuses an index opened in place")
        };

        *names = Some(built);
    }

    /// Index every record in `bytes`.
    pub fn build<R: Read>(bytes: R, format: Format) -> Result<Index> {
        let mut reader = Reader::new(bytes, format);
        let mut offsets = Vec::new();

        while reader.advance()? {
            offsets.push(reader.offset());
        }

        Ok(Index::new(format, offsets))
    }

    /// Index every record in `bytes` and name it, in one pass.
    pub fn build_named<R: Read>(bytes: R, format: Format) -> Result<Index> {
        let mut reader = Reader::new(bytes, format);
        let mut offsets = Vec::new();
        let mut names = NameBuilder::default();

        while reader.advance()? {
            offsets.push(reader.offset());

            // read out of the reader's buffer while the
            // record is still in it: naming afterwards
            // through name_with re-reads every record this
            // pass already walked, one pread at a time
            names.push(crate::seq::name_of(format, reader.record()))?;
        }

        let mut index = Index::new(format, offsets);
        index.nameable()?;
        index.attach(names.finish());

        Ok(index)
    }

    /// Read the whole index at `path` into memory.
    pub fn read(path: &Path) -> Result<Index> {
        let (index, _, _) = Index::read_from(BufReader::new(File::open(path)?))?;

        Ok(index)
    }

    /// Attach to the index at `path`, leaving its offsets and names in the
    /// file and reading one when it is asked for.
    ///
    /// Reads the header and nothing else, so opening an index of a billion
    /// records costs the same as opening one of ten.
    pub fn open(path: &Path) -> Result<Index> {
        let file = File::open(path)?;
        let (format, count, has_names) = read_head(&file)?;

        let names = if has_names {
            // the blob's length is the one thing the header
            // does not carry, and every run after it is
            // placed relative to where it ends
            let mut buf = [0u8; 8];
            let blob_len_at = HEADER + count as u64 * OFFSET;
            crate::pread::read_exact_at(&file, &mut buf, blob_len_at)?;
            let blob_len = u64::from_le_bytes(buf);

            let blob_at = blob_len_at + 8;
            let spans_at = blob_at + blob_len;
            let by_name_len_at = spans_at + count as u64 * SPAN;

            crate::pread::read_exact_at(&file, &mut buf, by_name_len_at)?;
            let by_name_len = narrow(u64::from_le_bytes(buf))?;

            Some(NameLayout {
                blob_at,
                spans_at,
                by_name_at: by_name_len_at + 8,
                by_name_len,
            })
        } else {
            None
        };

        Ok(Index {
            format,
            body: Body::File {
                file,
                layout: Layout { count, names },
            },
        })
    }

    /// Write this index to `path`, stamped with `data`'s length and
    /// modification time so [`is_current`] can tell it apart from a file that
    /// has since changed.
    pub fn write(&self, path: &Path, data: &Path) -> Result<()> {
        let metadata = std::fs::metadata(data)?;

        let mut temporary = OsString::from(path.as_os_str());
        temporary.push(".tmp");
        let temporary = PathBuf::from(temporary);

        // written aside and renamed, because rename is atomic
        // where a write is not: a reader never opens a half
        // written index
        let mut out = BufWriter::new(File::create(&temporary)?);
        self.write_to(&mut out, metadata.len(), metadata.modified()?)?;
        out.flush()?;
        drop(out);

        std::fs::rename(&temporary, path)?;

        Ok(())
    }

    fn write_to<W: Write>(&self, mut w: W, len: u64, mtime: SystemTime) -> Result<()> {
        let (secs, nanos) = epoch_parts(mtime);

        // an index left in its file is written back by
        // copying that file, not by this
        let Body::Memory { offsets, names } = &self.body else {
            bail_parse!("an index opened in place is already written")
        };

        w.write_all(&MAGIC)?;
        w.write_all(&[VERSION, format_byte(self.format)])?;
        w.write_all(&len.to_le_bytes())?;
        w.write_all(&secs.to_le_bytes())?;
        w.write_all(&nanos.to_le_bytes())?;
        w.write_all(&(offsets.len() as u64).to_le_bytes())?;
        w.write_all(&[u8::from(names.is_some())])?;

        for offset in offsets {
            w.write_all(&offset.start.to_le_bytes())?;
            w.write_all(&(offset.n_bytes as u64).to_le_bytes())?;
        }

        let Some(names) = names else {
            return Ok(());
        };

        w.write_all(&(names.blob.len() as u64).to_le_bytes())?;
        w.write_all(&names.blob)?;

        for span in &names.spans {
            let (present, start, len) = match span {
                Some(s) => (1u8, s.start(), s.len()),
                None => (0, 0, 0),
            };

            w.write_all(&[present])?;
            w.write_all(&start.to_le_bytes())?;
            w.write_all(&len.to_le_bytes())?;
        }

        w.write_all(&(names.by_name.len() as u64).to_le_bytes())?;
        for &n in &names.by_name {
            w.write_all(&n.to_le_bytes())?;
        }

        Ok(())
    }

    /// The index in `r`, with the length and modification time of the file it
    /// was written for.
    fn read_from<R: Read>(mut r: R) -> Result<(Index, u64, SystemTime)> {
        let magic: [u8; 4] = take(&mut r)?;
        if magic != MAGIC {
            bail_parse!("not a sail index: it opens {magic:?}, not {MAGIC:?}")
        }

        let [version, format] = take(&mut r)?;
        if version != VERSION {
            bail_parse!("a sail index written at version {version}, and this reads {VERSION}")
        }
        let format = format_from_byte(format)?;

        let len = u64::from_le_bytes(take(&mut r)?);
        let secs = u64::from_le_bytes(take(&mut r)?);
        let nanos = u32::from_le_bytes(take(&mut r)?);
        let count = u64::from_le_bytes(take(&mut r)?);
        let [has_names] = take(&mut r)?;

        // not reserved from count: it has been read but not
        // yet corroborated by anything, and a corrupt one
        // would allocate against a file that cannot fill it
        let mut offsets = Vec::new();
        for _ in 0..count {
            offsets.push(read_offset(&mut r)?);
        }

        let names = match has_names {
            0 => None,
            1 => Some(read_names(&mut r, offsets.len())?),
            other => bail_parse!("an index carries names or it does not, and this says {other}"),
        };

        Ok((
            Index {
                format,
                body: Body::Memory { offsets, names },
            },
            len,
            UNIX_EPOCH + std::time::Duration::new(secs, nanos),
        ))
    }
}

// ---

/// Where an index for `data` goes by default: `data`, plus `.saidx`.
pub fn path_for(data: impl AsRef<Path>) -> PathBuf {
    // appended rather than substituted, so proteins.fa
    // and proteins.sto do not both index to
    // proteins.saidx
    let mut name = OsString::from(data.as_ref().as_os_str());
    name.push(".");
    name.push(EXTENSION);

    PathBuf::from(name)
}

/// Whether the index at `path` still describes `data` as it is now.
///
/// `false` for an index that is missing, unreadable, or stamped against a
/// different state of `data`: a caller's answer to all three is to build one.
pub fn is_current(path: &Path, data: &Path) -> Result<bool> {
    let Ok(file) = File::open(path) else {
        return Ok(false);
    };

    let mut buf = [0u8; HEADER as usize];
    if crate::pread::read_exact_at(&file, &mut buf, 0).is_err() {
        return Ok(false);
    }

    if buf[..4] != MAGIC || buf[4] != VERSION {
        return Ok(false);
    }

    let len = u64::from_le_bytes(buf[6..14].try_into().expect("eight bytes"));
    let secs = u64::from_le_bytes(buf[14..22].try_into().expect("eight bytes"));
    let nanos = u32::from_le_bytes(buf[22..26].try_into().expect("four bytes"));

    // the data file is the caller's, so a failure to stat it
    // is theirs to see rather than an answer of false
    let metadata = std::fs::metadata(data)?;

    Ok(len == metadata.len() && (secs, nanos) == epoch_parts(metadata.modified()?))
}

/// The record count of the index at `path`, from its header alone.
pub fn count_in(path: &Path) -> Result<usize> {
    let (_, count, _) = read_head(&File::open(path)?)?;

    Ok(count)
}

/// An index file's format, record count, and whether it carries names.
fn read_head(file: &File) -> Result<(Format, usize, bool)> {
    let mut buf = [0u8; HEADER as usize];
    crate::pread::read_exact_at(file, &mut buf, 0)?;

    if buf[..4] != MAGIC {
        bail_parse!("not a sail index: it opens {:?}, not {MAGIC:?}", &buf[..4])
    }
    if buf[4] != VERSION {
        bail_parse!(
            "a sail index written at version {}, and this reads {VERSION}",
            buf[4]
        )
    }

    let format = format_from_byte(buf[5])?;
    let count = narrow(u64::from_le_bytes(
        buf[26..34].try_into().expect("eight bytes"),
    ))?;

    let has_names = match buf[34] {
        0 => false,
        1 => true,
        other => bail_parse!("an index carries names or it does not, and this says {other}"),
    };

    Ok((format, count, has_names))
}

fn read_offset<R: Read>(r: &mut R) -> Result<Offset> {
    let start = u64::from_le_bytes(take(r)?);
    let n_bytes = narrow(u64::from_le_bytes(take(r)?))?;

    Ok(Offset { start, n_bytes })
}

/// A length off the disk as a `usize`.
//
// usize on a 32-bit target is narrower than what was
// written, and `as` would truncate a record silently
fn narrow(n: u64) -> Result<usize> {
    usize::try_from(n)
        .map_err(|_| Error::parse(format!("a record of {n} bytes does not fit a usize")))
}

fn read_names<R: Read>(r: &mut R, count: usize) -> Result<Names> {
    let blob_len = u64::from_le_bytes(take(r)?);
    if blob_len > MAX_BLOB {
        bail_parse!("an index of {blob_len} bytes of names exceeds the {MAX_BLOB} one holds")
    }

    // take + read_to_end rather than a reserve: the length
    // is uncorroborated, so the buffer grows with what the
    // file actually holds and the shortfall is the error
    let mut blob = Vec::new();
    (&mut *r).take(blob_len).read_to_end(&mut blob)?;
    if blob.len() as u64 != blob_len {
        bail_parse!(
            "an index of {blob_len} bytes of names holds {} of them",
            blob.len()
        )
    }

    let mut spans = Vec::new();
    for _ in 0..count {
        spans.push(read_span(r, blob.len())?);
    }

    let ordered = u64::from_le_bytes(take(r)?);
    if ordered > count as u64 {
        bail_parse!("an index of {count} records orders {ordered} names")
    }

    let mut by_name = Vec::new();
    for _ in 0..ordered {
        let n = u32::from_le_bytes(take(r)?);

        // position reads a name per entry, so every entry
        // has to name a record that carries one
        match spans.get(n as usize) {
            Some(Some(_)) => by_name.push(n),
            _ => bail_parse!("an index orders record {n}, which carries no name"),
        }
    }

    Ok(Names {
        blob,
        spans,
        by_name,
    })
}

fn read_span<R: Read>(r: &mut R, blob_len: usize) -> Result<Option<Span>> {
    let [present] = take(r)?;
    let start = u64::from_le_bytes(take(r)?);
    let len = u32::from_le_bytes(take(r)?);

    match present {
        0 => Ok(None),
        1 => {
            if len as usize > MAX_NAME {
                bail_parse!("a record name of {len} bytes exceeds the {MAX_NAME} an index holds")
            }
            if start + len as u64 > blob_len as u64 {
                bail_parse!("a name at {start} of {len} bytes runs past the {blob_len} stored")
            }

            Ok(Some(Span::new(start, len)))
        }
        other => bail_parse!("a name is present or it is not, and this record says {other}"),
    }
}

fn take<R: Read, const N: usize>(r: &mut R) -> Result<[u8; N]> {
    let mut buf = [0; N];
    r.read_exact(&mut buf)?;

    Ok(buf)
}

/// Seconds and nanoseconds of `t` since the epoch.
fn epoch_parts(t: SystemTime) -> (u64, u32) {
    match t.duration_since(UNIX_EPOCH) {
        Ok(d) => (d.as_secs(), d.subsec_nanos()),

        // a file dated before 1970 encodes as the epoch. it
        // still compares equal to itself, because load_from
        // compares what would be written rather than the
        // SystemTime it was built from
        Err(_) => (0, 0),
    }
}

fn format_byte(format: Format) -> u8 {
    match format {
        Format::Fasta => 0,
        Format::Stockholm => 1,
        Format::Hmm => 2,
    }
}

fn format_from_byte(byte: u8) -> Result<Format> {
    match byte {
        0 => Ok(Format::Fasta),
        1 => Ok(Format::Stockholm),
        2 => Ok(Format::Hmm),
        other => bail_parse!("a sail index of no format this knows: {other}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixture;

    /// A one-byte initiating delimiter, standing in for FASTA's `>`.
    struct InitiatingDelim;

    impl Delimit for InitiatingDelim {
        const DELIM: &'static [u8] = b">";
        const DELIM_TYPE: Delimiter = Delimiter::Initiating;
    }

    /// A multi-byte terminating delimiter, standing in for Stockholm's `//`.
    struct TerminatingDelim;

    impl Delimit for TerminatingDelim {
        const DELIM: &'static [u8] = b"//";
        const DELIM_TYPE: Delimiter = Delimiter::Terminating;
    }

    /// Two bytes that differ, unlike `//`, where the byte a search leads with
    /// is also the byte that continues the match.
    struct PairDelim;

    impl Delimit for PairDelim {
        const DELIM: &'static [u8] = b"@>";
        const DELIM_TYPE: Delimiter = Delimiter::Initiating;
    }

    /// Four bytes of anchored pattern, so it can straddle a chunk boundary in
    /// three different ways rather than two.
    struct WideDelim;

    impl Delimit for WideDelim {
        const DELIM: &'static [u8] = b">>>";
        const DELIM_TYPE: Delimiter = Delimiter::Terminating;
    }

    fn positions_at<D: Delimit>(input: &[u8], chunk: usize) -> Vec<u64> {
        scan::<_, D>(input, chunk).unwrap().0
    }

    fn off(start: u64, n_bytes: usize) -> Offset {
        Offset { start, n_bytes }
    }

    /// Every anchored delimiter in `input`, by the definition: [`Delimit::DELIM`]
    /// at a line start, where byte 0 counts as one.
    fn positions_naive<D: Delimit>(input: &[u8]) -> Vec<u64> {
        // the specification scan is an optimization of: no
        // buffer, no carry and no chunk, so it cannot go
        // wrong at a boundary
        (0..input.len())
            .filter(|&d| (d == 0 || input[d - 1] == b'\n') && input[d..].starts_with(D::DELIM))
            .map(|d| d as u64)
            .collect()
    }

    /// A xorshift64 stream, so a failure reproduces forever.
    fn soup(seed: u64, len: usize) -> Vec<u8> {
        // why: a uniform draw over 256 values puts a '>'
        //      after a '\n' once per 65,536 bytes, so a
        //      short input would contain no boundary at all
        const ALPHABET: &[u8] = b"\n\n>//ACGT";

        // xorshift64 by hand rather than seeded from rand,
        // which does not promise StdRng is byte-identical
        // across versions
        let mut x = seed | 1;
        (0..len)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                ALPHABET[(x >> 33) as usize % ALPHABET.len()]
            })
            .collect()
    }

    /// Sweep short dense inputs at every small chunk size, against the
    /// definition.
    fn agrees_with_the_definition<D: Delimit>(delim: &str) {
        for seed in 1..=120u64 {
            let input = soup(seed, seed as usize * 3 % 200);

            for chunk in 1..=9 {
                assert_eq!(
                    positions_at::<D>(&input, chunk),
                    positions_naive::<D>(&input),
                    "{delim} disagreed at chunk {chunk}, seed {seed}, \
                     input {:?}",
                    String::from_utf8_lossy(&input)
                );
            }
        }
    }

    #[test]
    fn the_scan_finds_exactly_what_the_definition_says_it_should() {
        // short dense inputs at every small chunk size: what
        // buffering risks is a delimiter lost or doubled at a
        // boundary, and a chunk of 2 against a four-byte
        // pattern splits it every way it can be split
        agrees_with_the_definition::<InitiatingDelim>(">");
        agrees_with_the_definition::<TerminatingDelim>("//");
        agrees_with_the_definition::<PairDelim>("@>");
        agrees_with_the_definition::<WideDelim>(">>>");
    }

    /// A [`Read`] that hands over a few bytes at a time, and sometimes an
    /// error instead.
    struct Choppy<'a> {
        rest: &'a [u8],
        state: u64,
    }

    impl Read for Choppy<'_> {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if self.rest.is_empty() || buf.is_empty() {
                return Ok(0);
            }

            self.state = self
                .state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);

            // cannot starve the scan: the state advances on
            // every call, so a run of interrupts is bounded
            //
            // nothing else in these tests reaches read_some's
            // Interrupted arm, where an error taken for EOF
            // would index a prefix in silence
            if self.state >> 61 == 0 {
                return Err(std::io::ErrorKind::Interrupted.into());
            }

            // read size and chunk size are different axes:
            // scan picks the chunk, the source picks the
            // read, and this is deliberately never the whole
            // buffer
            let n = ((self.state >> 58) as usize % 7 + 1)
                .min(buf.len())
                .min(self.rest.len());
            buf[..n].copy_from_slice(&self.rest[..n]);
            self.rest = &self.rest[n..];

            Ok(n)
        }
    }

    #[test]
    fn an_interrupted_read_is_retried_rather_than_reported() {
        // note: read_some's filling loop looks load-bearing
        //       and is not. scan loops and carries bytes
        //       itself, so a short read costs a round and
        //       changes no answer -- the loop keeps rounds
        //       large, not correct. only the Interrupted arm
        //       is pinned here
        for seed in 1..=40u64 {
            let input = soup(seed, 700 + seed as usize);
            let choppy = Choppy {
                rest: &input,
                state: seed,
            };

            assert_eq!(
                scan::<_, TerminatingDelim>(choppy, 16).unwrap().0,
                positions_naive::<TerminatingDelim>(&input),
                "an interrupted read lost or invented a delimiter, seed {seed}"
            );
        }
    }

    #[test]
    fn the_scan_agrees_with_the_definition_on_an_input_longer_than_a_chunk() {
        // the small-input test never reaches CHUNK, so it
        // never exercises the buffer the shipped code uses
        let input = soup(0xB10, CHUNK * 2 + 517);

        assert_eq!(
            build_index::<_, InitiatingDelim>(&input[..]).unwrap().len(),
            positions_naive::<InitiatingDelim>(&input).len(),
            "one offset per delimiter"
        );
        for delim_positions in [
            (
                scan::<_, InitiatingDelim>(&input[..], CHUNK).unwrap().0,
                positions_naive::<InitiatingDelim>(&input),
            ),
            (
                scan::<_, TerminatingDelim>(&input[..], CHUNK).unwrap().0,
                positions_naive::<TerminatingDelim>(&input),
            ),
        ] {
            let (found, want) = delim_positions;
            assert!(want.len() > 1000, "the input has to be dense to be a test");
            assert_eq!(found, want);
        }
    }

    #[test]
    fn the_carry_is_one_byte_short_of_the_anchored_pattern() {
        // both derived, so an impl cannot set one without the other
        assert_eq!(InitiatingDelim::WIDTH, 2);
        assert_eq!(InitiatingDelim::CARRY, 1);
        assert_eq!(TerminatingDelim::WIDTH, 3);
        assert_eq!(TerminatingDelim::CARRY, 2);
    }

    #[test]
    fn scanning_finds_every_delimiter_in_order() {
        let input = b">a\nSEQ\n>b\nSEQSEQ\n>c\nQ\n";
        let (positions, total) = scan::<_, InitiatingDelim>(&input[..], 64).unwrap();
        assert_eq!(positions, vec![0, 7, 17]);
        assert_eq!(total, input.len() as u64);
    }

    #[test]
    fn the_chunk_size_never_changes_what_is_found() {
        // the point of the carry: a delimiter split across a
        // buffer boundary must be found anyway, so every
        // chunk size has to agree. a chunk of 1 puts a
        // boundary between every pair of bytes
        let input = b">a\nSEQ\n>b\nSEQSEQ\n>c\nQ\n";
        for chunk in 1..=input.len() + 1 {
            assert_eq!(
                positions_at::<InitiatingDelim>(input, chunk),
                vec![0, 7, 17],
                "chunk size {chunk} disagreed"
            );
        }
    }

    #[test]
    fn a_multibyte_delimiter_survives_every_chunk_boundary() {
        // three bytes of pattern means it can straddle a
        // boundary in two different ways
        let input = b"x AC\n//\ny GT\n//\n";
        let want = positions_at::<TerminatingDelim>(input, input.len() + 1);
        assert_eq!(want, vec![5, 13], "the two terminators");
        for chunk in 1..=input.len() + 1 {
            assert_eq!(
                positions_at::<TerminatingDelim>(input, chunk),
                want,
                "chunk size {chunk} disagreed"
            );
        }
    }

    #[test]
    fn a_delimiter_mid_line_is_not_a_record_boundary() {
        // why the anchor lives in the trait: Pfam writes
        // wikipedia links into #=GF WK lines, and the // in
        // that URL does not end the family
        let input = b"# STOCKHOLM 1.0\n#=GF WK   http://en.wikipedia.org/wiki/X\nx AC\n//\n";
        let (positions, total) = scan::<_, TerminatingDelim>(&input[..], 64).unwrap();
        assert_eq!(positions.len(), 1, "only the real terminator");
        assert_eq!(&input[positions[0] as usize..][..2], b"//");
        assert_eq!(frame::<TerminatingDelim>(&positions, total).len(), 1);
    }

    #[test]
    fn a_delimiter_at_byte_zero_needs_no_newline_in_front_of_it() {
        // what the seeded newline is for: the first record has
        // no line ending before it, and still matches through
        // the ordinary path
        assert_eq!(positions_at::<InitiatingDelim>(b">a\nAC\n", 64), vec![0]);

        // a lone delimiter is a whole record
        assert_eq!(positions_at::<InitiatingDelim>(b">", 1), vec![0]);
    }

    #[test]
    fn input_with_nothing_to_find_yields_nothing() {
        for input in [&b""[..], &b"ACGT\nACGT\n"[..], &b"\n"[..], &b"/"[..]] {
            assert!(positions_at::<InitiatingDelim>(input, 4).is_empty());
            assert!(positions_at::<TerminatingDelim>(input, 4).is_empty());
        }
        assert_eq!(scan::<_, InitiatingDelim>(&b""[..], 4).unwrap().1, 0);
    }

    #[test]
    fn an_initiating_delimiter_opens_the_record_it_starts() {
        assert_eq!(
            frame::<InitiatingDelim>(&[0, 7, 17], 22),
            vec![off(0, 7), off(7, 10), off(17, 5)]
        );
    }

    #[test]
    fn a_terminating_delimiter_closes_the_record_it_ends() {
        // "A\n//\nB\n//\n" -- each // line belongs to the
        // record it terminates, trailing newline included
        assert_eq!(
            frame::<TerminatingDelim>(&[2, 7], 10),
            vec![off(0, 5), off(5, 5)]
        );
    }

    #[test]
    fn a_final_terminator_without_a_trailing_newline_still_frames() {
        // "A\n//" -- no newline to include, so the record
        // clamps to the end rather than running past it
        assert_eq!(frame::<TerminatingDelim>(&[2], 4), vec![off(0, 4)]);
    }

    #[test]
    fn bytes_outside_any_record_are_not_indexed() {
        // junk ahead of the first initiating delimiter
        assert_eq!(frame::<InitiatingDelim>(&[4], 10), vec![off(4, 6)]);

        // a truncated record after the last terminator
        assert_eq!(frame::<TerminatingDelim>(&[2], 12), vec![off(0, 5)]);
    }

    #[test]
    fn no_delimiters_means_no_offsets_in_either_direction() {
        assert!(frame::<InitiatingDelim>(&[], 10).is_empty());
        assert!(frame::<TerminatingDelim>(&[], 10).is_empty());
    }

    #[test]
    fn the_offsets_tile_the_input_with_no_gap_or_overlap() {
        let input = b">a\nSEQ\n>b\nSEQSEQ\n>c\nQ\n";
        let offsets = build_index::<_, InitiatingDelim>(&input[..]).unwrap();
        assert_eq!(offsets.len(), 3);

        let mut next = 0u64;
        for o in &offsets {
            assert_eq!(o.start, next, "records must be contiguous");
            assert_eq!(
                &input[o.start as usize..][..1],
                b">",
                "starts at its delimiter"
            );
            next = o.start + o.n_bytes as u64;
        }
        assert_eq!(
            next as usize,
            input.len(),
            "the last record runs to the end"
        );
    }

    // ---

    /// The length and modification time [`sample`] is stamped with.
    const SOURCE_LEN: u64 = 4869;
    const SOURCE_MTIME: (u64, u32) = (1_757_000_000, 123_456_789);

    fn stamp() -> SystemTime {
        UNIX_EPOCH + std::time::Duration::new(SOURCE_MTIME.0, SOURCE_MTIME.1)
    }

    /// An index over a named record, an unnamed one, and one named `""`.
    fn sample() -> Index {
        let names = [Some(&b"SH3_1"[..]), None, Some(&b""[..])];
        let mut index = Index::new(
            Format::Stockholm,
            vec![off(0, 2400), off(2400, 2469), off(4869, 0)],
        );

        let mut n = 0;
        index
            .name_with(|_| {
                let name = names[n];
                n += 1;

                Ok(name.map(Cow::Borrowed))
            })
            .unwrap();

        index
    }

    fn written(index: &Index) -> Vec<u8> {
        let mut bytes = Vec::new();
        index.write_to(&mut bytes, SOURCE_LEN, stamp()).unwrap();

        bytes
    }

    fn read_back(bytes: &[u8]) -> Result<Index> {
        Index::read_from(bytes).map(|(index, _, _)| index)
    }

    /// Copy a fixture somewhere writable, since saving an index writes beside
    /// the file it describes.
    fn scratch(tag: &str, from: &str) -> PathBuf {
        let name = Path::new(from)
            .file_name()
            .unwrap()
            .to_string_lossy()
            .into_owned();
        let to = std::env::temp_dir().join(format!("libsail-{tag}-{}-{name}", std::process::id()));
        std::fs::copy(from, &to).unwrap();

        to
    }

    /// Everything two indexes are expected to agree on, since an index that
    /// may hold an open file cannot derive PartialEq.
    fn same(a: &Index, b: &Index) {
        assert_eq!(a.format(), b.format());
        assert_eq!(a.len(), b.len());
        assert_eq!(a.has_names(), b.has_names());

        for n in 0..a.len() {
            assert_eq!(a.offset(n).unwrap(), b.offset(n).unwrap(), "offset {n}");
            assert_eq!(a.name(n).unwrap(), b.name(n).unwrap(), "name {n}");
        }
    }

    fn clean(path: &Path) {
        std::fs::remove_file(path_for(path)).ok();
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn an_index_opened_in_place_answers_what_one_read_whole_does() {
        // the two backings of one file, asked everything the
        // type can be asked. an opened index reaches its
        // offsets, spans and ordering a read at a time, so a
        // disagreement here is arithmetic on the layout
        let path = scratch("opened", fixture!("families.sto"));
        let at = path_for(&path);

        let index = Index::build_named(File::open(&path).unwrap(), Format::Stockholm).unwrap();
        index.write(&at, &path).unwrap();

        let whole = Index::read(&at).unwrap();
        let opened = Index::open(&at).unwrap();

        same(&whole, &opened);
        assert!(whole.has_names(), "the fixture names its alignments");

        for n in 0..whole.len() {
            let Some(name) = whole.name(n).unwrap() else {
                continue;
            };

            assert_eq!(
                opened.position(&name).unwrap(),
                whole.position(&name).unwrap(),
                "looking up record {n} by name"
            );
        }

        assert_eq!(opened.position(b"not-a-family").unwrap(), None);
        assert_eq!(opened.offset(whole.len()).unwrap(), None, "past the end");

        clean(&path);
    }

    #[test]
    fn only_an_index_read_whole_lends_its_offsets() {
        let path = scratch("lends", fixture!("proteins.fa"));
        let at = path_for(&path);

        Index::build(File::open(&path).unwrap(), Format::Fasta)
            .unwrap()
            .write(&at, &path)
            .unwrap();

        assert!(Index::read(&at).unwrap().offsets().is_some());
        assert!(
            Index::open(&at).unwrap().offsets().is_none(),
            "an opened index has no slice to lend"
        );

        clean(&path);
    }

    #[test]
    fn an_index_opened_in_place_refuses_to_be_named() {
        let path = scratch("naming", fixture!("proteins.fa"));
        let at = path_for(&path);

        Index::build(File::open(&path).unwrap(), Format::Fasta)
            .unwrap()
            .write(&at, &path)
            .unwrap();

        let err = Index::open(&at)
            .unwrap()
            .name_with(|_| Ok(None))
            .expect_err("naming walks every record and rewrites the layer");
        assert!(err.to_string().contains("read whole"), "{err}");

        clean(&path);
    }

    #[test]
    fn a_count_is_read_from_the_header_alone() {
        let path = scratch("counted", fixture!("families.sto"));
        let at = path_for(&path);

        let index = Index::build(File::open(&path).unwrap(), Format::Stockholm).unwrap();
        index.write(&at, &path).unwrap();

        assert_eq!(count_in(&at).unwrap(), index.len());

        clean(&path);
    }

    #[test]
    fn a_written_index_reads_back_equal() {
        let index = sample();

        same(&read_back(&written(&index)).unwrap(), &index);
    }

    #[test]
    fn an_index_with_no_name_round_trips_the_absence_not_an_empty_string() {
        // a Stockholm alignment with no #=GF ID has no name
        // at all, which is a different fact from being named
        // "". conflating them would make an unnamed alignment
        // answer to the empty string
        let back = read_back(&written(&sample())).unwrap();

        assert_eq!(back.name(0).unwrap().as_deref(), Some(&b"SH3_1"[..]));
        assert_eq!(back.name(1).unwrap(), None);
        assert_eq!(back.name(2).unwrap().as_deref(), Some(&b""[..]));
    }

    #[test]
    fn an_unnamed_record_is_not_found_by_any_name() {
        let index = sample();

        assert_eq!(index.position(b"SH3_1").unwrap(), Some(0));
        assert_eq!(index.position(b"").unwrap(), Some(2));
        assert_eq!(index.position(b"PDZ").unwrap(), None);
    }

    #[test]
    fn position_returns_the_first_of_a_repeated_name() {
        // a repeated identifier is a malformed file, and
        // answering with the last would hand back a record
        // the file does not lead with -- the rule
        // collection::Indexable::index_by follows
        let mut index = Index::new(
            Format::Fasta,
            vec![off(0, 10), off(10, 10), off(20, 10), off(30, 10)],
        );

        let names = [&b"b"[..], b"a", b"a", b"c"];
        let mut n = 0;
        index
            .name_with(|_| {
                let name = names[n];
                n += 1;

                Ok(Some(Cow::Borrowed(name)))
            })
            .unwrap();

        assert_eq!(index.position(b"a").unwrap(), Some(1));
        assert_eq!(index.position(b"b").unwrap(), Some(0));
        assert_eq!(index.position(b"c").unwrap(), Some(3));
    }

    #[test]
    fn many_repeats_in_scrambled_order_each_answer_with_their_first_record() {
        // past the length an unstable sort finishes by
        // insertion, so the tie-break on position is what
        // keeps a repeat in file order
        let names: Vec<Vec<u8>> = (0..1000u32)
            .map(|n| format!("r{}", n * 7919 % 13).into_bytes())
            .collect();
        let mut index = Index::new(
            Format::Fasta,
            (0..names.len()).map(|n| off(n as u64 * 10, 10)).collect(),
        );

        index
            .name_with(|o| Ok(Some(Cow::Owned(names[o.start as usize / 10].clone()))))
            .unwrap();

        for name in &names {
            let first = names.iter().position(|n| n == name);
            assert_eq!(index.position(name).unwrap(), first);
        }
    }

    #[test]
    fn a_span_at_the_limits_unpacks_to_what_was_packed() {
        assert_eq!(std::mem::size_of::<Option<Span>>(), 8);

        let start = MAX_BLOB - MAX_NAME as u64;
        let span = Span::new(start, MAX_NAME as u32);
        assert_eq!((span.start(), span.len()), (start, MAX_NAME as u32));

        let empty = Span::new(0, 0);
        assert_eq!((empty.start(), empty.len()), (0, 0));
    }

    #[test]
    fn an_index_with_no_names_answers_no_lookup_and_still_round_trips() {
        let index = Index::new(Format::Fasta, vec![off(0, 10), off(10, 4)]);

        assert!(!index.has_names());
        assert_eq!(index.name(0).unwrap(), None);
        assert_eq!(index.position(b"a").unwrap(), None);
        same(&read_back(&written(&index)).unwrap(), &index);
    }

    #[test]
    fn the_offsets_come_back_in_file_order() {
        assert_eq!(
            read_back(&written(&sample())).unwrap().offsets(),
            Some(&[off(0, 2400), off(2400, 2469), off(4869, 0)][..])
        );
    }

    #[test]
    fn a_truncated_index_is_an_error_and_never_a_panic() {
        // cut at every length, so the header, the offsets,
        // the name blob, a span and the ordering each get
        // truncated
        let bytes = written(&sample());
        for cut in 0..bytes.len() {
            assert!(
                read_back(&bytes[..cut]).is_err(),
                "{cut} bytes read back as a whole index"
            );
        }
    }

    #[test]
    fn an_unrecognised_magic_is_rejected_before_anything_else_is_read() {
        let mut bytes = written(&sample());
        bytes[0] = b'X';

        assert!(read_back(&bytes).is_err());
    }

    #[test]
    fn a_version_this_does_not_know_is_rejected_rather_than_misread() {
        // what makes bumping VERSION a safe way to change the
        // layout: an older reader stops instead of reading
        // the new bytes as the old shape
        let mut bytes = written(&sample());
        bytes[4] = VERSION + 1;

        assert!(read_back(&bytes).is_err());
    }

    #[test]
    fn a_format_byte_naming_no_format_is_rejected() {
        let mut bytes = written(&sample());
        bytes[5] = 9;

        assert!(read_back(&bytes).is_err());
    }

    #[test]
    fn a_record_count_larger_than_the_offsets_is_an_error_not_an_allocation() {
        // the count is read before anything corroborates it;
        // reserving from it would let a corrupt index ask for
        // memory no file could fill
        let mut bytes = written(&sample());
        bytes[26..34].copy_from_slice(&u64::MAX.to_le_bytes());

        assert!(read_back(&bytes).is_err());
    }

    #[test]
    fn a_name_longer_than_an_index_holds_is_refused() {
        let mut index = Index::new(Format::Fasta, vec![off(0, 1)]);

        assert!(
            index
                .name_with(|_| Ok(Some(Cow::Owned(vec![b'x'; MAX_NAME + 1]))))
                .is_err()
        );
    }

    #[test]
    fn an_index_path_is_appended_so_two_formats_do_not_collide() {
        assert_eq!(
            path_for("data/proteins.fa"),
            PathBuf::from("data/proteins.fa.saidx")
        );
        assert_eq!(
            path_for("data/proteins.sto"),
            PathBuf::from("data/proteins.sto.saidx")
        );
    }

    #[test]
    fn an_index_saved_beside_a_file_addresses_what_a_fresh_scan_of_it_does() {
        let path = scratch("saved", fixture!("families.sto"));
        let offsets =
            build_index::<_, TerminatingDelim>(std::fs::File::open(&path).unwrap()).unwrap();

        let mut index = Index::new(Format::Stockholm, offsets.clone());
        index.name_with(|_| Ok(None)).unwrap();
        index.write(&path_for(&path), &path).unwrap();

        let back = Index::read(&path_for(&path)).unwrap();
        assert_eq!(back.format(), Format::Stockholm);
        assert_eq!(back.offsets(), Some(&offsets[..]));
        assert!(back.has_names());

        clean(&path);
    }

    #[test]
    fn an_index_saved_beside_a_file_is_refused_once_the_file_changes() {
        // size and mtime, not a digest: verifying a digest
        // means re-reading the file, which is the read the
        // index exists to avoid. a same-size edit inside one
        // mtime tick is the gap that leaves, and it is not
        // closed
        let path = scratch("stale", fixture!("families.sto"));
        let offsets =
            build_index::<_, TerminatingDelim>(std::fs::File::open(&path).unwrap()).unwrap();
        Index::new(Format::Stockholm, offsets)
            .write(&path_for(&path), &path)
            .unwrap();

        assert!(is_current(&path_for(&path), &path).unwrap());

        let mut bytes = std::fs::read(&path).unwrap();
        bytes.extend_from_slice(b"# STOCKHOLM 1.0\n//\n");
        std::fs::write(&path, &bytes).unwrap();

        assert!(!is_current(&path_for(&path), &path).unwrap());

        clean(&path);
    }
}

// ---

/// The delimiter data for a format: the same two constants [`Delimit`]
/// carries, reached through a value rather than a type parameter.
fn delim_of(format: Format) -> (&'static [u8], Delimiter) {
    use crate::seq::fasta::FastaDelimiter;
    use crate::seq::p7hmm::HmmDelimiter;
    use crate::seq::stockholm::StockholmDelimiter;

    match format {
        Format::Fasta => (FastaDelimiter::DELIM, FastaDelimiter::DELIM_TYPE),
        Format::Stockholm => (StockholmDelimiter::DELIM, StockholmDelimiter::DELIM_TYPE),
        Format::Hmm => (HmmDelimiter::DELIM, HmmDelimiter::DELIM_TYPE),
    }
}

/// Index every record in `bytes` as `format`.
fn build_for<R: Read>(bytes: R, format: Format) -> Result<Vec<Offset>> {
    use crate::seq::fasta::FastaDelimiter;
    use crate::seq::p7hmm::HmmDelimiter;
    use crate::seq::stockholm::StockholmDelimiter;

    match format {
        Format::Fasta => build_index::<_, FastaDelimiter>(bytes),
        Format::Stockholm => build_index::<_, StockholmDelimiter>(bytes),
        Format::Hmm => build_index::<_, HmmDelimiter>(bytes),
    }
}

/// How many records `bytes` holds, by one pass over them.
pub fn count<R: Read>(bytes: R, format: Format) -> Result<usize> {
    let mut reader = Reader::new(bytes, format);
    while reader.advance()? {}

    Ok(reader.framed)
}

/// One record's bytes at a time, framed out of a stream.
///
/// Unlike [`build_index`], nothing is kept per record: a pass over a file
/// costs the buffer and no more.
///
/// ```no_run
/// # use libsail::{format::Format, index::Reader};
/// # fn main() -> libsail::Result<()> {
/// let mut reader = Reader::new(std::fs::File::open("proteins.fa")?, Format::Fasta);
/// while reader.advance()? {
///     let record = reader.record();
///     // ...
/// }
/// # Ok(())
/// # }
/// ```
pub struct Reader<R: Read> {
    source: R,
    format: Format,
    delim: &'static [u8],
    kind: Delimiter,

    /// `'\n' + DELIM`, the anchored pattern [`scan`] also searches for.
    needle: Vec<u8>,

    /// A newline that is not in the file, then the file's bytes: `buf[i]` is at
    /// file offset `base + i - 1`.
    //
    // the same seeded sentinel scan uses, and for the same
    // reason: it stands in for the line ending before byte
    // 0, so a file opening with a delimiter needs no first
    // round special case
    //
    // **note: the sentinel is per stream, not per file. a
    //         Read::chain of two files loses the second
    //         file's first record when the first file has no
    //         final newline. Readers frames each file apart
    buf: Vec<u8>,
    filled: usize,
    base: u64,

    /// Where the next record begins in `buf`.
    cursor: usize,

    /// The current record, as a range of `buf`.
    start: usize,
    len: usize,

    /// The index of the current record.
    position: Option<usize>,

    /// How many records [`advance`](Self::advance) has framed.
    framed: usize,

    eof: bool,
    started: bool,
    done: bool,
}

impl<R: Read> Reader<R> {
    /// Frame `bytes` as `format`.
    pub fn new(bytes: R, format: Format) -> Self {
        let (delim, kind) = delim_of(format);

        let mut needle = Vec::with_capacity(1 + delim.len());
        needle.push(b'\n');
        needle.extend_from_slice(delim);

        let mut buf = vec![0u8; CHUNK + delim.len()];
        buf[0] = b'\n';

        Reader {
            source: bytes,
            format,
            delim,
            kind,
            needle,
            buf,
            filled: 1,
            base: 0,
            cursor: 1,
            start: 1,
            len: 0,
            position: None,
            framed: 0,
            eof: false,
            started: false,
            done: false,
        }
    }

    /// Frame `bytes` as `format`, with a buffer of `chunk` fresh bytes.
    #[cfg(test)]
    fn with_chunk(bytes: R, format: Format, chunk: usize) -> Self {
        let mut reader = Reader::new(bytes, format);

        // chunk changes how often the buffer refills and
        // nothing else, which is what lets the framing tests
        // vary it across every delimiter boundary
        reader.buf = vec![0u8; chunk.max(1) + reader.delim.len()];
        reader.buf[0] = b'\n';

        reader
    }

    /// Frame `bytes` as whatever format they open with.
    pub fn detect(bytes: R) -> Result<Self> {
        let mut reader = Reader::new(bytes, Format::Fasta);
        reader.fill()?;

        let format = crate::format::detect(&reader.buf[1..reader.filled])?;
        let (delim, kind) = delim_of(format);

        reader.needle.truncate(1);
        reader.needle.extend_from_slice(delim);
        reader.format = format;
        reader.delim = delim;
        reader.kind = kind;

        Ok(reader)
    }

    pub fn format(&self) -> Format {
        self.format
    }

    /// The current record's bytes; empty before the first [`advance`](Self::advance).
    /// An advance that returns `false` leaves the last record current.
    pub fn record(&self) -> &[u8] {
        &self.buf[self.start..self.start + self.len]
    }

    /// The index of the current record; `None` before the first
    /// [`advance`](Self::advance). An advance that returns `false` leaves the
    /// last record current, so `position() + 1` is then the record count.
    pub fn position(&self) -> Option<usize> {
        self.position
    }

    /// Where the current record sits in the stream.
    pub fn offset(&self) -> Offset {
        Offset {
            start: self.base + self.start as u64 - 1,
            n_bytes: self.len,
        }
    }

    /// Advance to the next record. `false` at the end of the stream.
    pub fn advance(&mut self) -> Result<bool> {
        let more = !self.done
            && match self.kind {
                Delimiter::Initiating => self.advance_initiating()?,
                Delimiter::Terminating => self.advance_terminating()?,
            };

        if more {
            self.position = Some(self.framed);
            self.framed += 1;
        }

        Ok(more)
    }

    /// Advance until record `n` is current. `false` if the stream ends first,
    /// with the last record left current.
    ///
    /// # Panics
    ///
    /// If record `n` is behind the current one.
    pub fn advance_to(&mut self, n: usize) -> Result<bool> {
        let floor = self.position.unwrap_or(self.framed);
        assert!(n >= floor, "advance_to({n}) is behind record {floor}");

        while self.framed <= n {
            if !self.advance()? {
                return Ok(false);
            }
        }

        Ok(true)
    }

    /// A record runs from its own delimiter to the next one, or to the end for
    /// the last. Bytes before the first delimiter belong to no record.
    fn advance_initiating(&mut self) -> Result<bool> {
        if !self.started {
            self.started = true;

            match self.find(true)? {
                Some(first) => self.cursor = first,
                None => {
                    self.done = true;

                    return Ok(false);
                }
            }
        }

        // the next delimiter ends this record, so a record is
        // only complete once the one after it is in hand --
        // the lookahead the terminating formats do not need
        let end = match self.find(false)? {
            Some(next) => next,
            None => {
                self.done = true;

                self.filled
            }
        };

        self.start = self.cursor;
        self.len = end - self.cursor;
        self.cursor = end;

        Ok(true)
    }

    /// A record runs from the end of the previous one through its delimiter's
    /// own line. Bytes after the last delimiter are a truncated record.
    fn advance_terminating(&mut self) -> Result<bool> {
        self.started = true;

        let Some(at) = self.find(true)? else {
            self.done = true;

            return Ok(false);
        };

        // + 1 for the newline ending the delimiter's line; a
        // final delimiter with no trailing newline clamps to
        // the end, exactly as frame does
        let want = at - self.cursor + self.delim.len() + 1;
        self.ensure(want)?;
        let end = (self.cursor + want).min(self.filled);

        self.start = self.cursor;
        self.len = end - self.cursor;
        self.cursor = end;

        Ok(true)
    }

    /// The buffer index of the next delimiter, refilling until one is found or
    /// the stream ends.
    fn find(&mut self, anchored: bool) -> Result<Option<usize>> {
        // `anchored` includes the byte before the cursor in
        // the search, which is the newline the needle matches
        // on. it is set for a record starting just past the
        // previous one's last byte, and clear for one
        // starting *at* its own delimiter, whose own match
        // would otherwise be found again every round
        let back = usize::from(anchored);

        loop {
            // recomputed every round: a refill compacts, and
            // compaction moves everything by the same amount
            // as the cursor
            let from = self.cursor - back;

            if let Some(j) = memmem::find(&self.buf[from..self.filled], &self.needle) {
                return Ok(Some(from + j + 1));
            }

            if self.eof {
                return Ok(None);
            }

            self.fill()?;
        }
    }

    /// Read until `buf` holds `need` bytes past the cursor, or the stream ends.
    fn ensure(&mut self, need: usize) -> Result<()> {
        while self.filled < self.cursor + need && !self.eof {
            self.fill()?;
        }

        Ok(())
    }

    /// Make room and read once.
    fn fill(&mut self) -> Result<()> {
        if self.eof {
            return Ok(());
        }

        // everything before the byte the cursor anchors on is
        // spent, except the current record, which stays
        // readable until the next one replaces it. the buffer
        // can then hold two records rather than one
        let keep = match self.position {
            Some(_) => self.start,
            None => self.cursor,
        };
        let spent = keep - 1;

        if spent > 0 {
            self.buf.copy_within(spent..self.filled, 0);
            self.filled -= spent;
            self.base += spent as u64;
            self.cursor -= spent;

            if self.position.is_some() {
                self.start -= spent;
            }
        } else if self.filled == self.buf.len() {
            // the record is longer than the buffer, so
            // compacting frees nothing. growing here is what
            // keeps a record spanning k chunks at O(S) rather
            // than the O(S^2/chunk) a fixed buffer would cost
            self.buf.resize(self.buf.len() * 2, 0);
        }

        let n = read_some(&mut self.source, &mut self.buf[self.filled..])?;
        self.filled += n;

        if n == 0 {
            self.eof = true;
        }

        Ok(())
    }
}

// ---

/// Several files' records, framed one file at a time and numbered as one run.
///
/// ```no_run
/// # use libsail::{format::Format, index::Readers};
/// # fn main() -> libsail::Result<()> {
/// let mut readers = Readers::new(["a.fa", "b.fa"], Format::Fasta);
/// while readers.advance()? {
///     let record = readers.record();
///     // ...
/// }
/// # Ok(())
/// # }
/// ```
pub struct Readers {
    paths: Vec<PathBuf>,
    format: Format,

    /// How many records each path holds, once [`count`](Self::count) has run.
    counts: Option<Vec<usize>>,

    /// The reader over `paths[part]`, while it has records left.
    current: Option<Reader<File>>,
    part: usize,

    /// The reader that framed the current record, once its path has run out.
    held: Option<Reader<File>>,

    /// The index of `paths[part]`'s first record.
    base: usize,

    /// The index of the current record.
    position: Option<usize>,
}

impl Readers {
    /// Frame each of `paths` as `format`, in order. Nothing is opened until it
    /// is read.
    pub fn new(paths: impl IntoIterator<Item = impl AsRef<Path>>, format: Format) -> Readers {
        Readers {
            paths: paths
                .into_iter()
                .map(|path| path.as_ref().to_path_buf())
                .collect(),
            format,
            counts: None,
            current: None,
            part: 0,
            held: None,
            base: 0,
            position: None,
        }
    }

    /// How many records the paths hold together, by one pass over each. Once
    /// counted, [`advance_to`](Self::advance_to) skips whole files by these
    /// counts.
    pub fn count(&mut self) -> Result<usize> {
        if self.counts.is_none() {
            let counts = self
                .paths
                .iter()
                .map(|path| count(File::open(path)?, self.format))
                .collect::<Result<Vec<_>>>()?;

            self.counts = Some(counts);
        }

        Ok(self.counts.as_deref().unwrap_or_default().iter().sum())
    }

    /// The current record's bytes; empty before the first [`advance`](Self::advance).
    /// An advance that returns `false` leaves the last record framed current.
    pub fn record(&self) -> &[u8] {
        let holder = match &self.current {
            Some(reader) if reader.position().is_some() => Some(reader),
            _ => self.held.as_ref(),
        };

        holder.map_or(&[], Reader::record)
    }

    /// The index of the current record across every path; `None` before the
    /// first [`advance`](Self::advance). An advance that returns `false` leaves
    /// the last record framed current. That is the last record overall unless
    /// [`advance_to`](Self::advance_to) skipped counted files to get there.
    pub fn position(&self) -> Option<usize> {
        self.position
    }

    /// Advance to the next record. `false` once every path has run out.
    pub fn advance(&mut self) -> Result<bool> {
        loop {
            self.open()?;

            let Some(reader) = self.current.as_mut() else {
                return Ok(false);
            };

            if reader.advance()? {
                self.position = reader.position().map(|at| self.base + at);
                self.held = None;

                return Ok(true);
            }

            self.close();
        }
    }

    /// Advance until record `n` is current. `false` if the paths run out first,
    /// with the last record framed left current.
    ///
    /// # Panics
    ///
    /// If record `n` is behind the current one.
    pub fn advance_to(&mut self, n: usize) -> Result<bool> {
        let floor = self.position.unwrap_or(0);
        assert!(n >= floor, "advance_to({n}) is behind record {floor}");

        if self.position == Some(n) {
            return Ok(true);
        }

        loop {
            self.skip_to(n);
            self.open()?;

            let Some(reader) = self.current.as_mut() else {
                return Ok(false);
            };

            if reader.advance_to(n - self.base)? {
                self.position = Some(n);
                self.held = None;

                return Ok(true);
            }

            self.close();
        }
    }

    /// Open `paths[part]` if no reader is open and a path remains.
    fn open(&mut self) -> Result<()> {
        if self.current.is_none() && self.part < self.paths.len() {
            let file = File::open(&self.paths[self.part])?;
            self.current = Some(Reader::new(file, self.format));
        }

        Ok(())
    }

    /// Move past `paths[part]` once its reader has run out.
    fn close(&mut self) {
        if let Some(reader) = self.current.take() {
            self.base += reader.framed;
            retire(reader, &mut self.held);
        }

        self.part += 1;
    }

    /// Move past every path that ends before record `n`, without opening it,
    /// once the counts are known.
    fn skip_to(&mut self, n: usize) {
        let Some(counts) = &self.counts else {
            return;
        };

        while self.part < counts.len() && self.base + counts[self.part] <= n {
            if let Some(reader) = self.current.take() {
                retire(reader, &mut self.held);
            }

            self.base += counts[self.part];
            self.part += 1;
        }
    }
}

/// Keep `reader` in `held` if it framed any record, since the last record
/// framed is the current one.
fn retire(reader: Reader<File>, held: &mut Option<Reader<File>>) {
    if reader.position().is_some() {
        *held = Some(reader);
    }
}

#[cfg(test)]
mod reader_tests {
    use super::*;
    use crate::fixture;
    use crate::seq::fasta::FastaDelimiter;
    use crate::seq::p7hmm::HmmDelimiter;
    use crate::seq::stockholm::StockholmDelimiter;

    /// Every record the reader frames, as (offset, bytes).
    fn framed(bytes: &[u8], format: Format, chunk: usize) -> Vec<(Offset, Vec<u8>)> {
        let mut reader = Reader::with_chunk(bytes, format, chunk);
        let mut out = Vec::new();

        while reader.advance().unwrap() {
            out.push((reader.offset(), reader.record().to_vec()));
        }

        out
    }

    /// What `build_index` frames, which is the answer the reader has to match.
    fn indexed(bytes: &[u8], format: Format) -> Vec<(Offset, Vec<u8>)> {
        let offsets = match format {
            Format::Fasta => build_index::<_, FastaDelimiter>(bytes),
            Format::Stockholm => build_index::<_, StockholmDelimiter>(bytes),
            Format::Hmm => build_index::<_, HmmDelimiter>(bytes),
        }
        .unwrap();

        offsets
            .into_iter()
            .map(|o| {
                let at = o.start as usize;
                let range = bytes[at..at + o.n_bytes].to_vec();

                (o, range)
            })
            .collect()
    }

    #[test]
    fn the_reader_frames_exactly_what_build_index_frames() {
        // two independent framings of the same bytes, one
        // buffered whole and one streamed: a carry that drops
        // or repeats a delimiter shows up here and nowhere
        // else
        for (path, format) in [
            (fixture!("proteins.fa"), Format::Fasta),
            (fixture!("families.sto"), Format::Stockholm),
            (fixture!("models.hmm"), Format::Hmm),
        ] {
            let bytes = std::fs::read(path).unwrap();
            let want = indexed(&bytes, format);

            assert!(!want.is_empty(), "{path} framed no records");

            // every chunk size from one byte up, so a
            // delimiter straddling a boundary is exercised at
            // every position it can straddle one
            for chunk in 1..=256 {
                assert_eq!(
                    framed(&bytes, format, chunk),
                    want,
                    "{path} at chunk {chunk}"
                );
            }
        }
    }

    #[test]
    fn a_record_longer_than_the_buffer_is_framed_whole() {
        // the case the compact-or-grow rule exists for: with
        // nothing spent to compact away, the buffer has to
        // grow or the record can never be completed
        let mut bytes = Vec::from(&b">big\n"[..]);
        bytes.extend(std::iter::repeat_n(b'A', 100_000));
        bytes.extend_from_slice(b"\n>small\nACGT\n");

        let want = indexed(&bytes, Format::Fasta);
        assert_eq!(want.len(), 2);

        for chunk in [1, 7, 64, 1024] {
            assert_eq!(framed(&bytes, Format::Fasta, chunk), want, "chunk {chunk}");
        }
    }

    #[test]
    fn an_empty_stream_frames_no_records() {
        for format in [Format::Fasta, Format::Stockholm, Format::Hmm] {
            assert!(framed(b"", format, 8).is_empty(), "{format}");
        }
    }

    #[test]
    fn bytes_before_the_first_fasta_record_belong_to_none() {
        // frame's rule for an initiating delimiter, which the
        // reader has to reproduce
        let bytes = b"junk before anything\n>a\nACGT\n";

        assert_eq!(
            framed(bytes, Format::Fasta, 4),
            indexed(bytes, Format::Fasta)
        );
    }

    #[test]
    fn a_truncated_terminating_record_is_not_framed() {
        // trailing bytes with no // after them are an
        // unfinished alignment, and frame drops them
        let bytes = b"# STOCKHOLM 1.0\nseq1 AC\n//\n# STOCKHOLM 1.0\nseq2 GT\n";

        let want = indexed(bytes, Format::Stockholm);
        assert_eq!(want.len(), 1);

        for chunk in 1..=64 {
            assert_eq!(
                framed(bytes, Format::Stockholm, chunk),
                want,
                "chunk {chunk}"
            );
        }
    }

    #[test]
    fn detect_reads_the_format_off_the_stream_and_still_frames_it() {
        for (path, format) in [
            (fixture!("proteins.fa"), Format::Fasta),
            (fixture!("families.sto"), Format::Stockholm),
            (fixture!("models.hmm"), Format::Hmm),
        ] {
            let bytes = std::fs::read(path).unwrap();
            let mut reader = Reader::detect(&bytes[..]).unwrap();

            assert_eq!(reader.format(), format, "{path}");

            let mut n = 0;
            while reader.advance().unwrap() {
                n += 1;
            }
            assert_eq!(n, indexed(&bytes, format).len(), "{path}");
        }
    }

    const FIXTURES: [(&str, Format); 3] = [
        (fixture!("proteins.fa"), Format::Fasta),
        (fixture!("families.sto"), Format::Stockholm),
        (fixture!("models.hmm"), Format::Hmm),
    ];

    #[test]
    fn position_names_the_record_the_reader_holds() {
        for (path, format) in FIXTURES {
            let bytes = std::fs::read(path).unwrap();
            let want = indexed(&bytes, format);

            for chunk in [1, 7, 64, 4096] {
                let mut reader = Reader::with_chunk(&bytes[..], format, chunk);
                assert_eq!(reader.position(), None, "{path} before the first advance");

                let mut n = 0;
                while reader.advance().unwrap() {
                    assert_eq!(reader.position(), Some(n), "{path} at chunk {chunk}");
                    assert_eq!(reader.record(), want[n].1, "{path} record {n}");
                    n += 1;
                }

                assert_eq!(n, want.len(), "{path} at chunk {chunk}");
                assert_eq!(reader.position(), Some(n - 1), "{path} past the end");
                assert_eq!(reader.record(), want[n - 1].1, "{path} past the end");
            }
        }
    }

    #[test]
    fn advance_to_lands_on_the_numbered_record() {
        for (path, format) in FIXTURES {
            let bytes = std::fs::read(path).unwrap();
            let want = indexed(&bytes, format);
            let last = want.len() - 1;

            for chunk in [1, 7, 64, 4096] {
                let mut reader = Reader::with_chunk(&bytes[..], format, chunk);

                // 0 on a fresh reader frames a record rather
                // than taking the empty one as record 0, and
                // the repeat is a no-op on the current record
                for n in [0, 0, 1, last, last] {
                    assert!(reader.advance_to(n).unwrap(), "{path} to {n}");
                    assert_eq!(reader.position(), Some(n), "{path} to {n}");
                    assert_eq!(reader.record(), want[n].1, "{path} to {n}");
                }
            }
        }
    }

    #[test]
    fn advance_to_past_the_end_leaves_the_last_record_current() {
        for (path, format) in FIXTURES {
            let bytes = std::fs::read(path).unwrap();
            let want = indexed(&bytes, format);
            let held = want.len();
            let mut reader = Reader::with_chunk(&bytes[..], format, 7);

            assert!(!reader.advance_to(held).unwrap(), "{path}");
            assert_eq!(reader.position(), Some(held - 1), "{path}");
            assert_eq!(reader.record(), want[held - 1].1, "{path}");
            assert_eq!(reader.offset(), want[held - 1].0, "{path}");

            assert!(!reader.advance_to(held + 10).unwrap(), "{path} once more");
            assert!(reader.advance_to(held - 1).unwrap(), "{path} the last");
        }
    }

    #[test]
    fn an_empty_stream_has_no_current_record_after_the_end() {
        let mut reader = Reader::new(&b""[..], Format::Fasta);

        assert!(!reader.advance().unwrap());
        assert_eq!(reader.position(), None);
        assert!(reader.record().is_empty());
    }

    #[test]
    #[should_panic(expected = "advance_to(1) is behind record 2")]
    fn advance_to_a_record_already_passed_panics() {
        let bytes = std::fs::read(fixture!("proteins.fa")).unwrap();
        let mut reader = Reader::new(&bytes[..], Format::Fasta);

        reader.advance_to(2).unwrap();
        reader.advance_to(1).unwrap();
    }

    #[test]
    fn the_last_record_survives_the_refill_that_finds_the_end() {
        // the bytes after the last delimiter make the final
        // advance refill the buffer, and a refill compacting
        // from the cursor would move the last record's bytes
        let bytes = b"# STOCKHOLM 1.0\na A\n//\ntrailing bytes past the end\n";
        let mut reader = Reader::with_chunk(&bytes[..], Format::Stockholm, 4);

        assert!(reader.advance().unwrap());
        let last = reader.record().to_vec();

        assert!(!reader.advance().unwrap());
        assert_eq!(reader.record(), last);
        assert_eq!(reader.position(), Some(0));
    }

    #[test]
    fn count_agrees_with_the_records_build_index_frames() {
        for (path, format) in FIXTURES {
            let bytes = std::fs::read(path).unwrap();

            assert_eq!(
                count(&bytes[..], format).unwrap(),
                indexed(&bytes, format).len(),
                "{path}"
            );
        }

        assert_eq!(count(&b""[..], Format::Fasta).unwrap(), 0);
    }

    /// `bytes` written to a file of its own under the temp directory.
    fn written(tag: &str, bytes: &[u8]) -> PathBuf {
        let path =
            std::env::temp_dir().join(format!("libsail-readers-{tag}-{}", std::process::id()));
        std::fs::write(&path, bytes).unwrap();

        path
    }

    /// Every record of `paths`, read one file at a time.
    fn each(paths: &[&str], format: Format) -> Vec<Vec<u8>> {
        paths
            .iter()
            .flat_map(|path| indexed(&std::fs::read(path).unwrap(), format))
            .map(|(_, bytes)| bytes)
            .collect()
    }

    const TWO_FILES: [&str; 2] = [fixture!("proteins.fa"), fixture!("proteins_extra.fa")];

    #[test]
    fn readers_number_the_records_of_every_file_as_one_run() {
        let want = each(&TWO_FILES, Format::Fasta);
        let mut readers = Readers::new(TWO_FILES, Format::Fasta);
        assert_eq!(readers.position(), None);

        let mut n = 0;
        while readers.advance().unwrap() {
            assert_eq!(readers.position(), Some(n));
            assert_eq!(readers.record(), want[n], "record {n}");
            n += 1;
        }

        assert_eq!(n, want.len());
        assert_eq!(readers.position(), Some(n - 1));
        assert_eq!(readers.record(), want[n - 1]);
    }

    #[test]
    fn a_file_without_a_final_newline_loses_no_record_of_the_next() {
        // chained, the unterminated last line of the first
        // file swallows the second file's first delimiter.
        // the empty file between them holds no records and
        // must not shift the numbering
        let a = written("unterminated", b">a1\nAC\n>a2\nGT");
        let empty = written("empty", b"");
        let b = written("next", b">b1\nTT\n>b2\nCC\n");

        let chained = std::fs::File::open(&a)
            .unwrap()
            .chain(std::fs::File::open(&b).unwrap());
        assert_eq!(count(chained, Format::Fasta).unwrap(), 3);

        let mut readers = Readers::new([&a, &empty, &b], Format::Fasta);
        let mut got = Vec::new();
        while readers.advance().unwrap() {
            got.push(readers.record().to_vec());
        }

        let want: [&[u8]; 4] = [b">a1\nAC\n", b">a2\nGT", b">b1\nTT\n", b">b2\nCC\n"];
        assert_eq!(got, want);
        assert_eq!(readers.count().unwrap(), 4);

        for path in [a, empty, b] {
            std::fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn count_mid_walk_neither_moves_the_walk_nor_counts_from_it() {
        let want = each(&TWO_FILES, Format::Fasta);
        let mut readers = Readers::new(TWO_FILES, Format::Fasta);

        readers.advance_to(3).unwrap();
        assert_eq!(readers.count().unwrap(), want.len());
        assert_eq!(readers.position(), Some(3));

        assert!(readers.advance().unwrap());
        assert_eq!(readers.record(), want[4]);
    }

    #[test]
    fn advance_to_lands_on_the_same_records_counted_or_not() {
        let want = each(&TWO_FILES, Format::Fasta);
        let first = indexed(&std::fs::read(TWO_FILES[0]).unwrap(), Format::Fasta).len();
        let last = want.len() - 1;

        for counted in [false, true] {
            let mut readers = Readers::new(TWO_FILES, Format::Fasta);
            if counted {
                readers.count().unwrap();
            }

            // either side of the boundary between the files,
            // and a repeat that leaves the record current
            for n in [0, 1, first - 1, first, first, last] {
                assert!(readers.advance_to(n).unwrap(), "to {n}, counted {counted}");
                assert_eq!(readers.position(), Some(n));
                assert_eq!(readers.record(), want[n], "to {n}, counted {counted}");
            }

            assert!(!readers.advance_to(last + 1).unwrap());
            assert_eq!(readers.position(), Some(last));
            assert_eq!(readers.record(), want[last]);
        }
    }

    #[test]
    fn readers_keep_the_last_record_across_empty_files_after_it() {
        let a = written("before-empties", b">a1\nAC\n");
        let empty = written("trailing-empty", b"");

        for counted in [false, true] {
            let mut readers = Readers::new([&a, &empty, &empty], Format::Fasta);
            if counted {
                readers.count().unwrap();
            }

            assert!(readers.advance().unwrap());
            assert!(!readers.advance().unwrap(), "counted {counted}");
            assert_eq!(readers.position(), Some(0));
            assert_eq!(readers.record(), b">a1\nAC\n");

            assert!(!readers.advance_to(5).unwrap(), "counted {counted}");
            assert_eq!(readers.record(), b">a1\nAC\n");
        }

        for path in [a, empty] {
            std::fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn a_counted_file_before_the_target_is_never_opened() {
        // removing the first file after the count proves the
        // walk steps over it: opening it would fail
        let a = written("skipped", b">a1\nAC\n>a2\nGT\n");
        let b = written("landed", b">b1\nTT\n");

        let mut readers = Readers::new([&a, &b], Format::Fasta);
        assert_eq!(readers.count().unwrap(), 3);
        std::fs::remove_file(&a).unwrap();

        assert!(readers.advance_to(2).unwrap());
        assert_eq!(readers.record(), b">b1\nTT\n");

        std::fs::remove_file(b).unwrap();
    }

    #[test]
    fn readers_over_no_paths_hold_no_records() {
        let mut readers = Readers::new(Vec::<PathBuf>::new(), Format::Fasta);

        assert_eq!(readers.count().unwrap(), 0);
        assert!(!readers.advance().unwrap());
        assert!(!readers.advance_to(0).unwrap());
    }

    #[test]
    #[should_panic(expected = "advance_to(1) is behind record 2")]
    fn readers_advance_to_a_record_already_passed_panics() {
        let mut readers = Readers::new(TWO_FILES, Format::Fasta);

        readers.advance_to(2).unwrap();
        readers.advance_to(1).unwrap();
    }
}

#[cfg(test)]
mod build_tests {
    use super::*;
    use crate::fixture;
    use crate::seq::fasta::IndexedFasta;

    #[test]
    fn a_built_index_addresses_the_same_records_as_a_scanned_one() {
        // the one-pass build against the scan-then-name
        // pass, offsets and names both
        for (path, format) in [
            (fixture!("proteins.fa"), Format::Fasta),
            (fixture!("families.sto"), Format::Stockholm),
            (fixture!("models.hmm"), Format::Hmm),
        ] {
            let bytes = std::fs::read(path).unwrap();

            let built = Index::build(&bytes[..], format).unwrap();
            let named = Index::build_named(&bytes[..], format).unwrap();

            assert_eq!(built.len(), named.len(), "{path}");
            assert_eq!(built.offsets(), named.offsets(), "{path}");
            assert!(!built.has_names(), "{path}");
            assert!(named.has_names(), "{path}");

            for n in 0..named.len() {
                let want = crate::seq::name_of(
                    format,
                    &bytes[named.offset(n).unwrap().unwrap().start as usize..]
                        [..named.offset(n).unwrap().unwrap().n_bytes],
                );
                assert_eq!(named.name(n).unwrap().as_deref(), want, "{path} record {n}");
            }
        }
    }

    #[test]
    fn a_built_index_matches_the_one_the_collection_builds() {
        // IndexedFasta::open plus name_index is the other
        // route to a named index, built from preads rather
        // than from the framing buffer
        let path = fixture!("proteins.fa");

        let mut scanned = IndexedFasta::open(path).unwrap();
        scanned.name_index().unwrap();

        let built = Index::build_named(File::open(path).unwrap(), Format::Fasta).unwrap();

        assert_eq!(built.offsets(), scanned.index().offsets());
        for n in 0..built.len() {
            assert_eq!(
                built.name(n).unwrap(),
                scanned.index().name(n).unwrap(),
                "record {n}"
            );
        }
    }

    #[test]
    fn a_built_index_saves_and_loads_like_any_other() {
        let path = std::env::temp_dir().join(format!("libsail-build-{}.fa", std::process::id()));
        std::fs::write(&path, b">a\nACGT\n>b\nGGTT\n").unwrap();

        let built = Index::build_named(File::open(&path).unwrap(), Format::Fasta).unwrap();
        built.write(&path_for(&path), &path).unwrap();

        let back = Index::read(&path_for(&path)).unwrap();
        assert_eq!(back.offsets(), built.offsets());
        assert_eq!(back.position(b"b").unwrap(), Some(1));

        std::fs::remove_file(&path).ok();
        std::fs::remove_file(path_for(&path)).ok();
    }
}

// ---

/// A file and where each of its records sits in it, with no parser attached.
/// The disk-backed collections each hold one and add their parser; a caller
/// with a parser of its own holds one directly and reads [`record`](Self::record).
#[cfg(any(unix, windows))]
pub struct Framed<S: Source = FileSource> {
    source: S,
    index: Index,
}

#[cfg(any(unix, windows))]
impl Framed<FileSource> {
    /// Index the file at `path` as `format`, parsing none of it.
    pub fn open(path: impl AsRef<Path>, format: Format) -> Result<Self> {
        let source = FileSource::open(path)?;
        let offsets = build_for(source.file(), format)?;

        Ok(Framed {
            source,
            index: Index::new(format, offsets),
        })
    }

    /// Index the file at `path` as whatever format it opens with.
    pub fn detect(path: impl AsRef<Path>) -> Result<Self> {
        let format = crate::format::detect_path(&path)?;

        Self::open(path, format)
    }

    /// Address the file at `path` with `index` instead of scanning for it.
    /// Nothing here checks that the index describes this file; [`is_current`]
    /// does.
    pub fn with_index(path: impl AsRef<Path>, index: Index) -> Result<Self> {
        let source = FileSource::open(path)?;

        Ok(Framed { source, index })
    }
}

#[cfg(feature = "mmap")]
impl Framed<MmapSource> {
    /// Map the file at `path` and index it as `format`, parsing none of it.
    ///
    /// # Safety
    ///
    /// The file must not be modified or truncated while the value is alive;
    /// see [`MmapSource::open`].
    pub unsafe fn map(path: impl AsRef<Path>, format: Format) -> Result<Self> {
        // SAFETY: the caller has promised, through this
        //         function's own contract, not to modify
        //         or truncate the file while the value is
        //         alive
        let source = unsafe { MmapSource::open(path)? };
        let offsets = build_for(source.bytes(), format)?;

        Ok(Framed {
            source,
            index: Index::new(format, offsets),
        })
    }
}

#[cfg(any(unix, windows))]
impl<S: Source> Framed<S> {
    /// The format the records were framed as.
    pub fn format(&self) -> Format {
        self.index.format()
    }

    /// Where each record sits, and what each is called once
    /// [`name_index`](Self::name_index) has run.
    pub fn index(&self) -> &Index {
        &self.index
    }

    /// The number of records framed.
    pub fn len(&self) -> usize {
        self.index.len()
    }

    /// Whether no record was framed.
    pub fn is_empty(&self) -> bool {
        self.index.is_empty()
    }

    /// Name every record, by a byte scan of the range the index frames rather
    /// than a parse of it.
    pub fn name_index(&mut self) -> Result<()> {
        let format = self.index.format();
        let offsets = self.index.nameable()?;
        let mut names = NameBuilder::with_capacity(offsets.len());

        // pushed straight from the range rather than through
        // name_with, whose closure has to return an owned
        // copy of every name
        for o in offsets {
            let bytes = self.source.range(o.start, o.n_bytes)?;
            names.push(crate::seq::name_of(format, &bytes))?;
        }

        self.index.attach(names.finish());

        Ok(())
    }

    /// The `n`th record's framed bytes, or `None` past the end.
    ///
    /// Borrowed from a mapped source and owned from a read one, exactly as
    /// [`Source::range`] hands them over.
    pub fn record(&self, n: usize) -> Result<Option<Cow<'_, [u8]>>> {
        let Some(o) = self.index.offset(n)? else {
            return Ok(None);
        };

        Ok(Some(self.source.range(o.start, o.n_bytes)?))
    }
}

#[cfg(all(test, any(unix, windows)))]
mod framed_tests {
    use super::*;
    use crate::collection::Indexable;
    use crate::fixture;
    use crate::seq::fasta::IndexedFasta;

    const FIXTURES: [(&str, Format); 3] = [
        (fixture!("proteins.fa"), Format::Fasta),
        (fixture!("families.sto"), Format::Stockholm),
        (fixture!("models.hmm"), Format::Hmm),
    ];

    #[test]
    fn a_framed_file_hands_back_the_bytes_the_reader_streams() {
        for (path, format) in FIXTURES {
            let framed = Framed::open(path, format).unwrap();
            let mut reader = Reader::new(File::open(path).unwrap(), format);

            let mut n = 0;
            while reader.advance().unwrap() {
                let bytes = framed.record(n).unwrap().unwrap();
                assert_eq!(bytes.as_ref(), reader.record(), "{path} record {n}");
                assert_eq!(
                    framed.index().offset(n).unwrap().unwrap(),
                    reader.offset(),
                    "{path} record {n}"
                );
                n += 1;
            }

            assert_eq!(n, framed.len(), "{path}");
            assert!(framed.record(n).unwrap().is_none(), "{path}: past the end");
        }
    }

    #[test]
    fn detect_frames_each_fixture_as_its_own_format() {
        for (path, format) in FIXTURES {
            let detected = Framed::detect(path).unwrap();
            let told = Framed::open(path, format).unwrap();

            assert_eq!(detected.format(), format, "{path}");
            assert_eq!(detected.index().offsets(), told.index().offsets(), "{path}");
        }
    }

    #[test]
    fn an_index_handed_over_sets_the_format_and_the_offsets() {
        for (path, format) in FIXTURES {
            let scanned = Framed::open(path, format).unwrap();
            let built = Index::build(File::open(path).unwrap(), format).unwrap();
            let handed = Framed::with_index(path, built).unwrap();

            assert_eq!(handed.format(), format, "{path}");
            assert_eq!(
                handed.index().offsets(),
                scanned.index().offsets(),
                "{path}"
            );
        }
    }

    #[test]
    fn naming_a_framed_file_agrees_with_a_named_build_on_every_format() {
        // name_index reads each record back through the
        // source; build_named names out of the framing
        // buffer. one name per record either way
        for (path, format) in FIXTURES {
            let mut framed = Framed::open(path, format).unwrap();
            framed.name_index().unwrap();

            let built = Index::build_named(File::open(path).unwrap(), format).unwrap();

            assert_eq!(framed.len(), built.len(), "{path}");
            for n in 0..built.len() {
                assert_eq!(
                    framed.index().name(n).unwrap(),
                    built.name(n).unwrap(),
                    "{path} record {n}"
                );
            }
        }
    }

    #[test]
    fn the_fasta_collection_is_a_framed_file_with_a_parser_added() {
        let path = fixture!("proteins.fa");
        let framed = Framed::open(path, Format::Fasta).unwrap();
        let fasta = IndexedFasta::open(path).unwrap();

        assert_eq!(framed.len(), fasta.len());
        for n in 0..fasta.len() {
            assert_eq!(
                framed.record(n).unwrap(),
                fasta.record(n).unwrap(),
                "record {n}"
            );
        }
    }

    #[cfg(feature = "mmap")]
    #[test]
    fn a_mapped_file_frames_the_same_bytes_as_a_read_one() {
        for (path, format) in FIXTURES {
            let read = Framed::open(path, format).unwrap();
            let mapped = unsafe { Framed::map(path, format) }.unwrap();

            assert_eq!(read.len(), mapped.len(), "{path}");
            for n in 0..read.len() {
                let from_map = mapped.record(n).unwrap().unwrap();
                assert!(matches!(from_map, Cow::Borrowed(_)), "{path} record {n}");
                assert_eq!(
                    read.record(n).unwrap().unwrap(),
                    from_map,
                    "{path} record {n}"
                );
            }
        }
    }
}

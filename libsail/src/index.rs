//! Byte offsets of the records in a file.
//!
//! Records in these formats are separated by a fixed pattern at the start of a
//! line — FASTA's `>` opens the record it introduces, Stockholm's `//` closes
//! the one it ends.
//!
//! [`scan`] reads bytes and reports delimiter positions, which is where the
//! chunk-boundary care lives. [`frame`] turns those into [`Offset`]s, which is
//! pure arithmetic and the only part that consults [`Delimiter`].

use std::borrow::Cow;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use memchr::memmem;

use crate::error::{Error, Result, bail_parse};
use crate::format::Format;

/// Where one record sits in the bytes it was indexed from.
#[derive(Default, Clone, Debug, PartialEq)]
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

    /// Fails to compile for a delimiter [`scan`] cannot handle.
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
/// position is tracked by hand, so this indexes a pipe as happily as a file.
pub fn build_index<R, D>(bytes: R) -> Result<Vec<Offset>>
where
    R: Read,
    D: Delimit,
{
    let (positions, total) = scan::<R, D>(bytes, CHUNK)?;

    Ok(frame::<D>(&positions, total))
}

// ---

/// The extension [`Index::path_for`] appends.
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
/// adds the names, and [`save_as`](Index::save_as) writes the pair beside the
/// file they describe.
#[derive(Clone, Debug, PartialEq)]
pub struct Index {
    format: Format,
    offsets: Vec<Offset>,

    /// Every record's name, absent until [`name_with`](Index::name_with) runs.
    //
    // a layer rather than a field per record: an index is
    // built to address a file, and naming it costs a pass
    // over every record's bytes that most callers never want
    names: Option<Names>,
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
#[derive(Clone, Copy, Debug, PartialEq)]
struct Span {
    start: u64,
    len: u32,
}

/// The name `span` frames, or `None` where a record carries none.
fn name_in(blob: &[u8], span: Option<Span>) -> Option<&[u8]> {
    let span = span?;

    Some(&blob[span.start as usize..][..span.len as usize])
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
            offsets,
            names: None,
        }
    }

    pub fn format(&self) -> Format {
        self.format
    }

    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    pub fn offset(&self, n: usize) -> Option<&Offset> {
        self.offsets.get(n)
    }

    pub fn offsets(&self) -> &[Offset] {
        &self.offsets
    }

    pub fn has_names(&self) -> bool {
        self.names.is_some()
    }

    /// The `n`th record's name, `None` when the index is unnamed or that
    /// record carries no name.
    pub fn name(&self, n: usize) -> Option<&[u8]> {
        self.names.as_ref()?.name(n)
    }

    /// The position of the first record called `name`.
    pub fn position(&self, name: &[u8]) -> Option<usize> {
        let names = self.names.as_ref()?;

        // by_name is ordered by name bytes, so the first
        // entry not ordered before `name` is the only one
        // that can equal it. every entry names a record that
        // carries one -- name_with builds it that way and
        // read_from refuses anything else -- so the fallback
        // is unreachable rather than a default
        let at = names
            .by_name
            .partition_point(|&n| names.name(n as usize).unwrap_or_default() < name);

        let n = *names.by_name.get(at)? as usize;

        (names.name(n) == Some(name)).then_some(n)
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
        if u32::try_from(self.offsets.len()).is_err() {
            bail_parse!(
                "an index of {} records is past the u32 a name lookup addresses",
                self.offsets.len()
            )
        }

        let mut blob = Vec::new();
        let mut spans = Vec::with_capacity(self.offsets.len());

        for offset in &self.offsets {
            let Some(name) = name_of(offset)? else {
                spans.push(None);
                continue;
            };

            if name.len() > MAX_NAME {
                bail_parse!(
                    "a record name of {} bytes exceeds the {MAX_NAME} an index holds",
                    name.len()
                )
            }

            spans.push(Some(Span {
                start: blob.len() as u64,
                len: name.len() as u32,
            }));
            blob.extend_from_slice(&name);
        }

        let mut by_name: Vec<u32> = (0..spans.len() as u32)
            .filter(|&n| spans[n as usize].is_some())
            .collect();

        // stable, so equal names stay in file order and
        // position answers with the first of a repeat -- the
        // rule collection::Indexable::index_by follows
        by_name.sort_by(|&a, &b| {
            name_in(&blob, spans[a as usize]).cmp(&name_in(&blob, spans[b as usize]))
        });

        self.names = Some(Names {
            blob,
            spans,
            by_name,
        });

        Ok(())
    }

    /// The path a file's index is written to: the file's own, plus
    /// [`EXTENSION`].
    pub fn path_for(source: impl AsRef<Path>) -> PathBuf {
        // appended rather than substituted, so proteins.fa
        // and proteins.sto do not both index to
        // proteins.saidx
        let mut name = OsString::from(source.as_ref().as_os_str());
        name.push(".");
        name.push(EXTENSION);

        PathBuf::from(name)
    }

    /// Write this index to `at`, stamped with `source`'s length and
    /// modification time so [`load_from`](Index::load_from) can tell it apart
    /// from a file that has since changed.
    pub fn save_as(&self, source: &Path, at: &Path) -> Result<()> {
        let metadata = std::fs::metadata(source)?;

        let mut temporary = OsString::from(at.as_os_str());
        temporary.push(".tmp");
        let temporary = PathBuf::from(temporary);

        // written aside and renamed, because rename is atomic
        // where a write is not: a reader never opens a half
        // written index
        let mut out = BufWriter::new(File::create(&temporary)?);
        self.write_to(&mut out, metadata.len(), metadata.modified()?)?;
        out.flush()?;
        drop(out);

        std::fs::rename(&temporary, at)?;

        Ok(())
    }

    /// [`save_as`](Index::save_as), to [`path_for`](Index::path_for).
    pub fn save(&self, source: &Path) -> Result<()> {
        self.save_as(source, &Index::path_for(source))
    }

    /// Read the index at `at`, and refuse it unless it still describes
    /// `source`.
    ///
    /// # Errors
    ///
    /// A missing index is the not-found [`Io`](crate::ErrorKind::Io) failure
    /// `File::open` reports.
    pub fn load_from(source: &Path, at: &Path) -> Result<Index> {
        let (index, len, mtime) = Index::read_from(BufReader::new(File::open(at)?))?;

        let metadata = std::fs::metadata(source)?;
        if len != metadata.len() || epoch_parts(mtime) != epoch_parts(metadata.modified()?) {
            bail_parse!(
                "{} is an index of {} as it was, and the file has changed since",
                at.display(),
                source.display()
            )
        }

        Ok(index)
    }

    /// [`load_from`](Index::load_from), from [`path_for`](Index::path_for).
    pub fn load(source: &Path) -> Result<Index> {
        Index::load_from(source, &Index::path_for(source))
    }

    fn write_to<W: Write>(&self, mut w: W, len: u64, mtime: SystemTime) -> Result<()> {
        let (secs, nanos) = epoch_parts(mtime);

        w.write_all(&MAGIC)?;
        w.write_all(&[VERSION, format_byte(self.format)])?;
        w.write_all(&len.to_le_bytes())?;
        w.write_all(&secs.to_le_bytes())?;
        w.write_all(&nanos.to_le_bytes())?;
        w.write_all(&(self.offsets.len() as u64).to_le_bytes())?;
        w.write_all(&[u8::from(self.names.is_some())])?;

        for offset in &self.offsets {
            w.write_all(&offset.start.to_le_bytes())?;
            w.write_all(&(offset.n_bytes as u64).to_le_bytes())?;
        }

        let Some(names) = &self.names else {
            return Ok(());
        };

        w.write_all(&(names.blob.len() as u64).to_le_bytes())?;
        w.write_all(&names.blob)?;

        for span in &names.spans {
            let (present, start, len) = match span {
                Some(s) => (1u8, s.start, s.len),
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
                offsets,
                names,
            },
            len,
            UNIX_EPOCH + std::time::Duration::new(secs, nanos),
        ))
    }
}

// ---

fn read_offset<R: Read>(r: &mut R) -> Result<Offset> {
    let start = u64::from_le_bytes(take(r)?);
    let n_bytes = u64::from_le_bytes(take(r)?);

    // usize on a 32-bit target is narrower than what was
    // written, and `as` would truncate a record silently
    let n_bytes = usize::try_from(n_bytes)
        .map_err(|_| Error::parse(format!("a record of {n_bytes} bytes does not fit a usize")))?;

    Ok(Offset { start, n_bytes })
}

fn read_names<R: Read>(r: &mut R, count: usize) -> Result<Names> {
    let blob_len = u64::from_le_bytes(take(r)?);

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

            Ok(Some(Span { start, len }))
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
        // short dense inputs at every small chunk size: the
        // risk buffering buys is a delimiter lost or doubled
        // at a boundary, and a chunk of 2 against a four-byte
        // pattern splits it every way it can be split
        agrees_with_the_definition::<InitiatingDelim>(">");
        agrees_with_the_definition::<TerminatingDelim>("//");
        agrees_with_the_definition::<PairDelim>("@>");
        agrees_with_the_definition::<WideDelim>(">>>");
    }

    /// A [`Read`] that hands over as little as it likes, and sometimes an
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
        // what the seeded newline buys: the first record has
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

    // --- Index

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

    fn clean(path: &Path) {
        std::fs::remove_file(Index::path_for(path)).ok();
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn a_written_index_reads_back_equal() {
        let index = sample();

        assert_eq!(read_back(&written(&index)).unwrap(), index);
    }

    #[test]
    fn an_index_with_no_name_round_trips_the_absence_not_an_empty_string() {
        // a Stockholm alignment with no #=GF ID has no name
        // at all, which is a different fact from being named
        // "". conflating them would make an unnamed alignment
        // answer to the empty string
        let back = read_back(&written(&sample())).unwrap();

        assert_eq!(back.name(0), Some(&b"SH3_1"[..]));
        assert_eq!(back.name(1), None);
        assert_eq!(back.name(2), Some(&b""[..]));
    }

    #[test]
    fn an_unnamed_record_is_not_found_by_any_name() {
        let index = sample();

        assert_eq!(index.position(b"SH3_1"), Some(0));
        assert_eq!(index.position(b""), Some(2));
        assert_eq!(index.position(b"PDZ"), None);
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

        assert_eq!(index.position(b"a"), Some(1));
        assert_eq!(index.position(b"b"), Some(0));
        assert_eq!(index.position(b"c"), Some(3));
    }

    #[test]
    fn an_index_with_no_names_answers_no_lookup_and_still_round_trips() {
        let index = Index::new(Format::Fasta, vec![off(0, 10), off(10, 4)]);

        assert!(!index.has_names());
        assert_eq!(index.name(0), None);
        assert_eq!(index.position(b"a"), None);
        assert_eq!(read_back(&written(&index)).unwrap(), index);
    }

    #[test]
    fn the_offsets_come_back_in_file_order() {
        assert_eq!(
            read_back(&written(&sample())).unwrap().offsets(),
            &[off(0, 2400), off(2400, 2469), off(4869, 0)]
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
            Index::path_for("data/proteins.fa"),
            PathBuf::from("data/proteins.fa.saidx")
        );
        assert_eq!(
            Index::path_for("data/proteins.sto"),
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
        index.save(&path).unwrap();

        let back = Index::load(&path).unwrap();
        assert_eq!(back.format(), Format::Stockholm);
        assert_eq!(back.offsets(), &offsets[..]);
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
        Index::new(Format::Stockholm, offsets).save(&path).unwrap();

        assert!(Index::load(&path).is_ok());

        let mut bytes = std::fs::read(&path).unwrap();
        bytes.extend_from_slice(b"# STOCKHOLM 1.0\n//\n");
        std::fs::write(&path, &bytes).unwrap();

        assert!(Index::load(&path).is_err());

        clean(&path);
    }
}

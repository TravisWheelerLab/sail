//! Byte offsets of the newline-delimited records in a file.
//!
//! The other framing. [`index`](crate::index) finds a fixed pattern at the start
//! of a line, which is how FASTA, Stockholm and p7hmm mark a record; a hit table
//! marks nothing, and every line that is not a comment is a row.
//!
//! That cannot be a [`Delimit`](crate::index::Delimit) impl, for two reasons
//! worth knowing before trying:
//!
//! - A row's delimiter *is* the newline, and `Delimit::VALID` rejects both ways
//!   of saying so — an empty `DELIM`, and a newline inside one.
//! - A comment line is dropped for what it *contains*, and
//!   [`frame`](crate::index) never sees bytes. That split is what makes the
//!   framing testable with no `Read` impl, and it is the thing in the way here.
//!
//! Almost none of `scan`'s difficulty carries over, though: its carry and its
//! `WIDTH - 1` proof exist because a multi-byte pattern can straddle a chunk
//! boundary, and a single byte cannot. This carries two scalars instead.

use std::io::Read;

use memchr::{memchr, memchr_iter};

use crate::error::Result;
use crate::index::Offset;

/// Fresh bytes read per round.
//
// matches index::CHUNK so the two framings issue the same
// number of reads over the same file
const CHUNK: usize = 64 * 1024;

// ---

/// One line at a time, framed out of a stream, into a buffer reused between
/// calls. The rows [`line_index`] would have found, without the offsets.
///
/// ```no_run
/// # use libsail::lines::Rows;
/// # fn main() -> libsail::Result<()> {
/// let mut rows = Rows::new(std::fs::File::open("hits.tbl")?, |b| b == b'#');
/// while rows.advance()? {
///     let row = rows.row();
///     // ...
/// }
/// # Ok(())
/// # }
/// ```
pub struct Rows<R: Read> {
    source: R,
    skip: fn(u8) -> bool,

    /// The file's bytes from `base`: `buf[i]` is at file offset `base + i`.
    buf: Vec<u8>,
    filled: usize,
    base: u64,

    /// Where the line after the current one begins in `buf`.
    cursor: usize,

    /// The current row as a range of `buf`, with no line ending.
    start: usize,
    len: usize,

    /// The whole line the row came from, its line ending included.
    span: usize,

    /// The line `cursor` sits on, counting from 1.
    at: u64,

    /// The line the current row came from.
    line: u64,

    eof: bool,
}

impl<R: Read> Rows<R> {
    /// Walk the lines of `bytes` that `skip` does not reject. `skip` sees a
    /// line's first byte, and is not called for an empty line.
    pub fn new(bytes: R, skip: fn(u8) -> bool) -> Self {
        Rows {
            source: bytes,
            skip,
            buf: vec![0u8; CHUNK],
            filled: 0,
            base: 0,
            cursor: 0,
            start: 0,
            len: 0,
            span: 0,
            at: 1,
            line: 0,
            eof: false,
        }
    }

    /// Walk the lines of `bytes` with a buffer of `chunk` bytes.
    #[cfg(test)]
    fn with_chunk(bytes: R, skip: fn(u8) -> bool, chunk: usize) -> Self {
        let mut rows = Rows::new(bytes, skip);
        rows.buf = vec![0u8; chunk.max(1)];

        rows
    }

    /// The current row's bytes, with no line ending. Empty before the first
    /// [`advance`](Self::advance).
    pub fn row(&self) -> &[u8] {
        &self.buf[self.start..self.start + self.len]
    }

    /// The line the current row came from, counting from 1. Zero before the
    /// first [`advance`](Self::advance).
    pub fn line(&self) -> u64 {
        self.line
    }

    /// Where the whole line sits in the file, its line ending included, which
    /// is the [`Offset`] [`line_index`] reports for it.
    pub fn offset(&self) -> Offset {
        Offset {
            start: self.base + self.start as u64,
            n_bytes: self.span,
        }
    }

    /// Advance to the next row. `false` at the end of the stream.
    pub fn advance(&mut self) -> Result<bool> {
        loop {
            match memchr(b'\n', &self.buf[self.cursor..self.filled]) {
                Some(i) => {
                    let start = self.cursor;
                    let end = start + i;

                    self.cursor = end + 1;

                    if self.take(start, end, self.cursor - start) {
                        return Ok(true);
                    }
                }

                // a file whose last line carries no newline
                None if self.eof => {
                    let (start, end) = (self.cursor, self.filled);
                    self.cursor = self.filled;

                    if start == end {
                        return Ok(false);
                    }

                    return Ok(self.take(start, end, end - start));
                }

                None => self.fill()?,
            }
        }
    }

    /// Hold the line at `start..end` as the current row, unless it is empty or
    /// `skip` rejects it. `span` counts the line ending too.
    fn take(&mut self, start: usize, end: usize, span: usize) -> bool {
        let line = self.at;
        self.at += 1;

        // exactly the line line_index keeps: an empty one is
        // no line at all, and skip decides the rest from the
        // first byte. trimming before this would drop a
        // \r\n line that line_index holds on to
        if end == start || (self.skip)(self.buf[start]) {
            return false;
        }

        // a table written on windows ends its lines \r\n, and
        // the \r belongs to no column
        let end = if self.buf[end - 1] == b'\r' {
            end - 1
        } else {
            end
        };

        self.start = start;
        self.len = end - start;
        self.span = span;
        self.line = line;

        true
    }

    /// Compact what is spent and read more.
    fn fill(&mut self) -> Result<()> {
        if self.eof {
            return Ok(());
        }

        if self.cursor > 0 {
            self.buf.copy_within(self.cursor..self.filled, 0);
            self.filled -= self.cursor;
            self.base += self.cursor as u64;
            self.cursor = 0;
        } else if self.filled == self.buf.len() {
            // a line longer than the buffer, so compacting
            // frees nothing. growing keeps a line spanning k
            // chunks at O(S) rather than O(S^2/chunk)
            self.buf.resize(self.buf.len() * 2, 0);
        }

        let n = crate::index::read_some(&mut self.source, &mut self.buf[self.filled..])?;
        self.filled += n;

        if n == 0 {
            self.eof = true;
        }

        Ok(())
    }
}

// ---

/// Index every line of `bytes` that `skip` does not reject, in file order.
/// `skip` sees a line's **first byte**, and is not called for an empty line.
pub fn line_index<R: Read>(bytes: R, skip: fn(u8) -> bool) -> Result<Vec<Offset>> {
    scan_lines(bytes, CHUNK, skip)
}

/// Index every line, reading `chunk` bytes per round.
fn scan_lines<R: Read>(mut bytes: R, chunk: usize, skip: fn(u8) -> bool) -> Result<Vec<Offset>> {
    // chunk changes how often the buffer refills and nothing
    // else, which is what the chunk-size sweep relies on
    let chunk = chunk.max(1);
    let mut buf = vec![0u8; chunk];
    let mut lines = Vec::new();

    // the whole of the state carried across a refill:
    //
    //     base  = file offset of buf[0]
    //     start = where the line in progress began
    //     first = its first byte, once a chunk has held it
    //
    // first is an Option because a line can begin in the
    // last byte of one chunk and end in the next, so the
    // byte deciding whether to skip it arrives a round
    // before the newline that ends it
    //
    // it is one byte and not the whole line because a line
    // straddling a boundary would then have to be buffered
    // to be handed over, which is the cost this avoids
    let mut base = 0u64;
    let mut start = 0u64;
    let mut first: Option<u8> = None;

    loop {
        let n = crate::index::read_some(&mut bytes, &mut buf)?;
        if n == 0 {
            break;
        }

        for i in memchr_iter(b'\n', &buf[..n]) {
            let end = base + i as u64 + 1;

            // a line that both began and ended in this chunk
            // has had no earlier round to be seen in
            if first.is_none() && start >= base {
                first = Some(buf[(start - base) as usize]);
            }

            keep(&mut lines, start, end, first, skip);

            start = end;
            first = None;
        }

        // a line still in progress: take its first byte while
        // the chunk holding it is still here
        if first.is_none() && start >= base && start - base < n as u64 {
            first = Some(buf[(start - base) as usize]);
        }

        base += n as u64;
    }

    // a last line that ended at EOF rather than at a newline
    keep(&mut lines, start, base, first, skip);

    Ok(lines)
}

/// Record the line `[start, end)` unless it is empty or `skip` rejects it.
fn keep(
    lines: &mut Vec<Offset>,
    start: u64,
    end: u64,
    first: Option<u8>,
    skip: fn(u8) -> bool,
) -> Option<()> {
    // an empty line is one newline and no content, so there
    // is no first byte to ask about
    let n_bytes = end.checked_sub(start)?;
    if n_bytes == 0 || first.is_none_or(|b| b == b'\n' || skip(b)) {
        return None;
    }

    lines.push(Offset {
        start,
        n_bytes: n_bytes as usize,
    });

    Some(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Nothing is a comment, so every non-empty line is a record.
    fn all(_: u8) -> bool {
        false
    }

    /// `#` opens a comment, as every table format in the crate spells it.
    fn hash(b: u8) -> bool {
        b == b'#'
    }

    /// Every line of `input` by the definition: the stretches between newlines,
    /// dropping the empty ones and those `skip` rejects.
    fn lines_naive(input: &[u8], skip: fn(u8) -> bool) -> Vec<Offset> {
        // the specification scan_lines is an optimization
        // of: no buffer and no chunk, so it cannot go wrong
        // at a boundary
        let mut out = Vec::new();
        let mut start = 0usize;

        for i in 0..=input.len() {
            let at_end = i == input.len();
            if !at_end && input[i] != b'\n' {
                continue;
            }

            let end = if at_end { i } else { i + 1 };
            if end > start && !skip(input[start]) && input[start] != b'\n' {
                out.push(Offset {
                    start: start as u64,
                    n_bytes: end - start,
                });
            }
            start = end;
        }

        out
    }

    /// A xorshift64 stream over bytes that make lines, comments and blanks all
    /// common.
    fn soup(seed: u64, len: usize) -> Vec<u8> {
        // why: newlines and '#' have to be frequent enough
        //      that a short input contains a comment, a blank
        //      line and a boundary between them
        const ALPHABET: &[u8] = b"\n\n\n##abc\t";

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

    #[test]
    fn the_scan_finds_exactly_what_the_definition_says_it_should() {
        // short dense inputs at every small chunk size: the
        // risk buffering buys is a line lost or doubled at a
        // boundary, and the first byte of a line arriving one
        // round before the newline that ends it
        for seed in 1..=120u64 {
            let input = soup(seed, seed as usize * 3 % 200);

            for chunk in 1..=9 {
                for skip in [all as fn(u8) -> bool, hash] {
                    assert_eq!(
                        scan_lines(&input[..], chunk, skip).unwrap(),
                        lines_naive(&input, skip),
                        "disagreed at chunk {chunk}, seed {seed}, input {:?}",
                        String::from_utf8_lossy(&input)
                    );
                }
            }
        }
    }

    #[test]
    fn the_scan_agrees_with_the_definition_past_a_chunk() {
        // the small-input sweep never reaches CHUNK, so it
        // never exercises the buffer the shipped code uses
        let input = soup(0xB10, CHUNK * 2 + 517);
        let found = line_index(&input[..], hash).unwrap();

        assert!(found.len() > 1000, "the input has to be dense to be a test");
        assert_eq!(found, lines_naive(&input, hash));
    }

    #[test]
    fn the_chunk_size_never_changes_what_is_found() {
        let input = b"# a header\nrow one\nrow two\n";
        let want = scan_lines(&input[..], input.len() + 1, hash).unwrap();

        assert_eq!(want.len(), 2, "the header is not a row");
        for chunk in 1..=input.len() + 1 {
            assert_eq!(
                scan_lines(&input[..], chunk, hash).unwrap(),
                want,
                "chunk size {chunk} disagreed"
            );
        }
    }

    #[test]
    fn a_comment_block_at_the_end_is_skipped_too() {
        // HMMER writes a "# Program:" block *after* the rows.
        // a rule that only skipped a header would read those
        // as rows and fail on the first one
        let input = b"row one\nrow two\n#\n# Program:  hmmsearch\n# [ok]\n";
        let found = line_index(&input[..], hash).unwrap();

        assert_eq!(found.len(), 2);
        assert_eq!(
            found[0],
            Offset {
                start: 0,
                n_bytes: 8
            }
        );
        assert_eq!(
            found[1],
            Offset {
                start: 8,
                n_bytes: 8
            }
        );
    }

    #[test]
    fn a_line_keeps_its_newline_the_way_a_framed_record_does() {
        // index::frame's terminating arm includes the
        // delimiter's line ending, and these agree so a parser
        // can trim either the same way
        let found = line_index(&b"ab\ncd\n"[..], all).unwrap();
        assert_eq!(found[0].n_bytes, 3, "\"ab\\n\"");
        assert_eq!(found[1].n_bytes, 3);
    }

    #[test]
    fn a_last_line_without_a_newline_is_still_a_line() {
        let found = line_index(&b"ab\ncd"[..], all).unwrap();
        assert_eq!(found.len(), 2);
        assert_eq!(
            found[1],
            Offset {
                start: 3,
                n_bytes: 2
            }
        );
    }

    #[test]
    fn blank_lines_are_not_records() {
        // BLAST separates query blocks with them, and a blank
        // line has no first byte for the predicate to judge
        let found = line_index(&b"\n\nrow\n\n\n"[..], hash).unwrap();
        assert_eq!(found.len(), 1);
        assert_eq!(
            found[0],
            Offset {
                start: 2,
                n_bytes: 4
            }
        );
    }

    #[test]
    fn input_with_no_rows_yields_nothing() {
        for input in [&b""[..], &b"\n\n\n"[..], &b"# only\n# comments\n"[..]] {
            assert!(line_index(input, hash).unwrap().is_empty(), "{input:?}");
        }
    }

    #[test]
    fn the_offsets_tile_the_rows_with_no_gap_or_overlap() {
        let input = b"one\ntwo\nthree\n";
        let found = line_index(&input[..], all).unwrap();

        let mut next = 0u64;
        for o in &found {
            assert_eq!(o.start, next, "rows must be contiguous here");
            assert_eq!(input[(o.start + o.n_bytes as u64 - 1) as usize], b'\n');
            next = o.start + o.n_bytes as u64;
        }
        assert_eq!(next as usize, input.len());
    }

    // ---

    #[test]
    fn rows_frame_what_line_index_finds_on_every_fixture() {
        // two independent readings of one file: line_index
        // collects every offset in a pass, Rows hands them
        // over one at a time, and a disagreement is a
        // framing bug in whichever moved
        for name in [
            "hits.blast.tbl",
            "hits.hmmer.tblout",
            "hits.hmmer.domtblout",
            "hits.nail.tbl",
        ] {
            let path = format!("{}/../fixtures/{name}", env!("CARGO_MANIFEST_DIR"));
            let bytes = std::fs::read(&path).unwrap();
            let want = line_index(&bytes[..], hash).unwrap();

            // every chunk size puts the boundary somewhere
            // else, and one byte at a time puts it between
            // every pair
            for chunk in [1, 2, 7, 64, 4096] {
                let mut rows = Rows::with_chunk(&bytes[..], hash, chunk);
                let mut got = Vec::new();

                while rows.advance().unwrap() {
                    got.push(rows.offset());
                }

                assert_eq!(got, want, "{name} framed differently at chunk {chunk}");
            }
        }
    }

    #[test]
    fn a_row_carries_neither_its_newline_nor_a_windows_return() {
        let input = b"a b\r\nc d\n";
        let mut rows = Rows::new(&input[..], hash);

        assert!(rows.advance().unwrap());
        assert_eq!(rows.row(), b"a b");
        assert_eq!(rows.offset().n_bytes, 5, "though the offset spans both");

        assert!(rows.advance().unwrap());
        assert_eq!(rows.row(), b"c d");
        assert!(!rows.advance().unwrap());
    }

    #[test]
    fn a_last_line_with_no_newline_is_still_a_row() {
        let mut rows = Rows::new(&b"a\nb"[..], hash);

        assert!(rows.advance().unwrap());
        assert!(rows.advance().unwrap());
        assert_eq!(rows.row(), b"b");
        assert!(!rows.advance().unwrap());
    }

    #[test]
    fn a_comment_partway_down_is_skipped_like_one_at_the_head() {
        // tables concatenated from shards repeat their
        // header, so a reader that skipped only a leading
        // block would parse those as rows
        let input = b"# head\na 1\n\n# again\nb 2\n";
        let mut rows = Rows::new(&input[..], hash);

        assert!(rows.advance().unwrap());
        assert_eq!(rows.row(), b"a 1");
        assert_eq!(rows.line(), 2);

        assert!(rows.advance().unwrap());
        assert_eq!(rows.row(), b"b 2");
        assert_eq!(rows.line(), 5, "the line counts lines and not rows");

        assert!(!rows.advance().unwrap());
    }

    #[test]
    fn a_row_longer_than_the_buffer_is_still_one_row() {
        let long = "x".repeat(10_000);
        let input = format!("a\n{long}\nb\n");
        let mut rows = Rows::with_chunk(input.as_bytes(), hash, 8);

        assert!(rows.advance().unwrap());
        assert!(rows.advance().unwrap());
        assert_eq!(rows.row().len(), 10_000);

        assert!(rows.advance().unwrap());
        assert_eq!(rows.row(), b"b");
    }

    #[test]
    fn nothing_rejected_means_every_line_is_a_row() {
        let mut rows = Rows::new(&b"# not a comment here\na\n"[..], all);

        assert!(rows.advance().unwrap());
        assert_eq!(rows.row(), b"# not a comment here");
    }
}

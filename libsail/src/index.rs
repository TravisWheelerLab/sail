//! Byte offsets of the records in a file.
//!
//! Records in these formats are separated by a fixed pattern at the start of a
//! line — FASTA's `>` opens the record it introduces, Stockholm's `//` closes
//! the one it ends.
//!
//! [`scan`] reads bytes and reports delimiter positions, which is where the
//! chunk-boundary care lives. [`frame`] turns those into [`Offset`]s, which is
//! pure arithmetic and the only part that consults [`Delimiter`].

use std::io::Read;

use memchr::memmem;

use crate::error::Result;

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

#[cfg(test)]
mod tests {
    use super::*;

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
}

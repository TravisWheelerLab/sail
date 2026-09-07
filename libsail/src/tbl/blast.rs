//! BLAST tabular layout.

use std::fmt::Display;
use std::io::Write;

use super::{Hit, HitColumns, WriteRow, fields_and_tail, number, row_text};
use crate::Record;
use crate::error::{Result, bail_parse};
use crate::parse::Parse;

/// The default `-outfmt 6` / `-outfmt 7` column set, also what
/// `diamond`, `vsearch`, and `mmseqs convertalis` emit by default.
///
/// ```text
/// | 0 | 1 |   2   |   3   |    4     |    5    |   6   |  7  |   8   |  9  |  10  |   11   |
/// qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
/// ```
///
/// A non-default `-outfmt "6 ..."` column set is a different layout; write an
/// impl of [`HitColumns`] for it.
pub struct BlastTable;

impl HitColumns for BlastTable {
    const QUERY: usize = 0;
    const TARGET: usize = 1;
    const Q_START: Option<usize> = Some(6);
    const Q_END: Option<usize> = Some(7);
    const T_START: Option<usize> = Some(8);
    const T_END: Option<usize> = Some(9);
    const SCORE: usize = 11;
    const E_VALUE: usize = 10;
    const N_FIELDS: usize = 12;
}

// ---

/// One row of [`BlastTable`], every column of it.
#[derive(Debug, Clone, PartialEq)]
pub struct BlastRow {
    pub query: String,       //  0  qseqid
    pub target: String,      //  1  sseqid
    pub identity: f32,       //  2  pident, a percentage
    pub length: usize,       //  3  length, gaps included
    pub mismatches: usize,   //  4  mismatch
    pub gap_opens: usize,    //  5  gapopen
    pub query_start: usize,  //  6  qstart
    pub query_end: usize,    //  7  qend
    pub target_start: usize, //  8  sstart
    pub target_end: usize,   //  9  send
    pub e_value: f64,        // 10  evalue
    pub bit_score: f32,      // 11  bitscore
}

impl BlastRow {
    /// Write this row in [`BlastTable`]'s column order, tab-separated rather
    /// than padded to the widths a tool aligns its columns to.
    pub fn write_to<W: Write>(&self, mut w: W) -> Result<()> {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:e}\t{}",
            self.query,
            self.target,
            self.identity,
            self.length,
            self.mismatches,
            self.gap_opens,
            self.query_start,
            self.query_end,
            self.target_start,
            self.target_end,
            self.e_value,
            self.bit_score,
        )?;

        Ok(())
    }
}

impl Display for BlastRow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut buf = Vec::new();
        self.write_to(&mut buf).map_err(|_| std::fmt::Error)?;
        f.write_str(std::str::from_utf8(&buf).map_err(|_| std::fmt::Error)?)
    }
}

impl Record for BlastRow {}

impl Parse for BlastRow {
    type Record = Self;

    /// The inverse of [`BlastRow::write_to`].
    fn parse(bytes: &[u8]) -> Result<Self> {
        let text = row_text(bytes)?;

        // the tail is dropped, not rejected: `-outfmt "6 std
        // qlen slen"` appends columns after the standard
        // twelve, and the twelve in front are still these
        let (f, _) = fields_and_tail(text, BlastTable::N_FIELDS);
        if f.len() < BlastTable::N_FIELDS {
            bail_parse!(
                "a BLAST row has {} columns, found {}: {:?}",
                BlastTable::N_FIELDS,
                f.len(),
                text.trim_end()
            )
        }

        Ok(BlastRow {
            query: f[0].to_string(),
            target: f[1].to_string(),
            identity: number(f[2], "pident")?,
            length: number(f[3], "length")?,
            mismatches: number(f[4], "mismatch")?,
            gap_opens: number(f[5], "gapopen")?,
            query_start: number(f[6], "qstart")?,
            query_end: number(f[7], "qend")?,
            target_start: number(f[8], "sstart")?,
            target_end: number(f[9], "send")?,
            e_value: number(f[10], "evalue")?,
            bit_score: number(f[11], "bitscore")?,
        })
    }
}

impl WriteRow for BlastRow {
    fn write_row<W: Write>(record: &Self, w: W) -> Result<()> {
        record.write_to(w)
    }
}

impl From<&BlastRow> for Hit {
    fn from(row: &BlastRow) -> Hit {
        Hit {
            query: row.query.clone(),
            target: row.target.clone(),
            query_start: Some(row.query_start),
            query_end: Some(row.query_end),
            target_start: Some(row.target_start),
            target_end: Some(row.target_end),
            score: row.bit_score,
            e_value: row.e_value,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::collection::Indexable;
    use crate::fixture;
    use crate::tbl::BlastRows;

    fn row() -> BlastRow {
        BlastRow {
            query: "q1".into(),
            target: "t1".into(),
            identity: 54.348,
            length: 46,
            mismatches: 21,
            gap_opens: 0,
            query_start: 1,
            query_end: 46,
            target_start: 160,
            target_end: 205,
            e_value: 4.7e-23,
            bit_score: 69.0,
        }
    }

    #[test]
    fn a_written_row_returns_every_column_it_was_given() {
        // unlike a Hit, a row holds all twelve columns, so
        // this round trip loses nothing -- and the reversed
        // coordinates are here because a minus-strand match
        // must survive both directions
        for r in [
            row(),
            BlastRow {
                target_start: 900,
                target_end: 400,
                identity: 100.0,
                e_value: 3.2e-88,
                ..row()
            },
        ] {
            let mut out = Vec::new();
            r.write_to(&mut out).unwrap();

            assert_eq!(BlastRow::parse(&out).unwrap(), r);
            assert_eq!(
                r.to_string().into_bytes(),
                out,
                "Display and write_to agree"
            );
        }
    }

    #[test]
    fn the_first_fixture_row_holds_every_column_the_tool_wrote() {
        // the only check that catches a column read into the
        // wrong field: real output, compared value by value
        let rows = BlastRows::open(fixture!("hits.blast.tbl")).unwrap();
        let first = rows.get(0).unwrap();

        assert_eq!(rows.len(), 4);
        assert_eq!(
            *first,
            BlastRow {
                query: "SH3_1".into(),
                target: "SEM5_CAEEL".into(),
                identity: 54.348,
                length: 46,
                mismatches: 21,
                gap_opens: 0,
                query_start: 1,
                query_end: 46,
                target_start: 160,
                target_end: 205,
                e_value: 4.7e-23,
                bit_score: 69.0,
            }
        );
    }

    #[test]
    fn a_row_too_narrow_for_the_format_says_so() {
        let err = BlastRow::parse(b"q1\tt1\t54.3\n").unwrap_err();

        assert!(err.to_string().contains("12 columns"), "{err}");
        assert!(err.to_string().contains("found 3"), "{err}");
    }

    #[test]
    fn a_column_that_is_not_a_number_names_itself_as_the_format_spells_it() {
        let text = row().to_string().replace("21", "twenty-one");

        let err = BlastRow::parse(text.as_bytes()).unwrap_err();
        assert!(err.to_string().contains("mismatch"), "{err}");
    }
}

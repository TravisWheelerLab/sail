//! nail table layout.

use std::fmt::Display;
use std::io::Write;

use super::{Hit, HitColumns, WriteRow, number, split_columns};
use crate::Record;
use crate::error::Result;
use crate::parse::Parse;

/// ```text
/// |   0   |   1  |   2  |  3  |  4  |  5  |  6  |  7 |   8  |  9  |
/// # target  query  tstart tend  qstart qend  score bias evalue cellfrac
/// ```
pub struct NailTable;

impl HitColumns for NailTable {
    const QUERY: usize = 1;
    const TARGET: usize = 0;
    const Q_START: Option<usize> = Some(4);
    const Q_END: Option<usize> = Some(5);
    const T_START: Option<usize> = Some(2);
    const T_END: Option<usize> = Some(3);
    const SCORE: usize = 6;
    const E_VALUE: usize = 8;
    const N_COLUMNS: usize = 10;
}

// ---

/// One row of [`NailTable`], every column of it.
#[derive(Debug, Clone, PartialEq)]
pub struct NailRow {
    pub target: String,      // 0  target
    pub query: String,       // 1  query
    pub target_start: usize, // 2  target start
    pub target_end: usize,   // 3  target end
    pub query_start: usize,  // 4  query start
    pub query_end: usize,    // 5  query end
    pub score: f32,          // 6  score
    pub bias: f32,           // 7  comp bias
    pub e_value: f64,        // 8  evalue

    /// Fraction of the dynamic-programming matrix the sparse pass computed.
    pub cell_fraction: f32, // 9  cell frac
}

impl NailRow {
    /// Write this row in [`NailTable`]'s column order, tab-separated rather
    /// than padded to the widths a tool aligns its columns to.
    pub fn write_to<W: Write>(&self, mut w: W) -> Result<()> {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:e}\t{}",
            self.target,
            self.query,
            self.target_start,
            self.target_end,
            self.query_start,
            self.query_end,
            self.score,
            self.bias,
            self.e_value,
            self.cell_fraction,
        )?;

        Ok(())
    }
}

impl Display for NailRow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut buf = Vec::new();
        self.write_to(&mut buf).map_err(|_| std::fmt::Error)?;
        f.write_str(std::str::from_utf8(&buf).map_err(|_| std::fmt::Error)?)
    }
}

impl Record for NailRow {}

impl Parse for NailRow {
    type Record = Self;

    /// The inverse of [`NailRow::write_to`].
    fn parse(bytes: &[u8]) -> Result<Self> {
        let (c, _) = split_columns(bytes, NailTable::N_COLUMNS, "a nail row")?;

        Ok(NailRow {
            target: c[0].to_string(),
            query: c[1].to_string(),
            target_start: number(c[2], "target start")?,
            target_end: number(c[3], "target end")?,
            query_start: number(c[4], "query start")?,
            query_end: number(c[5], "query end")?,
            score: number(c[6], "score")?,
            bias: number(c[7], "bias")?,
            e_value: number(c[8], "evalue")?,
            cell_fraction: number(c[9], "cell frac")?,
        })
    }
}

impl WriteRow for NailRow {
    fn write_row<W: Write>(record: &Self, w: W) -> Result<()> {
        record.write_to(w)
    }
}

impl From<&NailRow> for Hit {
    fn from(row: &NailRow) -> Hit {
        Hit {
            query: row.query.clone(),
            target: row.target.clone(),
            query_start: Some(row.query_start),
            query_end: Some(row.query_end),
            target_start: Some(row.target_start),
            target_end: Some(row.target_end),
            score: row.score,
            e_value: row.e_value,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::collection::Indexable;
    use crate::fixture;
    use crate::tbl::NailRows;

    fn row() -> NailRow {
        NailRow {
            target: "t1".into(),
            query: "q1".into(),
            target_start: 160,
            target_end: 205,
            query_start: 1,
            query_end: 46,
            score: 69.0,
            bias: 0.1,
            e_value: 4.7e-23,
            cell_fraction: 0.073,
        }
    }

    #[test]
    fn a_written_row_returns_every_column_it_was_given() {
        // the bias and the cell fraction are the columns a Hit
        // cannot hold, so they are what this is really pinning
        for r in [
            row(),
            NailRow {
                target_start: 900,
                target_end: 400,
                bias: 12.75,
                cell_fraction: 1.0,
                ..row()
            },
        ] {
            let mut out = Vec::new();
            r.write_to(&mut out).unwrap();

            assert_eq!(NailRow::parse(&out).unwrap(), r);
            assert_eq!(
                r.to_string().into_bytes(),
                out,
                "Display and write_to agree"
            );
        }
    }

    #[test]
    fn the_first_fixture_row_holds_every_column_the_tool_wrote() {
        // nail leads with the target, where BLAST leads with
        // the query: read positionally, the two are transposed
        let rows = NailRows::open(fixture!("hits.nail.tbl")).unwrap();

        assert_eq!(rows.len(), 4);
        assert_eq!(
            *rows.get(0).unwrap(),
            NailRow {
                target: "SEM5_CAEEL".into(),
                query: "SH3_1".into(),
                target_start: 160,
                target_end: 205,
                query_start: 1,
                query_end: 46,
                score: 69.0,
                bias: 0.1,
                e_value: 4.7e-23,
                cell_fraction: 0.073,
            }
        );
    }

    #[test]
    fn a_row_too_narrow_for_the_format_says_so() {
        let err = NailRow::parse(b"t1\tq1\t160\t205\n").unwrap_err();

        assert!(err.to_string().contains("10 columns"), "{err}");
        assert!(err.to_string().contains("found 4"), "{err}");
    }
}

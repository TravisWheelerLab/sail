//! HMMER table layouts.

use std::fmt::Display;
use std::io::Write;

use super::{Hit, HitColumns, WriteRow, number, split_columns};
use crate::Record;
use crate::error::Result;
use crate::parse::Parse;

/// `hmmsearch`/`hmmscan` `--tblout`: one row per sequence.
///
/// ```text
/// |      0      |    1     |     2      |     3    |    4   |  5   |  6  | .. |         18          |
/// # target name   accession   query name  accession  E-value  score  bias  ..  description of target
/// ```
///
/// The trailing description is free text containing spaces, so it is not a
/// column and is not captured.
pub struct HmmerTable;

impl HitColumns for HmmerTable {
    const QUERY: usize = 2;
    const TARGET: usize = 0;
    const Q_START: Option<usize> = None;
    const Q_END: Option<usize> = None;
    const T_START: Option<usize> = None;
    const T_END: Option<usize> = None;
    const SCORE: usize = 5;
    const E_VALUE: usize = 4;
    const N_COLUMNS: usize = 18;
}

/// `hmmsearch`/`hmmscan` `--domtblout`: one row per domain.
///
/// ```text
/// |      0      |  1  |  2  |     3     |  4  |  5  |    6    |  7  | .. | 15 | 16 | 17 | 18 | .. |
/// # target name  acc   tlen  query name  acc   qlen  E-value  score  ..  hmm from/to  ali from/to ..
/// ```
///
/// Coordinates are the hmm (query) and ali (target) spans.
pub struct HmmerDomTable;

impl HitColumns for HmmerDomTable {
    const QUERY: usize = 3;
    const TARGET: usize = 0;
    const Q_START: Option<usize> = Some(15);
    const Q_END: Option<usize> = Some(16);
    const T_START: Option<usize> = Some(17);
    const T_END: Option<usize> = Some(18);
    const SCORE: usize = 13;
    const E_VALUE: usize = 12;
    const N_COLUMNS: usize = 22;
}

// ---

/// A column's value, or `None` where HMMER wrote the `-` it fills a column
/// it has nothing for with.
fn reported(value: &str) -> Option<String> {
    match value {
        "" | "-" => None,
        v => Some(v.to_string()),
    }
}

/// The text a column's value is written as: `-` where there is none.
fn as_written(value: &Option<String>) -> &str {
    match value {
        Some(v) => v,
        None => "-",
    }
}

/// One row of [`HmmerTable`], every column of it plus the description.
#[derive(Debug, Clone, PartialEq)]
pub struct HmmerRow {
    pub target: String,                   //  0  target name
    pub target_accession: Option<String>, //  1  accession
    pub query: String,                    //  2  query name
    pub query_accession: Option<String>,  //  3  accession

    // full sequence
    pub e_value: f64, //  4  E-value
    pub score: f32,   //  5  score
    pub bias: f32,    //  6  bias

    // best 1 domain
    pub best_e_value: f64, //  7  E-value
    pub best_score: f32,   //  8  score
    pub best_bias: f32,    //  9  bias

    // domain number estimation, HMMER's own count of what it
    // found, which the cryptic column names do not give away
    pub exp: f32,   // 10  exp  expected number of domains
    pub reg: usize, // 11  reg  regions defined
    pub clu: usize, // 12  clu  regions that looked multidomain and were clustered
    pub ov: usize,  // 13  ov   envelopes overlapping another envelope
    pub env: usize, // 14  env  envelopes defined
    pub dom: usize, // 15  dom  domains defined
    pub rep: usize, // 16  rep  domains above the reporting threshold
    pub inc: usize, // 17  inc  domains above the inclusion threshold

    /// The target's description, which is free text and runs to end of line.
    pub description: Option<String>, // 18
}

impl HmmerRow {
    /// Write this row in [`HmmerTable`]'s column order, tab-separated rather
    /// than padded to the widths HMMER aligns its columns to.
    pub fn write_to<W: Write>(&self, mut w: W) -> Result<()> {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:e}\t{}\t{}\t{:e}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.target,
            as_written(&self.target_accession),
            self.query,
            as_written(&self.query_accession),
            self.e_value,
            self.score,
            self.bias,
            self.best_e_value,
            self.best_score,
            self.best_bias,
            self.exp,
            self.reg,
            self.clu,
            self.ov,
            self.env,
            self.dom,
            self.rep,
            self.inc,
            as_written(&self.description),
        )?;

        Ok(())
    }
}

impl Display for HmmerRow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut buf = Vec::new();
        self.write_to(&mut buf).map_err(|_| std::fmt::Error)?;
        f.write_str(std::str::from_utf8(&buf).map_err(|_| std::fmt::Error)?)
    }
}

impl Record for HmmerRow {}

impl Parse for HmmerRow {
    type Record = Self;

    /// The inverse of [`HmmerRow::write_to`].
    fn parse(bytes: &[u8]) -> Result<Self> {
        let (c, description) = split_columns(bytes, HmmerTable::N_COLUMNS, "a hmmer --tblout row")?;

        Ok(HmmerRow {
            target: c[0].to_string(),
            target_accession: reported(c[1]),
            query: c[2].to_string(),
            query_accession: reported(c[3]),
            e_value: number(c[4], "E-value")?,
            score: number(c[5], "score")?,
            bias: number(c[6], "bias")?,
            best_e_value: number(c[7], "best domain E-value")?,
            best_score: number(c[8], "best domain score")?,
            best_bias: number(c[9], "best domain bias")?,
            exp: number(c[10], "exp")?,
            reg: number(c[11], "reg")?,
            clu: number(c[12], "clu")?,
            ov: number(c[13], "ov")?,
            env: number(c[14], "env")?,
            dom: number(c[15], "dom")?,
            rep: number(c[16], "rep")?,
            inc: number(c[17], "inc")?,
            description: reported(description),
        })
    }
}

impl WriteRow for HmmerRow {
    fn write_row<W: Write>(record: &Self, w: W) -> Result<()> {
        record.write_to(w)
    }
}

impl From<&HmmerRow> for Hit {
    fn from(row: &HmmerRow) -> Hit {
        Hit {
            query: row.query.clone(),
            target: row.target.clone(),
            query_start: None,
            query_end: None,
            target_start: None,
            target_end: None,
            score: row.score,
            e_value: row.e_value,
        }
    }
}

// ---

/// One row of [`HmmerDomTable`], every column of it plus the description.
#[derive(Debug, Clone, PartialEq)]
pub struct HmmerDomRow {
    pub target: String,                   //  0  target name
    pub target_accession: Option<String>, //  1  accession
    pub target_len: usize,                //  2  tlen
    pub query: String,                    //  3  query name
    pub query_accession: Option<String>,  //  4  accession
    pub query_len: usize,                 //  5  qlen

    // full sequence
    pub e_value: f64, //  6  E-value
    pub score: f32,   //  7  score
    pub bias: f32,    //  8  bias

    // this domain
    pub domain_index: usize, //  9  #   this domain's number
    pub domain_count: usize, // 10  of  domains reported for the pair
    pub c_e_value: f64,      // 11  c-Evalue, conditional
    pub i_e_value: f64,      // 12  i-Evalue, independent
    pub domain_score: f32,   // 13  score
    pub domain_bias: f32,    // 14  bias

    pub hmm_from: usize, // 15  hmm coord from
    pub hmm_to: usize,   // 16  hmm coord to
    pub ali_from: usize, // 17  ali coord from
    pub ali_to: usize,   // 18  ali coord to
    pub env_from: usize, // 19  env coord from
    pub env_to: usize,   // 20  env coord to

    /// Mean posterior probability of the aligned residues.
    pub accuracy: f32, // 21  acc

    /// The target's description, which is free text and runs to end of line.
    pub description: Option<String>, // 22
}

impl HmmerDomRow {
    /// Write this row in [`HmmerDomTable`]'s column order, tab-separated
    /// rather than padded to the widths HMMER aligns its columns to.
    pub fn write_to<W: Write>(&self, mut w: W) -> Result<()> {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:e}\t{}\t{}\t{}\t{}\t{:e}\t{:e}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.target,
            as_written(&self.target_accession),
            self.target_len,
            self.query,
            as_written(&self.query_accession),
            self.query_len,
            self.e_value,
            self.score,
            self.bias,
            self.domain_index,
            self.domain_count,
            self.c_e_value,
            self.i_e_value,
            self.domain_score,
            self.domain_bias,
            self.hmm_from,
            self.hmm_to,
            self.ali_from,
            self.ali_to,
            self.env_from,
            self.env_to,
            self.accuracy,
            as_written(&self.description),
        )?;

        Ok(())
    }
}

impl Display for HmmerDomRow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut buf = Vec::new();
        self.write_to(&mut buf).map_err(|_| std::fmt::Error)?;
        f.write_str(std::str::from_utf8(&buf).map_err(|_| std::fmt::Error)?)
    }
}

impl Record for HmmerDomRow {}

impl Parse for HmmerDomRow {
    type Record = Self;

    /// The inverse of [`HmmerDomRow::write_to`].
    fn parse(bytes: &[u8]) -> Result<Self> {
        let (c, description) =
            split_columns(bytes, HmmerDomTable::N_COLUMNS, "a hmmer --domtblout row")?;

        Ok(HmmerDomRow {
            target: c[0].to_string(),
            target_accession: reported(c[1]),
            target_len: number(c[2], "tlen")?,
            query: c[3].to_string(),
            query_accession: reported(c[4]),
            query_len: number(c[5], "qlen")?,
            e_value: number(c[6], "E-value")?,
            score: number(c[7], "score")?,
            bias: number(c[8], "bias")?,
            domain_index: number(c[9], "domain number")?,
            domain_count: number(c[10], "domain count")?,
            c_e_value: number(c[11], "c-Evalue")?,
            i_e_value: number(c[12], "i-Evalue")?,
            domain_score: number(c[13], "domain score")?,
            domain_bias: number(c[14], "domain bias")?,
            hmm_from: number(c[15], "hmm from")?,
            hmm_to: number(c[16], "hmm to")?,
            ali_from: number(c[17], "ali from")?,
            ali_to: number(c[18], "ali to")?,
            env_from: number(c[19], "env from")?,
            env_to: number(c[20], "env to")?,
            accuracy: number(c[21], "acc")?,
            description: reported(description),
        })
    }
}

impl WriteRow for HmmerDomRow {
    fn write_row<W: Write>(record: &Self, w: W) -> Result<()> {
        record.write_to(w)
    }
}

impl From<&HmmerDomRow> for Hit {
    fn from(row: &HmmerDomRow) -> Hit {
        Hit {
            query: row.query.clone(),
            target: row.target.clone(),

            // hmm coordinates are the query's and ali the
            // target's, which is what HmmerDomTable indexes.
            // the env coordinates have nowhere to go: a Hit
            // holds one span per sequence
            query_start: Some(row.hmm_from),
            query_end: Some(row.hmm_to),
            target_start: Some(row.ali_from),
            target_end: Some(row.ali_to),

            // this row describes one domain, so the values
            // that vary across the rows of a query-target pair
            // are the per-domain ones. the independent E-value
            // is the one comparable across searches
            score: row.domain_score,
            e_value: row.i_e_value,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::collection::{Indexable, Iterable};
    use crate::fixture;
    use crate::tbl::{HmmerDomRows, HmmerRows};

    fn row() -> HmmerRow {
        HmmerRow {
            target: "t1".into(),
            target_accession: None,
            query: "q1".into(),
            query_accession: Some("PF00018.35".into()),
            e_value: 1.5e-34,
            score: 105.8,
            bias: 1.1,
            best_e_value: 1.1e-22,
            best_score: 67.8,
            best_bias: 0.1,
            exp: 3.0,
            reg: 3,
            clu: 0,
            ov: 0,
            env: 3,
            dom: 3,
            rep: 2,
            inc: 2,
            description: Some("a described target".into()),
        }
    }

    fn dom_row() -> HmmerDomRow {
        HmmerDomRow {
            target: "t1".into(),
            target_accession: None,
            target_len: 228,
            query: "q1".into(),
            query_accession: Some("PF00018.35".into()),
            query_len: 46,
            e_value: 1.5e-34,
            score: 105.8,
            bias: 1.1,
            domain_index: 1,
            domain_count: 2,
            c_e_value: 1.5e-12,
            i_e_value: 2.5e-12,
            domain_score: 34.7,
            domain_bias: 0.0,
            hmm_from: 2,
            hmm_to: 46,
            ali_from: 5,
            ali_to: 50,
            env_from: 4,
            env_to: 50,
            accuracy: 0.96,
            description: Some("a described target".into()),
        }
    }

    #[test]
    fn a_written_row_returns_every_column_it_was_given() {
        // the second case is the one worth having: an absent
        // accession and an absent description both write as `-`
        // and must read back as None, not as a name of "-"
        for r in [
            row(),
            HmmerRow {
                target_accession: Some("PF00595.30".into()),
                description: None,
                ..row()
            },
        ] {
            let mut out = Vec::new();
            r.write_to(&mut out).unwrap();

            assert_eq!(HmmerRow::parse(&out).unwrap(), r);
            assert_eq!(
                r.to_string().into_bytes(),
                out,
                "Display and write_to agree"
            );
        }
    }

    #[test]
    fn a_written_domain_row_returns_every_column_it_was_given() {
        for r in [
            dom_row(),
            HmmerDomRow {
                domain_index: 2,
                description: None,
                ..dom_row()
            },
        ] {
            let mut out = Vec::new();
            r.write_to(&mut out).unwrap();

            assert_eq!(HmmerDomRow::parse(&out).unwrap(), r);
            assert_eq!(
                r.to_string().into_bytes(),
                out,
                "Display and write_to agree"
            );
        }
    }

    #[test]
    fn a_description_with_spaces_stays_one_field() {
        // the reason these parsers do not just split the row:
        // the description runs to end of line, spaces and all,
        // and no split of the whole row hands it back
        let rows = HmmerRows::open(fixture!("hits.hmmer.tblout")).unwrap();

        assert_eq!(
            rows.get(0).unwrap().description.as_deref(),
            Some("Sex muscle abnormal protein 5")
        );

        // and a '-' is HMMER saying it had nothing to report,
        // not an accession named '-'
        assert!(rows.iter().all(|r| r.target_accession.is_none()));
        assert_eq!(
            rows.get(0).unwrap().query_accession.as_deref(),
            Some("PF00018.35"),
            "the query accession is Pfam's, and is reported"
        );
    }

    #[test]
    fn the_first_fixture_row_holds_every_column_the_tool_wrote() {
        // eighteen columns of real output, value by value. the
        // two E-value/score/bias triples are the pair a reader
        // of the header can most easily swap
        let rows = HmmerRows::open(fixture!("hits.hmmer.tblout")).unwrap();

        assert_eq!(rows.len(), 5);
        assert_eq!(
            *rows.get(0).unwrap(),
            HmmerRow {
                target: "SEM5_CAEEL".into(),
                target_accession: None,
                query: "SH3_1".into(),
                query_accession: Some("PF00018.35".into()),
                e_value: 1.5e-34,
                score: 105.8,
                bias: 1.1,
                best_e_value: 1.1e-22,
                best_score: 67.8,
                best_bias: 0.1,
                exp: 3.0,
                reg: 3,
                clu: 0,
                ov: 0,
                env: 3,
                dom: 3,
                rep: 2,
                inc: 2,
                description: Some("Sex muscle abnormal protein 5".into()),
            }
        );
    }

    #[test]
    fn the_first_fixture_domain_row_holds_every_column_the_tool_wrote() {
        // SEM5's first SH3 domain: 46 model positions against
        // a 228-residue protein, the weaker of its two hits
        let rows = HmmerDomRows::open(fixture!("hits.hmmer.domtblout")).unwrap();

        assert_eq!(rows.len(), 7);
        assert_eq!(
            *rows.get(0).unwrap(),
            HmmerDomRow {
                target: "SEM5_CAEEL".into(),
                target_accession: None,
                target_len: 228,
                query: "SH3_1".into(),
                query_accession: Some("PF00018.35".into()),
                query_len: 46,
                e_value: 1.5e-34,
                score: 105.8,
                bias: 1.1,
                domain_index: 1,
                domain_count: 2,
                c_e_value: 1.5e-12,
                i_e_value: 2.5e-12,
                domain_score: 34.7,
                domain_bias: 0.0,
                hmm_from: 2,
                hmm_to: 46,
                ali_from: 5,
                ali_to: 50,
                env_from: 4,
                env_to: 50,
                accuracy: 0.96,
                description: Some("Sex muscle abnormal protein 5".into()),
            }
        );
    }

    #[test]
    fn the_two_rows_of_a_tandem_domain_differ_only_where_the_domain_does() {
        // SEM5's two SH3 domains share every whole-sequence
        // column and differ in every per-domain one, which is
        // what makes --domtblout worth parsing separately
        let rows = HmmerDomRows::open(fixture!("hits.hmmer.domtblout")).unwrap();
        let (first, second) = (rows.get(0).unwrap(), rows.get(1).unwrap());

        assert_eq!(first.e_value, second.e_value);
        assert_eq!(first.score, second.score);
        assert_eq!(first.target_len, second.target_len);

        assert_eq!((first.domain_index, second.domain_index), (1, 2));
        assert_eq!(first.domain_count, second.domain_count);
        assert!(second.domain_score > first.domain_score);
        assert!(second.ali_from > first.ali_to, "the domains do not overlap");
    }

    #[test]
    fn a_row_too_narrow_for_the_format_says_so() {
        // a --tblout row read as a --domtblout one: eighteen
        // columns and a three-word description, so 21 fields
        // where 22 are needed
        let tblout = row().to_string();

        let err = HmmerDomRow::parse(tblout.as_bytes()).unwrap_err();
        assert!(err.to_string().contains("22 columns"), "{err}");
        assert!(err.to_string().contains("found 21"), "{err}");
    }
}

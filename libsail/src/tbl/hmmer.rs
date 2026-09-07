//! HMMER table layouts.

use super::HitColumns;

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
    const N_FIELDS: usize = 18;
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
    const N_FIELDS: usize = 22;
}

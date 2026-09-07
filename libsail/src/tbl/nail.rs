//! nail table layout.

use super::HitColumns;

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
    const N_FIELDS: usize = 10;
}

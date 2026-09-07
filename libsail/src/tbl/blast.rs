//! BLAST tabular layout.

use super::HitColumns;

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

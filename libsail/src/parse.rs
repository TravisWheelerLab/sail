//! Turning bytes back into records.

use crate::Record;
use crate::error::Result;

/// A trait to read one record of a format out of the bytes it occupies.
pub trait Parse {
    type Record: Record;

    /// The record held in `bytes`.
    //
    // no receiver, so a disk-backed Indexable::get -- which
    // takes &self -- can call it without a RefCell or a
    // Mutex
    //
    // an error rather than an Option: across millions of
    // records, knowing that one failed is not enough to
    // act on
    fn parse(bytes: &[u8]) -> Result<Self::Record>;
}

//! Sequence, alignment, and profile records: one type per format.
//!
//! A format module also owns its [`Delimit`](crate::index::Delimit) impl, since
//! where a record starts and stops in a file is part of the format, and its
//! in-memory and disk-backed collections beside them.

pub mod fasta;
pub mod p7hmm;
pub mod stockholm;

use crate::format::Format;

/// The name in one framed record's bytes, or `None` for a record carrying
/// none, such as a Stockholm alignment with no `#=GF ID`.
pub fn name_of(format: Format, record: &[u8]) -> Option<&[u8]> {
    match format {
        Format::Fasta => fasta::name_of(record),
        Format::Stockholm => stockholm::name_of(record),
        Format::Hmm => p7hmm::name_of(record),
    }
}

/// Decode `bytes`, naming what they were so a failure says where it happened.
pub(crate) fn as_str<'a>(bytes: &'a [u8], what: &str) -> crate::Result<&'a str> {
    match std::str::from_utf8(bytes) {
        Ok(s) => Ok(s),

        // an error rather than replacement characters: a
        // mangled name leaves the record unfindable by the
        // name it has
        Err(e) => crate::error::bail_parse!("{what} is not valid utf-8: {e}"),
    }
}

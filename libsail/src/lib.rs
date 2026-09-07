//! Sequence analysis IO for the profile-HMM and homology-search formats:
//! FASTA, Stockholm, HMMER3 p7hmm, and the BLAST/nail/HMMER hit tables.
//!
//! Organised by concern rather than by format: [`seq`] and [`tbl`] hold the
//! record types, [`index`] finds record boundaries in a file, [`collection`]
//! addresses records by position, and [`parse`] reads text back into records.

pub mod collection;
pub mod error;
pub mod index;
pub mod lines;
pub mod parse;
#[cfg(any(unix, windows))]
mod pread;
pub mod seq;
#[cfg(any(unix, windows))]
pub mod source;
pub mod tbl;

pub use error::{Error, ErrorKind, Result};

/// The path to a file in `fixtures/`, at the workspace root.
//
// cargo runs a test binary with its working directory set to
// the *package* root, so a relative "fixtures/..." would look
// inside libsail/ and find nothing. resolving from the
// manifest instead reaches the workspace root from either
// place, and does it at compile time
#[cfg(test)]
macro_rules! fixture {
    ($name:literal) => {
        concat!(env!("CARGO_MANIFEST_DIR"), "/../fixtures/", $name)
    };
}

#[cfg(test)]
pub(crate) use fixture;

/// One record of some format.
//
// at the root rather than in a module because both
// collection::Indexable and parse::Parse are bound by it
//
// Clone is what Indexable::cloned needs: get lends a record
// where the collection holds one, so asking for it owned has
// to copy it
pub trait Record: Clone {}

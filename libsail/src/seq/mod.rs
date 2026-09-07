//! Sequence, alignment, and profile records: one type per format.
//!
//! A format module also owns its [`Delimit`](crate::index::Delimit) impl, since
//! where a record starts and stops in a file is part of the format, and its
//! in-memory and disk-backed collections beside them.

pub mod fasta;
pub mod p7hmm;
pub mod stockholm;

//! Which format a file holds, decided from its first bytes.
//!
//! Each magic checked here is the one that format's own parser already
//! demands, so a file this accepts is a file that format can go on to read.

use std::fmt::Display;
use std::fs::File;
use std::path::Path;

use crate::error::{Result, bail_parse};
use crate::index::read_some;
use crate::seq::stockholm;

/// Bytes of a file [`detect_path`] reads to identify it.
//
// every magic sits in the first line, and the longest is
// "# STOCKHOLM 1.0"
const PREFIX: usize = 64;

/// One of the record formats this crate reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Stockholm,
    Hmm,
}

impl Format {
    /// The format's name, as a caller naming one would spell it.
    pub fn name(&self) -> &'static str {
        match self {
            Format::Fasta => "fasta",
            Format::Stockholm => "stockholm",
            Format::Hmm => "hmm",
        }
    }
}

impl Display for Format {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

// ---

/// The format that `bytes` opens, where `bytes` starts at the head of a file.
pub fn detect(bytes: &[u8]) -> Result<Format> {
    let start = bytes.trim_ascii_start();

    if start.starts_with(b">") {
        Ok(Format::Fasta)
    } else if start.starts_with(stockholm::HEADER.as_bytes()) {
        Ok(Format::Stockholm)
    } else if start.starts_with(b"HMMER3") {
        Ok(Format::Hmm)
    } else {
        bail_parse!(
            "not fasta, stockholm or hmm: the file opens {:?}",
            String::from_utf8_lossy(&bytes[..bytes.len().min(16)])
        )
    }
}

/// The format of the file at `path`, read from its head.
///
/// The file is opened and closed here, so a caller that goes on to read it
/// starts from byte 0 of its own handle.
pub fn detect_path(path: impl AsRef<Path>) -> Result<Format> {
    let mut head = [0; PREFIX];
    let filled = read_some(&mut File::open(path)?, &mut head)?;

    detect(&head[..filled])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixture;

    #[test]
    fn each_fixtures_own_format_is_read_from_its_first_bytes() {
        // the fixtures are real tool output, so this checks
        // the magics against what the tools emit rather than
        // against what the spec says
        for (path, want) in [
            (fixture!("proteins.fa"), Format::Fasta),
            (fixture!("families.sto"), Format::Stockholm),
            (fixture!("models.hmm"), Format::Hmm),
        ] {
            assert_eq!(detect_path(path).unwrap(), want, "{path}");
        }
    }

    #[test]
    fn a_file_shorter_than_the_prefix_is_still_identified() {
        // read_some fills what it can and stops at EOF, so a
        // file of four bytes must not be waiting on 64
        let path = std::env::temp_dir().join(format!("libsail-short-{}.fa", std::process::id()));
        std::fs::write(&path, b">a\nA").unwrap();

        assert_eq!(detect_path(&path).unwrap(), Format::Fasta);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn bytes_matching_none_of_the_three_magics_are_rejected() {
        assert!(detect(b"CLUSTAL W (1.81) multiple sequence alignment\n").is_err());
        assert!(detect(b"@read1\nACGT\n+\n!!!!\n").is_err());
    }

    #[test]
    fn an_empty_input_cannot_be_identified() {
        assert!(detect(b"").is_err());
    }

    #[test]
    fn leading_blank_lines_do_not_hide_the_format() {
        assert_eq!(detect(b"\n\n>a\nACGT\n").unwrap(), Format::Fasta);
    }

    #[test]
    fn a_file_that_is_not_there_fails_as_io_and_not_as_a_format() {
        let error = detect_path("no/such/file.fa").unwrap_err();

        assert!(matches!(error.kind(), crate::ErrorKind::Io { .. }));
    }
}

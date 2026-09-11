//! What an operation reads.
//!
//! An operation takes one or more files of a single format, or `-` for stdin,
//! and every one is buffered whole: the in-memory collections read to end
//! anyway, and a pipe cannot be identified without reading its head.

use std::io::{IsTerminal, Read};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use libsail::format::{self, Format};

use crate::cli::{Axis, FormatArg};

/// One file's bytes, and what to call it in a message.
pub struct Input {
    name: String,
    bytes: Vec<u8>,
}

impl Input {
    /// The path this was read from, or `<stdin>`.
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }
}

/// Everything one operation reads, all of a single format.
pub struct Inputs {
    format: Format,
    entries: Vec<Input>,
}

impl Inputs {
    /// Read every path, identify the format, and fail when the inputs disagree
    /// with each other or with `asserted`.
    pub fn read(paths: &[PathBuf], asserted: Option<FormatArg>) -> Result<Inputs> {
        if paths.is_empty() {
            bail!("no input given");
        }

        let mut entries: Vec<Input> = Vec::with_capacity(paths.len());
        let mut format = None;

        for path in paths {
            let entry = read_one(path)?;
            let found = format::detect(entry.bytes())
                .with_context(|| format!("identifying {}", entry.name()))?;

            // the first input settles it, so the message can
            // name the one that disagreed and the one it
            // disagreed with
            match format {
                Some(first) if first != found => bail!(
                    "{} reads as {found} and {} reads as {first}: \
                     one operation reads one format",
                    entry.name(),
                    entries[0].name()
                ),
                _ => format = Some(found),
            }

            entries.push(entry);
        }

        let format = format.expect("a non-empty path list settles a format");
        assert_format(format, asserted, "the input")?;

        Ok(Inputs { format, entries })
    }

    pub fn format(&self) -> Format {
        self.format
    }

    pub fn entries(&self) -> &[Input] {
        &self.entries
    }
}

// ---

/// The format of the file at `path`, for an operation that needs the file
/// rather than its bytes.
pub fn format_of_path(path: &Path, asserted: Option<FormatArg>) -> Result<Format> {
    let found =
        format::detect_path(path).with_context(|| format!("identifying {}", path.display()))?;
    assert_format(found, asserted, &path.display().to_string())?;

    Ok(found)
}

/// Fail when `asserted` names a format the input is not.
//
// one wording for every operation: --format is an assertion
// about a file and never a conversion of one, so a
// disagreement is the caller's mistake rather than a
// fallback to what they asked for
fn assert_format(found: Format, asserted: Option<FormatArg>, what: &str) -> Result<()> {
    if let Some(asserted) = asserted.map(Format::from)
        && asserted != found
    {
        bail!("--format {asserted} was given, but {what} reads as {found}");
    }

    Ok(())
}

/// Refuse stdin at a terminal, where the read would block with no sign of why.
//
// fifteen of the nineteen operations default their input to
// `-`, so a bare `sail count` at a prompt would otherwise
// sit on a read the caller never meant to start. taken as a
// bool rather than read here, so the decision is reachable
// from a test -- nothing can make a test binary's stdin a
// terminal
fn usable_stdin(is_terminal: bool) -> Result<()> {
    if is_terminal {
        bail!("no input: give a path, or pipe something in")
    }

    Ok(())
}

fn read_one(path: &Path) -> Result<Input> {
    if path.as_os_str() == "-" {
        usable_stdin(std::io::stdin().is_terminal())?;

        let mut bytes = Vec::new();
        std::io::stdin()
            .read_to_end(&mut bytes)
            .context("reading stdin")?;

        return Ok(Input {
            name: "<stdin>".to_string(),
            bytes,
        });
    }

    let bytes = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;

    Ok(Input {
        name: path.display().to_string(),
        bytes,
    })
}

// ---

/// What one record's SIZE is, which differs per format.
//
// here rather than in each operation because it is one rule
// three of them share: a FASTA record measures in residues,
// a profile in LENG, and an alignment in whichever axis
// `--by` named
#[macro_export]
macro_rules! size_of {
    (Fasta, $axis:expr) => {
        |record: &libsail::seq::fasta::FastaRecord| record.len()
    };
    (Stockholm, $axis:expr) => {
        move |record: &libsail::seq::stockholm::StockholmRecord| match $axis {
            $crate::cli::Axis::Depth => record.depth(),
            $crate::cli::Axis::Width => record.width(),
        }
    };
    (Hmm, $axis:expr) => {
        |record: &libsail::seq::p7hmm::HmmRecord| record.header.leng
    };
}

/// Run `body` against the in-memory collection the format calls for, and
/// against the size function that format's records are measured by.
//
// a macro rather than a generic function because the three
// collections are three unrelated types. Indexable's Yield
// GAT rules out `dyn`, so the type has to be chosen where it
// is written down -- and `body` is usually a call to a
// generic function, which is where the shared code lives
macro_rules! dispatch {
    ($format:expr, $input:expr, |$collection:ident| $body:expr) => {
        $crate::input::dispatch!(
            $format,
            $input,
            $crate::cli::Axis::Depth,
            |$collection, _size, _name| $body
        )
    };

    ($format:expr, $input:expr,
     |$collection:ident, $size:pat_param, $name:pat_param| $body:expr) => {
        $crate::input::dispatch!(
            $format,
            $input,
            $crate::cli::Axis::Depth,
            |$collection, $size, $name| $body
        )
    };

    ($format:expr, $input:expr, $axis:expr,
     |$collection:ident, $size:pat_param, $name:pat_param| $body:expr) => {
        match $format {
            libsail::format::Format::Fasta => {
                let $collection = libsail::seq::fasta::Fasta::new($input.bytes())?;
                let $size = $crate::size_of!(Fasta, $axis);
                let $name = $crate::output::name::fasta;
                $body
            }
            libsail::format::Format::Stockholm => {
                let $collection = libsail::seq::stockholm::Stockholm::new($input.bytes())?;
                let $size = $crate::size_of!(Stockholm, $axis);
                let $name = $crate::output::name::stockholm;
                $body
            }
            libsail::format::Format::Hmm => {
                let $collection = libsail::seq::p7hmm::Hmm::new($input.bytes())?;
                let $size = $crate::size_of!(Hmm, $axis);
                let $name = $crate::output::name::hmm;
                $body
            }
        }
    };
}

pub(crate) use dispatch;

/// Fail when `--by` was given for a format with only one axis to measure.
pub fn axis_for(format: Format, by: Option<Axis>) -> Result<Axis> {
    match (by, format) {
        (Some(_), Format::Fasta | Format::Hmm) => {
            bail!("--by names one of an alignment's two axes, and {format} records have one size")
        }
        (Some(axis), Format::Stockholm) => Ok(axis),
        (None, _) => Ok(Axis::Depth),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(bytes: &[u8], tag: &str, extension: &str) -> PathBuf {
        let path =
            std::env::temp_dir().join(format!("sail-in-{tag}-{}.{extension}", std::process::id()));
        std::fs::write(&path, bytes).unwrap();

        path
    }

    fn fixture(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../fixtures")
            .join(name)
    }

    #[test]
    fn the_format_is_read_from_the_bytes_rather_than_the_extension() {
        // the file is named .txt and holds an alignment, so
        // only its content can be answering
        let path = write(b"# STOCKHOLM 1.0\nseq1 AC\n//\n", "ext", "txt");
        let inputs = Inputs::read(std::slice::from_ref(&path), None).unwrap();

        assert_eq!(inputs.format(), Format::Stockholm);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn several_files_of_one_format_read_as_one_operations_input() {
        let inputs = Inputs::read(&[fixture("proteins.fa"), fixture("proteins.fa")], None).unwrap();

        assert_eq!(inputs.format(), Format::Fasta);
        assert_eq!(inputs.entries().len(), 2);
    }

    #[test]
    fn a_file_of_another_format_is_refused_and_the_message_names_both() {
        // an operation reads one format, so the mismatch has
        // to be caught before any of it is parsed
        let Err(error) = Inputs::read(&[fixture("proteins.fa"), fixture("models.hmm")], None)
        else {
            panic!("two formats read as one operation's input")
        };
        let error = error.to_string();

        assert!(error.contains("models.hmm"), "{error}");
        assert!(error.contains("proteins.fa"), "{error}");
    }

    #[test]
    fn a_mismatched_format_assertion_is_refused() {
        assert!(Inputs::read(&[fixture("proteins.fa")], Some(FormatArg::Stockholm)).is_err());
        assert!(Inputs::read(&[fixture("proteins.fa")], Some(FormatArg::Fasta)).is_ok());
    }

    #[test]
    fn a_path_is_identified_without_reading_all_of_it() {
        assert_eq!(
            format_of_path(&fixture("models.hmm"), None).unwrap(),
            Format::Hmm
        );
        assert!(format_of_path(&fixture("models.hmm"), Some(FormatArg::Fasta)).is_err());
    }

    #[test]
    fn no_input_at_all_is_an_error_rather_than_an_empty_run() {
        assert!(Inputs::read(&[], None).is_err());
    }

    #[test]
    fn stdin_at_a_terminal_is_refused_rather_than_read() {
        // the read would block with nothing on screen to say
        // why. a pipe and a redirect are not terminals, so
        // only the bare interactive call changes
        assert!(usable_stdin(true).is_err());
        assert!(usable_stdin(false).is_ok());
    }
}

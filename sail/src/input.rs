//! What an operation reads.
//!
//! An operation takes one or more files of a single format, or `-` for stdin.
//!
//! A file is identified from its first 64 bytes and then handed to the
//! collection as an open [`File`], so it is read once rather than buffered here
//! and copied again by the collection. A pipe cannot be reopened, so stdin is
//! read whole and held here.

use std::fs::File;
use std::io::{Cursor, IsTerminal, Read};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use libsail::format::{self, Format};

use crate::cli::{Axis, FormatArg, Mode, ReadArgs};

/// One input, identified but not yet read.
pub struct Input {
    name: String,
    source: Source,
}

/// Where an input's bytes come from.
enum Source {
    File(PathBuf),

    // read here rather than left where it came from,
    // because a pipe cannot be reopened
    Stdin(Vec<u8>),
}

impl Input {
    /// The path this was read from, or `<stdin>`.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Whether this input is the standard input.
    pub fn is_stdin(&self) -> bool {
        matches!(self.source, Source::Stdin(_))
    }

    /// The path this input names, or `None` for stdin.
    pub fn path(&self) -> Option<&Path> {
        match &self.source {
            Source::File(path) => Some(path),
            Source::Stdin(_) => None,
        }
    }

    /// A reader over this input, positioned at byte 0.
    pub fn reader(&self) -> Result<Box<dyn Read + '_>> {
        // boxed so the two arms have one type: one virtual
        // call per 64 KiB read, against a second copy of the
        // whole input
        match &self.source {
            Source::File(path) => Ok(Box::new(
                File::open(path).with_context(|| format!("reading {}", path.display()))?,
            )),
            Source::Stdin(bytes) => Ok(Box::new(Cursor::new(&bytes[..]))),
        }
    }
}

/// How an operation's input is actually read, with [`Mode::Auto`] settled.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Backend {
    Stream,
    Memory,
    Indexed,
}

/// What an operation needs of its input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Needs {
    /// one forward pass over the records, in order
    Pass,
    /// the record count before the first record is written
    Count,
    /// a few records, by name or by position
    Sparse,
    /// every record at once, in an order the file does not hold
    Whole,
}

/// Everything one operation reads, all of a single format.
pub struct Inputs {
    format: Format,
    entries: Vec<Input>,
    backend: Backend,
}

impl Inputs {
    /// [`plan`](Self::plan) for an operation that holds every record.
    pub fn read(paths: &[PathBuf], asserted: Option<FormatArg>) -> Result<Inputs> {
        // Whole rather than Pass: these are the operations
        // with no streaming body, and a backend they never
        // consult would be a claim the code does not keep
        Inputs::plan(paths, asserted, ReadArgs::default(), Needs::Whole)
    }

    pub fn backend(&self) -> Backend {
        self.backend
    }

    /// Identify every path, settle one format, and settle one backend.
    pub fn plan(
        paths: &[PathBuf],
        asserted: Option<FormatArg>,
        read: ReadArgs,
        needs: Needs,
    ) -> Result<Inputs> {
        let mut inputs = Inputs::identify(paths, asserted)?;
        inputs.backend = resolve(read.mode, needs, &inputs.entries)?;

        Ok(inputs)
    }

    /// Identify every path, settle one format, and fail when the inputs
    /// disagree with each other or with `asserted`.
    fn identify(paths: &[PathBuf], asserted: Option<FormatArg>) -> Result<Inputs> {
        if paths.is_empty() {
            bail!("no input given");
        }

        let mut entries: Vec<Input> = Vec::with_capacity(paths.len());
        let mut format = None;

        for path in paths {
            let (entry, found) = identify(path)?;

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

        Ok(Inputs {
            format,
            entries,
            backend: Backend::Memory,
        })
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

/// Settle `--read auto` against what the operation needs and what it was given.
fn resolve(mode: Mode, needs: Needs, entries: &[Input]) -> Result<Backend> {
    let piped = entries.iter().any(Input::is_stdin);

    // one table rather than a decision in each operation:
    // "fastest by default" is a claim about these arms, and
    // spread over seventeen run methods it is not a claim
    // anyone can check
    match (mode, needs) {
        // --read asserts a method the way --format asserts a
        // content type. neither is a hint and neither falls
        // back: falling back to memory on a pipe does the one
        // thing the caller was avoiding, and falling back to
        // streaming cannot serve an operation that orders
        // every record
        (Mode::Indexed, _) if piped => bail!(
            "--read indexed addresses records at byte offsets in a file, and <stdin> has none: \
             give a path, or --read stream"
        ),
        (Mode::Stream, Needs::Whole) => bail!(
            "this operation puts the records in an order the file does not hold, so it cannot \
             read them one at a time: --read memory, or --read indexed for a file larger than memory"
        ),

        (Mode::Stream, _) => Ok(Backend::Stream),
        (Mode::Memory, _) => Ok(Backend::Memory),
        (Mode::Indexed, _) => Ok(Backend::Indexed),

        // the defaults: the fastest backend that can answer,
        // given what the operation needs and whether the
        // input is seekable
        (Mode::Auto, Needs::Pass) => Ok(Backend::Stream),
        (Mode::Auto, Needs::Whole) => Ok(Backend::Memory),
        (Mode::Auto, Needs::Count | Needs::Sparse) if piped => Ok(Backend::Stream),
        (Mode::Auto, Needs::Count | Needs::Sparse) => Ok(Backend::Indexed),
    }
}

/// Fail when `asserted` names a format the input is not.
fn assert_format(found: Format, asserted: Option<FormatArg>, what: &str) -> Result<()> {
    if let Some(asserted) = asserted.map(Format::from)
        && asserted != found
    {
        // --format is an assertion about a file and never a
        // conversion of one, so a disagreement is the
        // caller's mistake rather than a fallback to what
        // they asked for
        bail!("--format {asserted} was given, but {what} reads as {found}");
    }

    Ok(())
}

/// Refuse stdin at a terminal, where the read would block with no sign of why.
fn usable_stdin(is_terminal: bool) -> Result<()> {
    // fifteen of the nineteen operations default their
    // input to `-`, so a bare `sail count` at a prompt
    // would otherwise sit on a read the caller never meant
    // to start
    if is_terminal {
        bail!("no input: give a path, or pipe something in")
    }

    Ok(())
}

/// One input and the format it opens with.
fn identify(path: &Path) -> Result<(Input, Format)> {
    if path.as_os_str() == "-" {
        // is_terminal read here and passed in, so the
        // refusal is reachable from a test -- nothing can
        // make a test binary's stdin a terminal
        usable_stdin(std::io::stdin().is_terminal())?;

        let mut bytes = Vec::new();
        std::io::stdin()
            .read_to_end(&mut bytes)
            .context("reading stdin")?;

        let found = format::detect(&bytes).context("identifying <stdin>")?;

        return Ok((
            Input {
                name: "<stdin>".to_string(),
                source: Source::Stdin(bytes),
            },
            found,
        ));
    }

    let found =
        format::detect_path(path).with_context(|| format!("identifying {}", path.display()))?;

    Ok((
        Input {
            name: path.display().to_string(),
            source: Source::File(path.to_path_buf()),
        },
        found,
    ))
}

// ---

/// What one record's SIZE is, which differs per format.
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

/// One framed record's SIZE, the streaming counterpart of [`size_of!`].
pub fn size_framed(format: Format, record: &[u8], axis: Axis) -> Result<usize> {
    use libsail::parse::Parse;

    Ok(match format {
        // FASTA reads its length off the bytes; the other
        // two parse, because their size is a property of a
        // structure rather than a byte count, and a file
        // holds thousands of those where it holds millions
        // of FASTA records
        Format::Fasta => libsail::seq::fasta::len_of(record),
        Format::Stockholm => {
            let record = libsail::seq::stockholm::StockholmParser::parse(record)?;

            match axis {
                Axis::Depth => record.depth(),
                Axis::Width => record.width(),
            }
        }
        Format::Hmm => libsail::seq::p7hmm::HmmParser::parse(record)?.header.leng,
    })
}

/// Run `body` against the in-memory collection the format calls for, and
/// against the size function that format's records are measured by.
macro_rules! dispatch {
    ($format:expr, $input:expr, |$collection:ident| $body:expr) => {
        $crate::input::dispatch!(
            $format,
            $input,
            $crate::cli::Axis::Depth,
            |$collection, _size, _name, _write| $body
        )
    };

    ($format:expr, $input:expr,
     |$collection:ident, $size:pat_param, $name:pat_param, $write:pat_param| $body:expr) => {
        $crate::input::dispatch!(
            $format,
            $input,
            $crate::cli::Axis::Depth,
            |$collection, $size, $name, $write| $body
        )
    };

    ($format:expr, $input:expr, $axis:expr,
     |$collection:ident, $size:pat_param, $name:pat_param, $write:pat_param| $body:expr) => {
        // a macro rather than a generic function: the three
        // collections are unrelated types, and Indexable's
        // Yield GAT rules out `dyn`, so the type has to be
        // chosen where it is written down. `body` is
        // usually a call to a generic function, which is
        // where the shared code lives
        match $format {
            libsail::format::Format::Fasta => {
                let $collection = libsail::seq::fasta::Fasta::new($input.reader()?)?;
                let $size = $crate::size_of!(Fasta, $axis);
                let $name = $crate::output::name::fasta;
                let $write = $crate::output::write::fasta
                    as $crate::output::Writer<libsail::seq::fasta::FastaRecord>;
                $body
            }
            libsail::format::Format::Stockholm => {
                let $collection = libsail::seq::stockholm::Stockholm::new($input.reader()?)?;
                let $size = $crate::size_of!(Stockholm, $axis);
                let $name = $crate::output::name::stockholm;
                let $write = $crate::output::write::stockholm
                    as $crate::output::Writer<libsail::seq::stockholm::StockholmRecord>;
                $body
            }
            libsail::format::Format::Hmm => {
                let $collection = libsail::seq::p7hmm::Hmm::new($input.reader()?)?;
                let $size = $crate::size_of!(Hmm, $axis);
                let $name = $crate::output::name::hmm;
                let $write = $crate::output::write::hmm
                    as $crate::output::Writer<libsail::seq::p7hmm::HmmRecord>;
                $body
            }
        }
    };
}

pub(crate) use dispatch;

/// Run `body` against the disk-backed collection the format calls for.
macro_rules! indexed {
    ($format:expr, $path:expr, |$collection:ident| $body:expr) => {
        // `dispatch!` for `--read indexed`: the same three
        // unrelated types, addressed on disk rather than
        // held in memory
        match $format {
            libsail::format::Format::Fasta => {
                let $collection = libsail::seq::fasta::IndexedFasta::open($path)?;
                $body
            }
            libsail::format::Format::Stockholm => {
                let $collection = libsail::seq::stockholm::IndexedStockholm::open($path)?;
                $body
            }
            libsail::format::Format::Hmm => {
                let $collection = libsail::seq::p7hmm::IndexedHmm::open($path)?;
                $body
            }
        }
    };
}

pub(crate) use indexed;

/// The path behind an input, for an operation reading it through an index.
pub fn indexed_path(entry: &Input) -> Result<&Path> {
    entry
        .path()
        // resolve already refuses --read indexed on a pipe,
        // so this is unreachable rather than a second policy
        .context("--read indexed needs a file, and <stdin> has no byte offsets")
}

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

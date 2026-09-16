//! Indexing a FASTA, against a reference scan of the same file.
//!
//! ```text
//! index reference <path>   read the file, counting '\n' with memchr
//! index run       <path>   build_index over the file
//! ```
//!
//! The reference counts newlines with `memchr` rather than
//! `iter().filter(|&&b| b == b'\n').count()`, which is ~50% slower over this
//! fixture for a reason unrelated to scanning: `count()` returns `usize`, so
//! the reduction widens 8→16→32→64 before it can accumulate. It is also the
//! harder comparison, since the fixture holds a newline every 52 bytes and a
//! `>` every 384.

use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;
use std::process::ExitCode;

// the reference is a memchr scan, not a stray dependency:
// referencing libsail's memchr search against a hand-rolled
// loop would price the search technology instead
use memchr::memchr_iter;

use libsail::index::build_index;
use libsail::seq::fasta::FastaDelimiter;

/// Bytes per read.
//
// mirrors the private CHUNK in libsail's src/index.rs, so
// the reference issues the same reads as the scan it bounds
const CHUNK: usize = 64 * 1024;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [sub, path] = args.as_slice() else {
        eprintln!("usage: index <reference|run> <path>");
        return ExitCode::FAILURE;
    };

    let path = Path::new(path);

    // std's boxed error rather than anyhow: this binary's
    // whole cost is measured, so its deps stay at libsail
    let line: Result<String, Box<dyn std::error::Error>> = match sub.as_str() {
        "reference" => count_newlines(path)
            .map(|n| format!("lines={n}"))
            .map_err(Into::into),
        "run" => index(path)
            .map(|n| format!("records={n}"))
            .map_err(Into::into),
        other => {
            eprintln!("unknown subcommand {other:?}");
            return ExitCode::FAILURE;
        }
    };

    match line {
        Ok(line) => {
            // observed from outside the process, so the
            // optimiser cannot skip the work
            let _ = writeln!(std::io::stdout().lock(), "{line}");
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("{sub} {}: {e}", path.display());
            ExitCode::FAILURE
        }
    }
}

/// Every newline in the file, found the way libsail finds a delimiter.
fn count_newlines(path: &Path) -> std::io::Result<u64> {
    let mut file = File::open(path)?;
    let mut buf = vec![0u8; CHUNK];
    let mut lines = 0u64;

    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            return Ok(lines);
        }

        lines += memchr_iter(b'\n', &buf[..n]).count() as u64;
    }
}

/// The thing being measured.
fn index(path: &Path) -> libsail::Result<usize> {
    let file = File::open(path)?;

    Ok(build_index::<_, FastaDelimiter>(file)?.len())
}

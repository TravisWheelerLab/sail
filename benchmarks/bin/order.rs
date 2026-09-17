//! Reading a random sample scattered against reading it forward.
//!
//! ```text
//! order reference <path>   IndexedFasta, sample(100_000), sum the lengths
//! order run       <path>   IndexedFasta, sample_in_order(100_000), sum the lengths
//! ```
//!
//! Same collection, same backend, same count, same seed: the only difference is
//! the order the `pread`s go out in. `sample` addresses the file through a
//! permutation, so consecutive reads land independently across it;
//! `sample_in_order` names its positions ascending, so the kernel's readahead
//! has something to read ahead into.
//!
//! **note: the driver's eviction is not enough on its own here.
//!         `IndexedFasta::open` scans the whole file to build the offset
//!         table, which pulls the gigabyte straight back into the page cache,
//!         and both variants then read it warm and time the same. So this
//!         evicts again itself, after the index is built and before the first
//!         sampled read. Both variants pay the same scan and the same
//!         eviction, so the difference between them is the read order alone.

use std::fs::File;
use std::io::Write;
use std::os::fd::AsRawFd;
use std::path::Path;
use std::process::ExitCode;

use libsail::collection::{Indexable, Iterable};
use libsail::seq::fasta::IndexedFasta;

/// Records drawn, the same count the `sample` benchmark uses so the two boards
/// can be read against each other.
const SAMPLE: usize = 100_000;

/// The seed both variants draw under, so they read the same records and only
/// the order differs.
const SEED: u64 = 0x5A11;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [sub, path] = args.as_slice() else {
        eprintln!("usage: order <reference|run> <path>");
        return ExitCode::FAILURE;
    };

    let path = Path::new(path);

    let line: Result<String, Box<dyn std::error::Error>> = IndexedFasta::open(path)
        .map_err(Into::into)
        .and_then(|c| {
            // between building the index and reading through
            // it: what is being timed is the reads
            evict(path)?;

            match sub.as_str() {
                "reference" => scattered(&c),
                "run" => forward(&c),
                other => Err(format!("unknown subcommand {other:?}").into()),
            }
        })
        .map(|n| format!("residues={n}"));

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

/// Drop the file's pages, so the reads that follow reach the disk.
fn evict(path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open(path)?;

    // SAFETY: the fd is open for the length of the call, and
    //         posix_fadvise reads no memory through it --
    //         0, 0 names the whole file
    let rc = unsafe { libc::posix_fadvise(file.as_raw_fd(), 0, 0, libc::POSIX_FADV_DONTNEED) };

    if rc != 0 {
        return Err(format!("posix_fadvise failed with {rc}").into());
    }

    Ok(())
}

/// Total residues over a draw read in the permutation's order.
fn scattered(c: &IndexedFasta) -> Result<u64, Box<dyn std::error::Error>> {
    let total = enough(c)?;

    // sample only chose positions, so every pread happens in
    // this sum
    let residues: u64 = c.sample(SAMPLE, SEED).iter().map(|r| r.len() as u64).sum();

    check(residues, total)
}

/// Total residues over the same draw read in the file's order.
fn forward(c: &IndexedFasta) -> Result<u64, Box<dyn std::error::Error>> {
    let total = enough(c)?;

    let residues: u64 = c
        .sample_in_order(SAMPLE, SEED)
        .map(|r| r.len() as u64)
        .sum();

    check(residues, total)
}

/// The record count, once it is known to be enough to draw from.
//
// both combinators panic when asked for more than there is,
// and a panic mid-measurement is the worse diagnostic
fn enough(c: &IndexedFasta) -> Result<usize, Box<dyn std::error::Error>> {
    let total = c.len();

    if total < SAMPLE {
        return Err(format!("{total} records, need {SAMPLE} to sample").into());
    }

    Ok(total)
}

/// The residue total, once it is known the reads came back with something.
//
// a draw that read nothing would otherwise land in the
// table as a very fast run
fn check(residues: u64, total: usize) -> Result<u64, Box<dyn std::error::Error>> {
    if residues < SAMPLE as u64 {
        return Err(format!("{SAMPLE} of {total} records summed to {residues} residues").into());
    }

    Ok(residues)
}

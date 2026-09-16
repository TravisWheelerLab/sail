//! A memory-mapped indexed reader against the `pread` one libsail ships.
//!
//! ```text
//! mmap reference <path>   IndexedFasta::open, sample(100_000), sum the lengths
//! mmap run       <path>   the same over an mmap, indexed and parsed in place
//! ```
//!
//! **note: this is an A/B, not a floor. Every other benchmark references
//!         against a task libsail should not beat; here the reference *is*
//!         libsail's current backend and the question is whether the run beats
//!         it. `sample` already sets the precedent, its reference losing 7.5x.
//!
//! Both rows are libsail's own, so this measures what ships rather than a
//! reimplementation of it. It began as a spike with a local `Indexable` over an
//! `Mmap`, which is what decided the library should grow `source::MmapSource`;
//! the shipped type replaced it and reproduced the number.
//!
//! **note: `max rss` does not mean the same thing in the two rows. Mapped pages
//!         count as resident once touched, so the run measures how much of the
//!         file was paged in -- reclaimable page cache the kernel charges to
//!         the process, not anonymous bytes an allocator holds. The reference
//!         puts the same file in the page cache and is not billed for
//!         it. Do not read the column as a like-for-like.
//!
//! The run's resident set is driven by **page granularity, not bytes read**: a
//! ~400 byte record costs a whole page, and 100,000 records drawn at random
//! over a 1 GiB file touch about 78% of its pages. Indexing is not the cause —
//! indexing through the `File` instead of the map was measured and moved it
//! only 1099 -> 890 MB, while costing wall time.

use std::borrow::Borrow;
use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

use libsail::collection::{Indexable, Iterable};
use libsail::seq::fasta::{FastaRecord, IndexedFasta, MmapFasta};

/// Records drawn, ten times what `sample` draws.
//
// sized so the gets dominate, the way the fixture is
// sized so startup does not. this benchmark is about the
// cost of a get, and both rows index the same 1 GiB
// first -- at sample's 100_000 that index is half the
// row, and the difference reads as 14% rather than the
// 1.66x per get it is
const SAMPLE: usize = 1_000_000;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [sub, path] = args.as_slice() else {
        eprintln!("usage: mmap <reference|run> <path>");
        return ExitCode::FAILURE;
    };

    let path = Path::new(path);

    // std's boxed error rather than anyhow: this binary's
    // whole cost is measured, so its deps stay at libsail
    // and the mmap crate the run is about
    let line: Result<String, Box<dyn std::error::Error>> = match sub.as_str() {
        "reference" => IndexedFasta::open(path)
            .map_err(Into::into)
            .and_then(draw)
            .map(|n| format!("records={n}")),
        // SAFETY: the fixture is written by `synth`, read by
        // this process alone, and removed only once the run
        // has exited. a benchmark is one of the few settings
        // where that contract is plainly met
        "run" => unsafe { MmapFasta::map(path) }
            .map_err(Into::into)
            .and_then(draw)
            .map(|n| format!("records={n}")),
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

/// Sample [`SAMPLE`] records, walk them, and report how many the collection
/// held. Copied from `sample.rs`, since a measured binary is self-contained.
fn draw<C>(c: C) -> Result<usize, Box<dyn std::error::Error>>
where
    C: Indexable<Record = FastaRecord>,
{
    // the collection's count is the only number the
    // generator knows, so the only one the gate can check.
    // the sample's own size and residue total are checked
    // below instead, since a draw of nothing would read as
    // a very fast run
    let total = c.len();

    // sample panics when asked for more than there is, and a
    // panic mid-measurement is the worse diagnostic
    if total < SAMPLE {
        return Err(format!("{total} records, need {SAMPLE} to sample").into());
    }

    let drawn = c.sample(SAMPLE);

    if drawn.len() != SAMPLE {
        return Err(format!("sampled {} records, expected {SAMPLE}", drawn.len()).into());
    }

    // the line that does the work: sample only chose
    // indices, so every get happens here. borrow() because
    // C is generic: Fasta lends a record and IndexedFasta
    // owns one, and this has to read the same for both
    let residues: u64 = drawn.iter().map(|r| r.borrow().len() as u64).sum();

    // every record has at least one residue, so a sum below
    // the count means the gets came back empty
    if residues < SAMPLE as u64 {
        return Err(format!("{SAMPLE} records summed to {residues} residues").into());
    }

    Ok(total)
}

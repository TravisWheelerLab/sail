//! Drawing a random sample of records, on disk against in memory.
//!
//! ```text
//! sample reference <path>   Fasta::open, sample(100_000), sum the lengths
//! sample run       <path>   IndexedFasta::open, sample(100_000), sum the lengths
//! ```
//!
//! Only the backend differs: the draw is `n` calls to `get`, a clone out of
//! memory on `Fasta` and a `pread` and a parse on `IndexedFasta`. `sample`
//! returns a `View` without touching a record, so summing the residue
//! lengths is what makes the `get`s happen — the cheapest walk that has to see
//! every sampled record.
//!
//! **note: the reference loses here, badly. `Fasta::open` parses every record
//!         before it can return one, where `IndexedFasta` parses only those
//!         asked for, so a sparse sample is faster and lighter on disk. See
//!         CLAUDE.md for where the crossover sits.

use std::borrow::Borrow;
use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

use libsail::collection::{Indexable, Iterable};
use libsail::seq::fasta::{Fasta, FastaRecord, IndexedFasta};

/// Records drawn, enough that the `get`s dominate `sample`'s pass over the
/// indices.
const SAMPLE: usize = 100_000;

/// The seed every run draws under, so two runs read the same records and the
/// numbers are comparable.
const SEED: u64 = 0x5A11;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [sub, path] = args.as_slice() else {
        eprintln!("usage: sample <reference|run> <path>");
        return ExitCode::FAILURE;
    };

    let path = Path::new(path);

    // std's boxed error rather than anyhow: this binary's
    // whole cost is measured, so its deps stay at libsail
    let line: Result<String, Box<dyn std::error::Error>> = match sub.as_str() {
        "reference" => Fasta::open(path)
            .map_err(Into::into)
            .and_then(draw)
            .map(|n| format!("records={n}")),
        "run" => IndexedFasta::open(path)
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

/// Sample [`SAMPLE`] records, walk them, and report how many the collection held.
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

    let drawn = c.sample(SAMPLE, SEED);

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

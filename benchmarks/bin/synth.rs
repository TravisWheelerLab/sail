//! Writes the FASTA every benchmark runs on, the counts they are checked
//! against, and the name list `compare fetch` asks for.
//!
//! ```text
//! synth <path>   write <path>, <path>.counts and <path>.names
//! ```
//!
//! Its own binary so no measured one links a random number generator into a
//! figure meant to describe libsail.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::process::ExitCode;

use rand::{RngExt, SeedableRng, rngs::StdRng};

/// Residues per line, matching what most FASTA writers emit.
const LINE_WIDTH: usize = 60;

/// The only size the suite measures.
//
// every measurement is a whole process, so ~1.2 ms of
// startup lands inside it. at 1 GiB that is under 1% of
// the cheapest row; where the work took microseconds it
// would be all of it
const SIZE: usize = 1 << 30;

/// Fixed, so every run measures the same bytes.
const SEED: u64 = 1;

// log-normal record lengths, fitted to real proteomes:
//
//     median = exp(MU)             = 300 aa
//     mean   = exp(MU + SIGMA^2/2) ~ 383 aa
//     p99                          ~ 1530 aa
//
// **note: the distribution is load-bearing. one fixed
//         length puts a boundary at a fixed stride, which
//         lets the branch predictor learn the
//         `>`-after-newline branch and flatters the scan.
//         a long right tail gives it no period to learn
const MU: f64 = 5.704; // ln(300)
const SIGMA: f64 = 0.7;

const MIN_LEN: usize = 10;

/// Length of titin, clamping a freak draw from the tail to something a real
/// proteome could hold.
const MAX_LEN: usize = 35_000;

const AMINO: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// Names written to `<path>.names`, for the operations that work from a list
/// of them rather than from the whole file.
//
// the same scale as the sample benchmark draws, so a fetch
// is dominated by the retrieval rather than by startup
const FETCH_NAMES: usize = 100_000;

/// Seeds the name draw.
//
// its own generator, and its own seed: drawing the names
// from the rng that picks record lengths would shift that
// stream and change the fixture's bytes
const NAMES_SEED: u64 = 2;

/// Size of the pre-drawn residue pool: a whole number of [`LINE_WIDTH`] lines.
const POOL: usize = 64 * 1024 / LINE_WIDTH * LINE_WIDTH;

/// What the generator wrote.
struct Counts {
    bytes: u64,
    records: u64,
    lines: u64,
    residues: u64,
    fetch_names: u64,
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [path] = args.as_slice() else {
        eprintln!("usage: synth <path>");
        return ExitCode::FAILURE;
    };

    let path = Path::new(path);

    match synth(path)
        .and_then(|mut c| {
            c.fetch_names = write_names(path, c.records)?;
            Ok(c)
        })
        .and_then(|c| write_counts(path, &c))
    {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("synth {}: {e}", path.display());
            ExitCode::FAILURE
        }
    }
}

/// Write a FASTA file of roughly [`SIZE`] bytes to `path`.
fn synth(path: &Path) -> std::io::Result<Counts> {
    let mut rng = StdRng::seed_from_u64(SEED);

    // why: an rng call per byte would cost more than every
    //      measurement taken over the fixture, and buys
    //      nothing -- the scan looks for "\n>" and the parser
    //      filters whitespace, so which residue a byte holds
    //      never changes a branch. only the lengths matter
    let pool: Vec<u8> = (0..POOL)
        .map(|_| AMINO[rng.random_range(0..AMINO.len())])
        .collect();
    let mut at = 0usize;

    // created fresh every run rather than cached, so a
    // change here cannot leave a stale fixture to be measured
    let mut out = BufWriter::with_capacity(1 << 20, File::create(path)?);
    let (mut bytes, mut records, mut lines, mut residues) = (0u64, 0u64, 0u64, 0u64);

    while (bytes as usize) < SIZE {
        let header = format!(">bench{records}\n");
        out.write_all(header.as_bytes())?;
        bytes += header.len() as u64;
        records += 1;
        lines += 1;

        let mut left = lognormal(&mut rng);
        while left > 0 {
            let line = LINE_WIDTH.min(left);

            // the pool is orders of magnitude longer than a
            // line, so this is one comparison per line
            if at + line > pool.len() {
                at = 0;
            }
            out.write_all(&pool[at..at + line])?;
            out.write_all(b"\n")?;
            at += line;

            bytes += line as u64 + 1;
            residues += line as u64;
            lines += 1;
            left -= line;
        }
    }

    out.flush()?;

    Ok(Counts {
        bytes,
        records,
        lines,
        residues,
        // filled in by write_names, which needs the record
        // count this function is still computing
        fetch_names: 0,
    })
}

/// One draw from the log-normal.
fn lognormal(rng: &mut StdRng) -> usize {
    // ln(0) is -inf, so a u1 of exactly zero has to be nudged
    let u1: f64 = rng.random::<f64>().max(f64::MIN_POSITIVE);
    let u2: f64 = rng.random::<f64>();

    // Box-Muller, so no rand_distr dependency:
    //     z   = sqrt(-2 ln u1) * cos(2 pi u2)
    //     len = exp(MU + SIGMA * z)
    let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();

    ((MU + SIGMA * z).exp() as usize).clamp(MIN_LEN, MAX_LEN)
}

/// Write [`FETCH_NAMES`] distinct record names to `<path>.names`, and report
/// how many were written.
//
// distinct, because a repeated name is fetched twice by a
// tool that honours the order asked and once by a tool that
// emits in file order, which would leave the two writing
// different numbers of records
fn write_names(path: &Path, records: u64) -> std::io::Result<u64> {
    let wanted = FETCH_NAMES.min(records as usize);
    let mut rng = StdRng::seed_from_u64(NAMES_SEED);
    let mut seen = HashSet::with_capacity(wanted);
    let mut drawn = Vec::with_capacity(wanted);

    // with replacement and then topped up, rather than a
    // shuffle of every index: at 100,000 draws over 2.6M
    // records the collision rate is under 2%, and a shuffle
    // would allocate the whole range to throw most of it away
    while drawn.len() < wanted {
        let n = rng.random_range(0..records);
        if seen.insert(n) {
            drawn.push(n);
        }
    }

    // in the order drawn, not sorted: a fetch of a sorted
    // list is a sequential read, which is not the access
    // pattern a name list describes
    let mut out = BufWriter::with_capacity(1 << 16, File::create(names_path(path))?);
    for n in &drawn {
        writeln!(out, "bench{n}")?;
    }
    out.flush()?;

    Ok(drawn.len() as u64)
}

/// Write what was generated to `<path>.counts`.
fn write_counts(path: &Path, c: &Counts) -> std::io::Result<()> {
    // generation happens here and not in the driver, so the
    // counts have to cross a process boundary for the board
    // header to say what was measured. a file beside the
    // fixture also outlives the run that wrote it
    let mut out = File::create(counts_path(path))?;

    // key=value, which is what Fixture::counts reads back
    writeln!(out, "bytes={}", c.bytes)?;
    writeln!(out, "records={}", c.records)?;
    writeln!(out, "lines={}", c.lines)?;
    writeln!(out, "residues={}", c.residues)?;
    // named for what asks for it rather than "names", which
    // would read as the number of named records -- every
    // record here has a name
    writeln!(out, "fetch_names={}", c.fetch_names)
}

/// `bench.fa` -> `bench.fa.counts`, an extension rather than a sibling name so
/// the two sort together.
fn counts_path(path: &Path) -> std::path::PathBuf {
    beside(path, ".counts")
}

/// `bench.fa` -> `bench.fa.names`.
fn names_path(path: &Path) -> std::path::PathBuf {
    beside(path, ".names")
}

fn beside(path: &Path, extension: &str) -> std::path::PathBuf {
    let mut name = path.as_os_str().to_os_string();
    name.push(extension);
    name.into()
}

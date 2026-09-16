//! Two halves, one shape.
//!
//! [`api`] measures libsail against the alternative it was chosen over.
//! [`cli`] measures `sail` against the tools people already reach for. Both
//! draw the same [`board`], take the same flags, and share the one generated
//! fixture.
//!
//! ```text
//! cargo run --release -- api            every libsail benchmark
//! cargo run --release -- api mmap       one of them
//! cargo run --release -- cli grep       one comparison
//! cargo run --release -- cli --repeat 9
//! ```
//!
//! Every measurement is a whole process, because that is the only thing
//! `getrusage` can be asked about — michi's `Closure` measures wall clock only,
//! since "a thread has no `wait4` to ask". So peak memory comes for free, and
//! micro-benchmarks do not: process startup lands inside every row, which is
//! also why the variants of one benchmark belong to the same binary.
//!
//! The fixture is a gigabyte, generated once into `runs/` and kept. Nothing
//! else is written: measured runs send their stdout to `/dev/null`, and the
//! only files a pass leaves are the tables beside it.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::Context;
use clap::{Args, Parser, Subcommand};
use michi::{Cmd, PipelineBuilder, Step};

mod api;
mod board;
mod cli;

// ---

#[derive(Parser)]
#[command(about = "libsail benchmarks")]
struct Cli {
    #[command(subcommand)]
    bench: Bench,
}

#[derive(Subcommand)]
enum Bench {
    /// libsail's own choices, each against the alternative it was chosen over.
    Api(Which),
    /// sail against the tools people already reach for.
    Cli(Which),
}

#[derive(Args, Clone)]
struct Which {
    /// Which one to run. Omit for all of them.
    which: Option<String>,

    #[command(flatten)]
    opts: Opts,
}

#[derive(Args, Clone)]
pub(crate) struct Opts {
    /// Runs of each measurement. michi runs a command once, so repetition is the
    /// driver's job; the report takes the minimum.
    #[arg(long, default_value_t = 5)]
    repeat: usize,

    /// Where the tables, the captured stderr, and the fixture go. Not
    /// `$TMPDIR`, which the OS is free to purge.
    #[arg(long)]
    out: Option<PathBuf>,

    /// Generate the fixture again even if one is already there.
    #[arg(long)]
    pub(crate) fresh: bool,
}

impl Opts {
    pub(crate) fn out(&self) -> PathBuf {
        self.out
            .clone()
            .unwrap_or_else(|| Path::new(env!("CARGO_MANIFEST_DIR")).join("runs"))
    }
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().bench {
        Bench::Api(w) => api::run(w.which.as_deref(), &w.opts),
        Bench::Cli(w) => cli::run(w.which.as_deref(), &w.opts),
    }
}

/// Remove every index a previous repeat left.
//
// without this the second run of `seqkit faidx` reads the
// index the first one wrote and is timed doing nothing. it
// also keeps the scanning operations honest, since a
// stray .ssi changes what esl-sfetch would do
pub(crate) fn reset(fasta: &Path) -> Cmd {
    Cmd::new("/bin/rm")
        .name("reset")
        .flag("-f")
        .path(Fixture::beside_path(fasta, "fai"))
        .path(Fixture::beside_path(fasta, "ssi"))
        .path(Fixture::beside_path(fasta, "saidx"))
}

/// Generate the fixture if it is not already there.
//
// once per invocation rather than once per operation: it is
// the same gigabyte every time, and writing it nine times
// was most of what a comparison did
pub(crate) fn ensure_fixture(fixture: &Fixture, out: &Path, o: &Opts) -> anyhow::Result<()> {
    if fixture.cached() && !o.fresh {
        return Ok(());
    }

    println!("generating {} ...", fixture.fasta.display());

    PipelineBuilder::new()
        .stderr_dir(out.join("stderr"))
        .step(Step::serial([build(&["synth"])]).name("build synth"))
        .step(Step::serial([synth(&fixture.fasta)]).name("fixture"))
        .build()?
        .run()?;

    Ok(())
}

// --- pipeline pieces

/// Where a built binary lands.
pub(crate) fn binary(name: &str) -> PathBuf {
    // cargo run builds the driver and not its siblings, so
    // the pipeline builds them itself and runs them by path
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("target/release")
        .join(name)
}

/// `cargo build --release` for the binaries a run needs.
pub(crate) fn build(bins: &[&str]) -> Cmd {
    let mut cmd = Cmd::new(option_env!("CARGO").unwrap_or("cargo"))
        .name(format!("cargo build {}", bins.join(" ")))
        .dir(env!("CARGO_MANIFEST_DIR"))
        .flag("build")
        .flag("--release");

    for bin in bins {
        cmd = cmd.arg("--bin", *bin);
    }

    cmd
}

/// Generate the fixture, as a step so its cost is visible rather than hidden in
/// the driver's startup.
pub(crate) fn synth(fasta: &Path) -> Cmd {
    Cmd::new(binary("synth")).name("synth").path(fasta)
}

// --- the fixture

/// The generated FASTA, kept between runs.
//
// generating a gigabyte per operation was most of what a
// comparison wrote. it is the same file every time, so the
// first run makes it and the rest find it
pub(crate) struct Fixture {
    pub(crate) fasta: PathBuf,
}

impl Fixture {
    pub(crate) fn at(out: &Path) -> anyhow::Result<Fixture> {
        std::fs::create_dir_all(out).with_context(|| format!("creating {}", out.display()))?;

        Ok(Fixture {
            // not named for its size: the size lives in synth,
            // and a 100mb.fa holding a gigabyte is the drift
            // a filename invites
            fasta: out.join("bench.fa"),
        })
    }

    /// Whether the fixture is already there, in which case it is taken as it
    /// is. A damaged one is deleted and made again, which is a person's job.
    pub(crate) fn cached(&self) -> bool {
        self.fasta.is_file()
    }

    fn counts_path(&self) -> PathBuf {
        self.beside("counts")
    }

    fn beside(&self, extension: &str) -> PathBuf {
        Fixture::beside_path(&self.fasta, extension)
    }

    /// `path` with `extension` appended, which is how every index beside the
    /// fixture is named.
    pub(crate) fn beside_path(path: &Path, extension: &str) -> PathBuf {
        let mut name = path.to_path_buf().into_os_string();
        name.push(".");
        name.push(extension);
        name.into()
    }

    /// What the generator says it wrote, as the same `key=value` lines the
    /// measured binaries print.
    pub(crate) fn counts(&self) -> anyhow::Result<HashMap<String, u64>> {
        let at = self.counts_path();
        let text =
            std::fs::read_to_string(&at).with_context(|| format!("reading {}", at.display()))?;

        text.lines()
            .map(|line| {
                let (key, value) = line
                    .split_once('=')
                    .with_context(|| format!("{}: {line:?} is not key=value", at.display()))?;
                Ok((key.to_string(), value.trim().parse()?))
            })
            .collect()
    }
}

impl Drop for Fixture {
    // a guard rather than a call at the end of the happy
    // path, because the paths that skip the end are the
    // likely ones: the pipeline failing, or the gate
    // rejecting a count
    //
    // **note: a signal still leaks a gigabyte -- Drop does
    //         not run on SIGINT -- but the next run
    //         overwrites it
    fn drop(&mut self) {
        // the fixture stays; these do not. compare's own
        // reset clears them before each repeat, which leaves
        // the last one standing
        for ext in ["fai", "ssi", "saidx"] {
            let mut beside = self.fasta.clone().into_os_string();
            beside.push(".");
            beside.push(ext);
            std::fs::remove_file(PathBuf::from(beside)).ok();
        }
    }
}

// --- formatting

pub(crate) fn ms(s: f64) -> String {
    format!("{:.2} ms", s * 1000.0)
}

/// Group digits so 104857705 reads as 104,857,705.
pub(crate) fn commas(n: u64) -> String {
    let d = n.to_string();
    let mut out = String::new();

    for (i, c) in d.chars().enumerate() {
        if i > 0 && (d.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }

    out
}

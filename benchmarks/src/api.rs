//! libsail's own choices, each measured against the alternative it was chosen
//! over.
//!
//! Every benchmark here is two ways of doing one thing, and the first is what
//! the second has to beat: a floor built from `memchr` for indexing, the
//! in-memory collection for sampling, `pread` for reading a mapped file. A
//! number here means something only against the row above it, which is why the
//! reference is listed first and every ratio is against it.
//!
//! The binaries doing the work live in `bin/`, one per benchmark, each taking
//! `<variant> <path>`. They print one `key=value` line so the optimiser cannot
//! throw the work away; nothing reads it.

use std::path::Path;

use anyhow::bail;
use michi::{Cmd, PipelineBuilder, Step, Table};

use crate::board::Board;
use crate::{Fixture, Opts, commas, ensure_fixture, reset};

// ---

/// One way of doing the benchmark's work.
struct Variant {
    /// The row label.
    name: &'static str,

    /// What the binary is asked for, which is the same two words everywhere.
    sub: &'static str,
}

/// One question, and the two answers being timed against each other.
struct Bench {
    /// The row group, and the binary under `bin/` that runs it.
    name: &'static str,

    /// What both variants do, printed under the board. A reader who cannot
    /// check the claim cannot read the numbers.
    same_work: &'static str,

    /// The reference first, because every ratio is against it.
    variants: &'static [Variant],
}

fn benches() -> Vec<Bench> {
    vec![
        Bench {
            name: "index",
            same_work: "cross a 1 GiB FASTA and count what it holds",
            variants: &[
                Variant {
                    name: "memchr floor",
                    sub: "reference",
                },
                Variant {
                    name: "build_index",
                    sub: "run",
                },
            ],
        },
        Bench {
            name: "sample",
            same_work: "draw 100,000 records and total their residues",
            variants: &[
                Variant {
                    name: "in memory",
                    sub: "reference",
                },
                Variant {
                    name: "on disk",
                    sub: "run",
                },
            ],
        },
        Bench {
            name: "mmap",
            same_work: "read 1,000,000 records out of an indexed collection",
            variants: &[
                Variant {
                    name: "pread",
                    sub: "reference",
                },
                Variant {
                    name: "mmap",
                    sub: "run",
                },
            ],
        },
    ]
}

// --- the run

pub(crate) fn run(which: Option<&str>, o: &Opts) -> anyhow::Result<()> {
    let all = benches();
    let chosen: Vec<&Bench> = match which {
        None => all.iter().collect(),
        Some(name) => {
            let found: Vec<&Bench> = all.iter().filter(|b| b.name == name).collect();
            if found.is_empty() {
                let known: Vec<&str> = all.iter().map(|b| b.name).collect();
                bail!(
                    "no benchmark called {name:?}; there is {}",
                    known.join(", ")
                )
            }

            found
        }
    };

    let out = o.out();
    let fixture = Fixture::at(&out)?;
    ensure_fixture(&fixture, &out, o)?;

    let counts = fixture.counts()?;
    let get = |k: &str| counts.get(k).copied().unwrap_or(0);
    let header = format!(
        "bench.fa - {} bytes, {} records, {} residues - min of {}",
        commas(get("bytes")),
        commas(get("records")),
        commas(get("residues")),
        o.repeat
    );

    let plan: Vec<(&'static str, &'static str, Vec<&'static str>)> = chosen
        .iter()
        .map(|b| {
            (
                b.name,
                b.same_work,
                b.variants.iter().map(|v| v.name).collect(),
            )
        })
        .collect();

    let board = Board::new(&plan, header, "base", o.repeat);
    let painter = board.start();

    let ran = chosen
        .iter()
        .try_for_each(|bench| one(bench, &fixture, &board, o));

    // the thread lets go of stdout before anything else
    // writes to it, including the error that stopped the run
    painter.stop();
    ran?;

    board.finish();

    Ok(())
}

fn one(bench: &Bench, fixture: &Fixture, board: &Board, o: &Opts) -> anyhow::Result<()> {
    let out = o.out();

    let mut pipeline = PipelineBuilder::new()
        .stderr_dir(out.join("stderr"))
        .step(Step::serial([crate::build(&[bench.name])]).name("build"))
        .step(warmup(bench, &fixture.fasta));

    for variant in bench.variants {
        pipeline =
            pipeline.step(Step::serial(runs(bench, variant, &fixture.fasta, o)).name(variant.name));
    }

    pipeline
        .sink(Table::new(out.join(format!("api-{}.tbl", bench.name))))
        .sink(board.clone())
        .build()?
        .run()?;

    Ok(())
}

/// One unrecorded run of each variant, so a measured one is not the first to
/// touch the file.
fn warmup(bench: &Bench, fasta: &Path) -> Step {
    let mut cmds = Vec::new();
    for variant in bench.variants {
        cmds.push(reset(fasta));
        cmds.push(cmd(bench, variant, fasta).name(format!("warm {}", variant.name)));
    }

    Step::serial(cmds).name("warmup")
}

/// `repeat` runs of one variant, each from a clean slate.
fn runs(bench: &Bench, variant: &Variant, fasta: &Path, o: &Opts) -> Vec<Cmd> {
    let mut cmds = Vec::new();

    for i in 1..=o.repeat {
        // an index the run before it left behind would
        // otherwise be there to be found, which is a
        // different measurement
        cmds.push(reset(fasta));
        cmds.push(
            cmd(bench, variant, fasta)
                .name(format!("{}#{i}", variant.sub))
                .field("op", bench.name)
                .field("bench", variant.name),
        );
    }

    cmds
}

fn cmd(bench: &Bench, variant: &Variant, fasta: &Path) -> Cmd {
    Cmd::new(crate::binary(bench.name))
        .sub(variant.sub)
        .path(fasta)
}

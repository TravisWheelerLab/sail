//! One operation, several tools, one row each.
//!
//! The sibling of [`api`](crate::api), and deliberately not the same thing.
//! There the first row is a floor libsail should beat. Here every row is a
//! peer, and losing to one is a result rather than a bug. Nothing is
//! subtracted from anything in either.
//!
//! ```text
//! cargo run --release -- cli stats
//! cargo run --release -- cli --repeat 9
//! ```
//!
//! **The hard part is not the plumbing, it is what counts as the same work.**
//! Each [`Op`] carries the argument for why its tools are comparable, and an
//! operation that cannot make that argument does not belong here. Three have
//! been turned away on that test so far: `sail count` against the `stats`
//! peers, which do strictly more; `esl-sfetch -f` against the `fetch` peers,
//! which scan where it needs an `.ssi` first; and `seqkit sample -n`, which
//! draws about n records where `sail sample -n` draws exactly n.
//!
//! Every measured run sends its stdout to `/dev/null`. Capturing it cost about
//! 45 GiB a pass and, worse, was measured: half of `reformat`'s wall clock was
//! writing the capture rather than doing the work.

use std::path::{Path, PathBuf};

use anyhow::bail;
use michi::{Cmd, PipelineBuilder, Step, Table};

use crate::board::Board;
use crate::{Fixture, Opts, commas, ensure_fixture, evict, reset};

// ---

/// One tool's turn at an operation.
struct Tool {
    /// The row label, and the field that groups a tool's repeats.
    name: &'static str,

    /// What is invoked. Resolved through [`tool`] before the run, so a missing
    /// tool is reported rather than measured as a fast failure.
    program: &'static str,

    /// The command, for this fixture.
    argv: fn(&Path) -> Cmd,
}

/// One comparable operation.
struct Op {
    name: &'static str,

    /// Why these tools are doing the same work. Printed above the rows,
    /// because a reader who cannot check the claim cannot read the numbers.
    same_work: &'static str,

    tools: &'static [Tool],
}

// --- the operations

fn ops() -> Vec<Op> {
    vec![
        stats(),
        index(),
        names(),
        reformat(),
        fetch(),
        grep(),
        filter(),
        sort(),
        dedup(),
        skim(),
        sample(),
    ]
}

/// Summary statistics over a whole FASTA.
fn stats() -> Op {
    Op {
        name: "stats",
        // the fair three-way: all of these walk the file and
        // report count, total, min, max and mean. `sail
        // count` is NOT in here, because counting records is
        // strictly less work than the other two are doing
        same_work: "count, total, min, max and mean residues over every record",
        tools: &[
            Tool {
                name: "sail",
                program: SAIL,
                argv: |fasta| sail().sub("stats").path(fasta),
            },
            Tool {
                name: "seqkit",
                program: "seqkit",
                argv: |fasta| tool("seqkit").sub("stats").path(fasta),
            },
            Tool {
                name: "esl-seqstat",
                program: "esl-seqstat",
                argv: |fasta| tool("esl-seqstat").path(fasta),
            },
        ],
    }
}

/// Building a record index beside the file.
fn index() -> Op {
    Op {
        name: "index",
        same_work: "scan the file and write an index of every record beside it",
        tools: &[
            Tool {
                name: "sail",
                program: SAIL,
                argv: |fasta| sail().sub("index").path(fasta),
                // "2660597 fasta records indexed to ..."
            },
            Tool {
                name: "seqkit",
                program: "seqkit",
                argv: |fasta| tool("seqkit").sub("faidx").path(fasta),
                // faidx prints progress and no count, so the
                // .fai it wrote becomes the capture and the
                // gate counts its lines -- one per record
            },
            Tool {
                name: "esl-sfetch",
                program: "esl-sfetch",
                argv: |fasta| tool("esl-sfetch").flag("--index").path(fasta),
                // "Indexed 2660597 sequences (2660597 names)."
            },
        ],
    }
}

/// Every record's name, one per line.
fn names() -> Op {
    Op {
        name: "names",
        // two rows, not three: easel has no command that
        // prints names alone. esl-seqstat -a prints a name
        // *and* a length per record, which is more work
        same_work: "print every record's name and nothing else (no easel equivalent)",
        tools: &[
            Tool {
                name: "sail",
                program: SAIL,
                argv: |fasta| sail().sub("names").path(fasta),
            },
            Tool {
                name: "seqkit",
                program: "seqkit",
                argv: |fasta| tool("seqkit").sub("seq").flag("-n").flag("-i").path(fasta),
            },
        ],
    }
}

/// Parsing every record and writing it back out.
fn reformat() -> Op {
    Op {
        name: "reformat",
        // the strongest three-way in the suite: over this
        // fixture all three write byte-identical output, so
        // "the same work" is not an argument here but a
        // checked fact. 60 columns is what sail emits and
        // what esl-reformat defaults to, and seqkit is told
        // to match
        same_work: "parse every record and write it back as 60-column FASTA (byte-identical output)",
        tools: &[
            Tool {
                name: "sail",
                program: SAIL,
                argv: |fasta| sail().sub("reformat").path(fasta),
            },
            Tool {
                name: "seqkit",
                program: "seqkit",
                argv: |fasta| tool("seqkit").sub("seq").arg("-w", LINE_WIDTH).path(fasta),
            },
            Tool {
                name: "esl-reformat",
                program: "esl-reformat",
                argv: |fasta| tool("esl-reformat").flag("fasta").path(fasta),
            },
        ],
    }
}

/// Extracting the records named in a file.
fn fetch() -> Op {
    Op {
        name: "fetch",
        // two rows, not three: esl-sfetch -f refuses to run
        // without an .ssi beside the file ("Failed to open
        // SSI index"), so it answers from an index where
        // these two scan. that is a different operation, and
        // pre-building the index for it would hide the cost
        // that makes it different
        same_work: "read a list of 100,000 names and write out those records (easel needs an index first)",
        tools: &[
            Tool {
                name: "sail",
                program: SAIL,
                argv: |fasta| {
                    sail()
                        .sub("fetch")
                        .path(fasta)
                        .arg("-f", Fixture::beside_path(fasta, "names"))
                },
            },
            Tool {
                name: "seqkit",
                program: "seqkit",
                // grep -f is seqkit's fetch-by-name-list; by
                // default it matches the ID, which is the
                // whole name in this fixture
                argv: |fasta| {
                    tool("seqkit")
                        .sub("grep")
                        .arg("-f", Fixture::beside_path(fasta, "names"))
                        .path(fasta)
                },
            },
        ],
    }
}

/// Keeping the records whose name matches a pattern.
fn grep() -> Op {
    Op {
        name: "grep",
        // two rows: easel has no name-pattern filter.
        // seqkit is given -n so it matches the full header
        // rather than the ID, which is the same string here
        // because synth writes no descriptions
        same_work: "keep every record whose name matches a regex (no easel equivalent)",
        // the share of names ending in a given digit is
        // about a tenth, but "about" is not a number the
        // gate can use: the exact count depends on how many
        // records the 1 GiB happened to hold
        tools: &[
            Tool {
                name: "sail",
                program: SAIL,
                argv: |fasta| sail().sub("grep").flag(PATTERN).path(fasta),
            },
            Tool {
                name: "seqkit",
                program: "seqkit",
                argv: |fasta| {
                    tool("seqkit")
                        .sub("grep")
                        .flag("-n")
                        .flag("-r")
                        .arg("-p", PATTERN)
                        .path(fasta)
                },
            },
        ],
    }
}

/// Keeping the records of at least a given length.
fn filter() -> Op {
    Op {
        name: "filter",
        same_work: "keep every record of at least 500 residues (no easel equivalent)",
        // as with grep, the survivors are a property of the
        // fixture's length draw rather than of anything the
        // generator wrote down
        tools: &[
            Tool {
                name: "sail",
                program: SAIL,
                argv: |fasta| sail().sub("filter").path(fasta).arg("--min", MIN_RESIDUES),
            },
            Tool {
                name: "seqkit",
                program: "seqkit",
                argv: |fasta| {
                    tool("seqkit")
                        .sub("seq")
                        .arg("-m", MIN_RESIDUES)
                        .path(fasta)
                },
            },
        ],
    }
}

/// Ordering every record by name.
fn sort() -> Op {
    Op {
        name: "sort",
        // both default to ordering by name, and both hold
        // the whole file to do it. easel has no sort
        same_work: "order every record by name and write them all back (no easel equivalent)",
        tools: &[
            Tool {
                name: "sail",
                program: SAIL,
                argv: |fasta| sail().sub("sort").path(fasta),
            },
            Tool {
                name: "seqkit",
                program: "seqkit",
                argv: |fasta| tool("seqkit").sub("sort").path(fasta),
            },
        ],
    }
}

/// Dropping records with a repeated name.
/// Drawing a fixed number of records at random, keeping the file's order.
fn skim() -> Op {
    Op {
        name: "skim",
        // **note: seqkit sample -p is NOT in here. it is
        //         bernoulli, so it draws n +/- sqrt(n) and
        //         not n, and a tool that may return 988
        //         records where another returns 1000 is
        //         doing different work. sample2 -2 is the
        //         fair peer: its own help promises an exact
        //         count with equal probability per record,
        //         which is what skim promises
        same_work: "draw exactly 10,000 records at random and write them in the file's order",
        tools: &[
            Tool {
                name: "sail",
                program: SAIL,
                argv: |fasta| {
                    sail()
                        .sub("skim")
                        .arg("-n", "10000")
                        .arg("--seed", "7")
                        .path(fasta)
                },
            },
            Tool {
                name: "seqkit",
                program: "seqkit",
                // -2 is two-pass, which is what sail does
                // without an index beside the file, and
                // reset deletes that index before every run
                argv: |fasta| {
                    tool("seqkit")
                        .sub("sample2")
                        .arg("-n", "10000")
                        .flag("-2")
                        .arg("--rand-seed", "7")
                        .path(fasta)
                },
            },
        ],
    }
}

/// Drawing a fixed number of records at random, in a random order.
fn sample() -> Op {
    Op {
        name: "sample",
        // separate from skim rather than another row on it:
        // these two write the records in a random order and
        // skim's two write them in the file's, so one board
        // would be comparing four tools doing two jobs
        same_work: "draw exactly 10,000 records at random and write them in a random order",
        tools: &[
            Tool {
                name: "sail",
                program: SAIL,
                argv: |fasta| {
                    sail()
                        .sub("sample")
                        .arg("-n", "10000")
                        .arg("--seed", "7")
                        .path(fasta)
                },
            },
            Tool {
                name: "seqkit",
                program: "seqkit",
                argv: |fasta| {
                    tool("seqkit")
                        .sub("sample2")
                        .arg("-n", "10000")
                        .arg("--rand-seed", "7")
                        .path(fasta)
                },
            },
        ],
    }
}

fn dedup() -> Op {
    Op {
        name: "dedup",
        // **note: synth names records bench0..benchN, so the
        //         fixture holds no duplicate at all. what
        //         this measures is the cost of hashing every
        //         name and emitting every record -- the path
        //         a deduplication takes when it finds
        //         nothing, which is the common case on real
        //         data. it does not measure the removal
        same_work: "hash every record's name and emit the first of each (no duplicates in this fixture)",
        tools: &[
            Tool {
                name: "sail",
                program: SAIL,
                argv: |fasta| sail().sub("dedup").path(fasta),
            },
            Tool {
                name: "seqkit",
                program: "seqkit",
                argv: |fasta| tool("seqkit").sub("rmdup").flag("-n").path(fasta),
            },
        ],
    }
}

// --- what the operations are parameterised by

/// Residues per line for [`reformat`], which is sail's default and easel's.
const LINE_WIDTH: usize = 60;

/// The name pattern [`grep`] matches.
//
// anchored at the end, so it is a regex both tools have to
// run rather than a substring search either could special-
// case. synth names records bench0..benchN, so this keeps
// roughly a tenth of them
const PATTERN: &str = "7$";

/// The length [`filter`] keeps records at or above.
//
// the fixture's lengths are log-normal with a median of
// 300, so this keeps roughly the top quarter -- enough
// output that the writing is not free, few enough that the
// filter is doing something
const MIN_RESIDUES: usize = 500;

// --- output readers

// --- the peer tools

/// Where `setup.sh` puts the tools sail is measured against.
fn tools_bin() -> &'static Path {
    // a OnceLock would buy nothing: this is two joins, and
    // it is called once per Cmd built
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/tools/bin"))
}

/// A command for `program` as `setup.sh` installed it, falling back to a bare
/// name for `PATH` to resolve.
//
// tools/bin first, so a run measures the versions setup.sh
// pinned even on a machine carrying its own seqkit at
// another one. the fallback is what lets someone who
// already has these tools skip setup.sh, at the cost of
// numbers that are theirs rather than ours
fn tool(program: &'static str) -> Cmd {
    let local = tools_bin().join(program);

    if local.is_file() {
        Cmd::new(local).name(program)
    } else {
        Cmd::new(program)
    }
}

// --- the sail binary

/// Stands in for a path only known at run time.
const SAIL: &str = "<sail>";

fn workspace() -> &'static Path {
    // benchmarks/ is deliberately outside the workspace, so
    // the CLI it measures is built through the parent
    // manifest rather than this one
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("benchmarks/ has a parent")
}

fn sail_bin() -> PathBuf {
    workspace().join("target/release/sail")
}

fn sail() -> Cmd {
    Cmd::new(sail_bin())
}

fn build_sail() -> Cmd {
    Cmd::new(option_env!("CARGO").unwrap_or("cargo"))
        .name("cargo build -p sail")
        .dir(workspace())
        .flag("build")
        .flag("--release")
        .arg("--package", "sail")
}

// --- pipeline pieces

/// One unrecorded run of each tool, so a measured one is never the first to
/// load a binary and never the first to fail.
//
// this was the page cache's warmup until eviction took that
// over, and eviction runs before every measured run rather
// than once at the top. what is left is the dynamic loader,
// and a bad invocation surfacing before anything is timed
fn dry_run(op: &Op, fasta: &Path) -> Step {
    let mut cmds = Vec::new();
    for tool in op.tools {
        // no bench field, so a dry run reaches the board no
        // more than it reaches the report
        cmds.push(reset(fasta));
        cmds.push((tool.argv)(fasta).name(format!("dry {}", tool.name)));
    }

    Step::serial(cmds).name("dry run")
}

/// `repeat` runs of one tool, each from a clean slate and each capturing its
/// stdout.
fn runs(op: &Op, tool: &Tool, fasta: &Path, o: &Opts) -> Vec<Cmd> {
    let mut cmds = Vec::new();

    for i in 1..=o.repeat {
        cmds.push(reset(fasta));

        // after reset, so the pages reset dirtied are gone
        // too, and immediately before the command, so
        // nothing between them pulls the file back in
        if !o.warm {
            cmds.push(evict(fasta));
        }

        cmds.push(
            (tool.argv)(fasta)
                .name(format!("{}#{i}", tool.name))
                .field("op", op.name)
                .field("bench", tool.name),
        );
    }

    cmds
}

/// Fail before any measurement if a tool is not installed.
//
// a missing program would otherwise exec-fail in
// milliseconds and land in the table as the fastest row
fn preflight(op: &Op) -> anyhow::Result<()> {
    let missing: Vec<&str> = op
        .tools
        .iter()
        .filter(|t| t.program != SAIL)
        .filter(|t| !tools_bin().join(t.program).is_file() && !on_path(t.program))
        .map(|t| t.program)
        .collect();

    if !missing.is_empty() {
        bail!(
            "cli {} needs {}, which {} in {} nor on PATH -- run benchmarks/setup.sh",
            op.name,
            missing.join(" and "),
            if missing.len() == 1 {
                "is neither"
            } else {
                "are neither"
            },
            tools_bin().display(),
        )
    }

    Ok(())
}

fn on_path(program: &str) -> bool {
    std::env::var_os("PATH")
        .map(|path| std::env::split_paths(&path).any(|dir| dir.join(program).is_file()))
        .unwrap_or(false)
}

// --- the run

pub(crate) fn run(which: Option<&str>, o: &Opts) -> anyhow::Result<()> {
    let all = ops();
    let chosen: Vec<&Op> = match which {
        None => all.iter().collect(),
        Some(name) => {
            let found: Vec<&Op> = all.iter().filter(|op| op.name == name).collect();
            if found.is_empty() {
                let known: Vec<&str> = all.iter().map(|op| op.name).collect();
                bail!(
                    "no comparison called {name:?}; there is {}",
                    known.join(", ")
                )
            }

            found
        }
    };

    // every tool checked before anything is measured, so a
    // missing one is reported rather than drawn as a gap
    for op in &chosen {
        preflight(op)?;
    }

    let out = o.out();
    let fixture = Fixture::at(&out)?;
    ensure_fixture(&fixture, &out, o)?;

    let counts = fixture.counts()?;
    let get = |k: &str| counts.get(k).copied().unwrap_or(0);
    let header = format!(
        "bench.fa - {} bytes, {} records, {} residues - min of {} - {}",
        commas(get("bytes")),
        commas(get("records")),
        commas(get("residues")),
        o.repeat,
        if o.warm { "warm" } else { "cold" }
    );

    let plan: Vec<(&'static str, &'static str, Vec<&'static str>)> = chosen
        .iter()
        .map(|op| {
            (
                op.name,
                op.same_work,
                op.tools.iter().map(|t| t.name).collect(),
            )
        })
        .collect();

    let board = Board::new(&plan, header, "sail", o.repeat);
    let painter = board.start();

    let ran = chosen
        .iter()
        .try_for_each(|op| one(op, &fixture, &board, o));

    // the thread lets go of stdout before anything else
    // writes to it, including the error that stopped the run
    painter.stop();
    ran?;

    board.finish();

    Ok(())
}

fn one(op: &Op, fixture: &Fixture, board: &Board, o: &Opts) -> anyhow::Result<()> {
    let out = o.out();

    let mut pipeline = PipelineBuilder::new()
        .stderr_dir(out.join("stderr"))
        .step(Step::serial([build_sail(), crate::build(&["evict"])]).name("build"))
        .step(dry_run(op, &fixture.fasta));

    for tool in op.tools {
        pipeline = pipeline.step(Step::serial(runs(op, tool, &fixture.fasta, o)).name(tool.name));
    }

    pipeline
        .sink(Table::new(out.join(format!("compare-{}.tbl", op.name))))
        .sink(board.clone())
        .build()?
        .run()?;

    Ok(())
}

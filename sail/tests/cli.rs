//! What every subcommand writes today, pinned byte for byte.
//!
//! These exist so a change of read backend cannot quietly change anyone's
//! output. Several of the commands pinned here have no unit tests at all, and
//! those are the ones such a change would alter without failing anything.
//!
//! They drive the real binary rather than the `*Args` structs, because what has
//! to stay fixed is what a caller sees: stdout, the exit code, and the first
//! line of stderr, clap's own parsing included.
//!
//! ```text
//! SAIL_BLESS=1 cargo test -p sail --test cli    # rewrite the goldens
//! ```

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

// ---

struct Run {
    stdout: Vec<u8>,
    stderr: String,
    ok: bool,
}

/// Run `sail` with `args`, feeding it `stdin`.
fn sail(args: &[&str], stdin: &[u8]) -> Run {
    use std::io::Write;

    // CARGO_BIN_EXE_sail is set for an integration test of
    // a package with a binary target, so this needs no dev
    // dependency and no guess at target/
    let mut child = Command::new(env!("CARGO_BIN_EXE_sail"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("sail is built before its integration tests run");

    child
        .stdin
        .take()
        .expect("stdin was piped")
        .write_all(stdin)
        .ok();

    let out = child.wait_with_output().expect("sail ran to completion");

    Run {
        stdout: out.stdout,
        stderr: String::from_utf8_lossy(&out.stderr).into_owned(),
        ok: out.status.success(),
    }
}

// ---

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn blessing() -> bool {
    std::env::var_os("SAIL_BLESS").is_some()
}

/// Compare `got` against the golden for `case`, or write it when blessing.
/// The report for a changed case, and `None` when it matches.
#[must_use]
fn golden(case: &str, got: &[u8]) -> Option<String> {
    let at = golden_dir().join(format!("{case}.out"));

    if blessing() {
        std::fs::create_dir_all(golden_dir()).expect("the golden directory is writable");
        std::fs::write(&at, got).expect("the golden is writable");

        return None;
    }

    let want = std::fs::read(&at).unwrap_or_else(|e| {
        panic!(
            "{}: {e}. run SAIL_BLESS=1 cargo test -p sail --test cli",
            at.display()
        )
    });

    if got == want {
        return None;
    }

    // compared as bytes, reported as text: a record name
    // need not be utf-8, so from_utf8_lossy shapes the
    // message and never the comparison
    Some(format!(
        "{case} changed\n--- want ---\n{}\n--- got ---\n{}",
        head(&want),
        head(got),
    ))
}

/// The first few lines of a golden, so a report over many cases stays readable.
fn head(bytes: &[u8]) -> String {
    let text = String::from_utf8_lossy(bytes);
    let mut out: String = text.lines().take(6).collect::<Vec<_>>().join("\n");

    if text.lines().count() > 6 {
        out.push_str("\n…");
    }

    out
}

/// Fail once, naming every case that changed.
fn report(changed: Vec<String>) {
    // one assert for the whole run, so a change touching
    // many operations names all of them rather than
    // stopping at the first
    assert!(
        changed.is_empty(),
        "{} case(s) changed:\n\n{}",
        changed.len(),
        changed.join("\n\n")
    );
}

// ---

struct Fixture {
    /// What the case is called, and the golden's stem.
    stem: &'static str,
    /// The fixture file's name, under `fixtures/`.
    file: &'static str,
    /// A record the file holds, for `fetch`.
    name: &'static str,
    /// A pattern matching at least one name, for `grep`.
    pattern: &'static str,
    /// A size at least one record clears, for `filter`.
    min: &'static str,
}

const FIXTURES: [Fixture; 3] = [
    Fixture {
        stem: "proteins",
        file: "proteins.fa",
        name: "DLG4_HUMAN",
        pattern: "HUMAN",
        min: "400",
    },
    Fixture {
        stem: "families",
        file: "families.sto",
        name: "PDZ",
        pattern: "PDZ",
        min: "6",
    },
    Fixture {
        stem: "models",
        file: "models.hmm",
        name: "PDZ",
        pattern: "PDZ",
        min: "46",
    },
];

fn fixture(name: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../fixtures")
        .join(name)
        .display()
        .to_string()
}

/// Every invocation pinned, as (case name, argv).
fn cases(f: &Fixture) -> Vec<(String, Vec<String>)> {
    let path = fixture(f.file);
    let case = |op: &str| format!("{op}-{}", f.stem);
    let argv = |args: &[&str]| args.iter().map(|a| a.to_string()).collect::<Vec<_>>();

    // split and index are not here: both write files rather
    // than a stream, so they get their own tests below
    vec![
        (case("cat"), argv(&["cat", &path])),
        (case("cat-twice"), argv(&["cat", &path, &path])),
        (case("count"), argv(&["count", &path])),
        (case("count-twice"), argv(&["count", &path, &path])),
        (case("dedup"), argv(&["dedup", &path])),
        (case("fetch"), argv(&["fetch", &path, f.name])),
        (case("filter"), argv(&["filter", &path, "--min", f.min])),
        (case("get"), argv(&["get", &path, "1"])),
        (case("grep"), argv(&["grep", f.pattern, &path])),
        (case("grep-invert"), argv(&["grep", f.pattern, &path, "-v"])),
        (case("head"), argv(&["head", &path, "-n", "1"])),
        (case("names"), argv(&["names", &path])),
        (case("reformat"), argv(&["reformat", &path])),
        (case("rename"), argv(&["rename", &path, "--prefix", "x_"])),
        (
            case("sample"),
            argv(&["sample", &path, "-n", "1", "--seed", "7"]),
        ),
        (case("shuffle"), argv(&["shuffle", &path, "--seed", "7"])),
        (case("sort"), argv(&["sort", &path])),
        (case("sort-size"), argv(&["sort", &path, "--key", "size"])),
        (case("stats"), argv(&["stats", &path])),
        (case("tail"), argv(&["tail", &path, "-n", "1"])),
    ]
}

// ---

#[test]
fn every_command_writes_what_it_wrote_before() {
    let mut changed = Vec::new();

    for f in &FIXTURES {
        for (case, args) in cases(f) {
            let args: Vec<&str> = args.iter().map(String::as_str).collect();
            let run = sail(&args, b"");

            assert!(run.ok, "{case} failed: {}", run.stderr);
            changed.extend(golden(&case, &run.stdout));
        }
    }

    report(changed);
}

#[test]
fn a_file_on_stdin_reads_the_same_as_the_file_itself() {
    // identifying a pipe means reading its head and putting
    // it back, which is where a streaming reader is most
    // likely to break
    for f in &FIXTURES {
        let bytes = std::fs::read(fixture(f.file)).unwrap();

        for op in ["cat", "count", "names", "stats", "dedup"] {
            let from_file = sail(&[op, &fixture(f.file)], b"");
            let from_pipe = sail(&[op, "-"], &bytes);

            assert!(from_pipe.ok, "{op} on stdin failed: {}", from_pipe.stderr);
            assert_eq!(
                from_pipe.stdout, from_file.stdout,
                "{op} {} differs between a path and a pipe",
                f.stem
            );
        }
    }
}

#[test]
fn validate_accepts_every_fixture() {
    for f in &FIXTURES {
        let run = sail(&["validate", &fixture(f.file)], b"");

        assert!(run.ok, "validate {} failed: {}", f.stem, run.stderr);
        assert!(
            String::from_utf8_lossy(&run.stdout).contains("ok"),
            "validate {} did not report ok",
            f.stem
        );
    }
}

#[test]
fn the_refusals_stay_refusals() {
    // each of these is a deliberate error with a sentence
    // attached, and each is a thing a backend change could
    // turn into a silent success
    let two_formats = sail(
        &["cat", &fixture("proteins.fa"), &fixture("models.hmm")],
        b"",
    );
    assert!(!two_formats.ok);
    assert!(
        two_formats
            .stderr
            .contains("one operation reads one format")
    );

    let wrong_format = sail(&["count", &fixture("proteins.fa"), "--format", "hmm"], b"");
    assert!(!wrong_format.ok);
    assert!(wrong_format.stderr.contains("--format hmm"));

    let missing_name = sail(&["fetch", &fixture("proteins.fa"), "NO_SUCH_RECORD"], b"");
    assert!(!missing_name.ok);
    assert!(missing_name.stderr.contains("NO_SUCH_RECORD"));

    let skipped = sail(
        &[
            "fetch",
            &fixture("proteins.fa"),
            "NO_SUCH_RECORD",
            "--skip-missing",
        ],
        b"",
    );
    assert!(skipped.ok, "{}", skipped.stderr);
    assert!(skipped.stdout.is_empty());

    for op in ["validate", "index"] {
        let piped = sail(&[op, "-"], b">a\nACGT\n");
        assert!(!piped.ok, "{op} accepted stdin");
    }
}

// ---

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        // the edge inputs are written per test rather than
        // into fixtures/, which is the repo's checked-in
        // input and not a scratch directory
        let at = std::env::temp_dir().join(format!("sail-cli-{tag}-{}", std::process::id()));
        std::fs::create_dir_all(&at).expect("a scratch directory");

        Scratch(at)
    }

    fn write(&self, name: &str, bytes: &[u8]) -> String {
        let at = self.0.join(name);
        std::fs::write(&at, bytes).expect("the scratch file is writable");

        at.display().to_string()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).ok();
    }
}

#[test]
fn the_awkward_inputs_keep_their_answers() {
    let s = Scratch::new("edge");
    let mut changed = Vec::new();

    let cases: [(&str, &[u8]); 6] = [
        ("empty.fa", b""),
        ("one.fa", b">only\nACGT\n"),
        // no trailing newline: the last record's extent is
        // EOF rather than the next delimiter
        ("no-final-newline.fa", b">a\nACGT\n>b\nGGTT"),
        ("crlf.fa", b">a desc\r\nACGT\r\n>b\r\nGGTT\r\n"),
        // the same name twice: dedup keeps the first, fetch
        // picks the first, sort is stable
        ("repeated.fa", b">dup\nAAAA\n>other\nCCCC\n>dup\nTTTT\n"),
        // a name that is not utf-8, which names and grep
        // both support on purpose
        ("not-utf8.fa", b">a\xff\xfeb\nACGT\n"),
    ];

    for (name, bytes) in cases {
        let path = s.write(name, bytes);
        let stem = name.trim_end_matches(".fa");

        for op in ["count", "names", "cat", "reformat", "dedup", "sort"] {
            let run = sail(&[op, &path], b"");

            // an empty file has no format to detect, so every
            // operation refuses it. this pins that answer
            // rather than endorsing it
            if stem == "empty" {
                assert!(!run.ok, "{op} accepted an empty file");
                continue;
            }

            assert!(run.ok, "{op} {stem} failed: {}", run.stderr);
            changed.extend(golden(&format!("{op}-{stem}"), &run.stdout));
        }
    }

    report(changed);
}

#[test]
fn an_alignment_with_no_id_still_has_an_answer() {
    let mut changed = Vec::new();
    let s = Scratch::new("noid");

    // Stockholm does not require #=GF ID, so a record can
    // have no name at all -- the case behind name_of
    // returning Option
    let path = s.write("no-id.sto", b"# STOCKHOLM 1.0\nseq1 ACGT\nseq2 ACGT\n//\n");

    for op in ["count", "names", "cat", "dedup", "sort"] {
        let run = sail(&[op, &path], b"");

        assert!(
            run.ok,
            "{op} on an unnamed alignment failed: {}",
            run.stderr
        );
        changed.extend(golden(&format!("{op}-no-id"), &run.stdout));
    }

    // --keep-unnamed is the flag that makes the absence
    // visible rather than skipped
    let kept = sail(&["names", &path, "--keep-unnamed"], b"");
    assert!(kept.ok, "{}", kept.stderr);
    changed.extend(golden("names-no-id-kept", &kept.stdout));

    report(changed);
}

// ---

#[test]
fn split_writes_the_parts_it_says_it_does() {
    let s = Scratch::new("split");
    let path = s.write("in.fa", b">a\nAAAA\n>b\nCCCC\n>c\nGGGG\n>d\nTTTT\n");
    let prefix = s.0.join("part").display().to_string();

    let run = sail(&["split", &path, "-n", "2", "--prefix", &prefix], b"");
    assert!(run.ok, "{}", run.stderr);

    let mut parts: Vec<PathBuf> = std::fs::read_dir(&s.0)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.file_name().unwrap().to_string_lossy().starts_with("part"))
        .collect();
    parts.sort();

    assert_eq!(parts.len(), 2, "expected two parts, got {parts:?}");

    let mut joined = Vec::new();
    for part in &parts {
        joined.extend(std::fs::read(part).unwrap());
    }
    report(golden("split-parts-joined", &joined).into_iter().collect());
}

#[test]
fn index_writes_an_index_file_that_loads_back() {
    let s = Scratch::new("index");
    let path = s.write("in.fa", b">a\nAAAA\n>b\nCCCC\n");

    let run = sail(&["index", &path], b"");
    assert!(run.ok, "{}", run.stderr);

    let index_file = format!("{path}.saidx");
    assert!(
        Path::new(&index_file).exists(),
        "index wrote no index file beside the input"
    );

    // matched on the count rather than the wording, so
    // rephrasing the summary line does not fail this
    assert!(
        String::from_utf8_lossy(&run.stdout).contains('2'),
        "index did not report two records: {}",
        String::from_utf8_lossy(&run.stdout)
    );
}

// ---

/// Operations whose three backends must write the same bytes, with the
/// arguments to run each under.
//
// grep, filter, dedup and fetch are not here: their
// streaming and indexed paths copy a record's bytes through
// where the in-memory path re-wraps them, which is the one
// deliberate difference between the backends. --rewrap puts
// them back in step, which the test below checks
const AGREE: [&[&str]; 7] = [
    &["count"],
    &["names"],
    &["stats"],
    &["cat"],
    &["reformat"],
    &["head", "-n", "2"],
    &["tail", "-n", "2"],
];

#[test]
fn every_backend_of_an_operation_writes_the_same_bytes() {
    // a streaming or indexed path that disagrees with the
    // in-memory one is a bug, and no test that runs one
    // backend at a time would catch it
    for f in &FIXTURES {
        let path = fixture(f.file);

        for op in AGREE {
            let args = |mode: &str| {
                let mut argv: Vec<String> = op.iter().map(|a| a.to_string()).collect();
                argv.insert(1, path.clone());
                argv.push("--read".to_string());
                argv.push(mode.to_string());

                argv
            };

            let want = {
                let a = args("memory");
                let a: Vec<&str> = a.iter().map(String::as_str).collect();
                let run = sail(&a, b"");

                assert!(run.ok, "{op:?} --read memory failed: {}", run.stderr);
                run.stdout
            };

            for mode in ["stream", "indexed", "auto"] {
                let a = args(mode);
                let a: Vec<&str> = a.iter().map(String::as_str).collect();
                let run = sail(&a, b"");

                assert!(run.ok, "{op:?} --read {mode} failed: {}", run.stderr);
                assert_eq!(
                    run.stdout, want,
                    "{} {op:?} --read {mode} differs from --read memory",
                    f.stem
                );
            }
        }
    }
}

#[test]
fn rewrap_puts_the_selecting_operations_back_in_step() {
    for f in &FIXTURES {
        let path = fixture(f.file);

        // written out per operation rather than assembled,
        // because grep takes its pattern before its input and
        // the others take the input first
        let ops: [Vec<String>; 4] = [
            vec!["grep".into(), f.pattern.into(), path.clone()],
            vec!["filter".into(), path.clone(), "--min".into(), f.min.into()],
            vec!["dedup".into(), path.clone()],
            vec!["fetch".into(), path.clone(), f.name.into()],
        ];

        for op in ops {
            let run = |extra: &[&str]| {
                let mut argv: Vec<String> = op.clone();
                argv.extend(extra.iter().map(|a| a.to_string()));

                let a: Vec<&str> = argv.iter().map(String::as_str).collect();
                let out = sail(&a, b"");

                assert!(out.ok, "{a:?} failed: {}", out.stderr);
                out.stdout
            };

            assert_eq!(
                run(&["--rewrap"]),
                run(&["--read", "memory"]),
                "{} {op:?} --rewrap differs from the in-memory path",
                f.stem
            );
        }
    }
}

#[test]
fn indexed_refuses_a_pipe_rather_than_quietly_streaming() {
    // --read asserts a method the way --format asserts a
    // content type, so a pipe is an error and not a
    // downgrade -- a benchmark of a backend that never ran
    // is worse than no benchmark
    let bytes = std::fs::read(fixture("proteins.fa")).unwrap();
    let run = sail(&["count", "-", "--read", "indexed"], &bytes);

    assert!(!run.ok, "--read indexed accepted stdin");
    assert!(
        run.stderr.contains("byte offsets"),
        "the refusal does not say why: {}",
        run.stderr
    );
}

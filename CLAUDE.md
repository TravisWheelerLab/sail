# sail

A cargo workspace of two crates. `libsail` reads the file formats of profile-HMM
homology search and addresses the records inside them; `sail` is the
command-line tool over it. `fixtures/` at the workspace root holds a small file
for each format, shared by both crates' tests.

Design rationale that the code does not state lives in foam's memories, and
pending work lives in foam's issues. `foam prime` prints both at session start.

## the boundary between the crates

**The CLI never builds its own offset, index, or record data structures.** Those
belong in `libsail`. `sail`'s job is to choose among the backends `libsail`
already offers and to expose that choice as flags.

When the CLI needs something `libsail` lacks, the answer is to add it to
`libsail`, even when the CLI-side version would be shorter. Before changing
`libsail`'s public API, write the proposal up as a foam issue, a `decision`
where a choice is open, and wait for review rather than editing on the
assumption it will be accepted. The same goes for a change to `sail`'s command
line. Every issue carries the label of the crate it changes, `libsail` or
`sail`, and a milestone that spans both carries both.

The rule concerns state proportional to the *input*. It is not a rule about
memory in general. The CLI may hold what it is about to write or what the caller
asked for: a set of requested names, a ring buffer sized by `-n`, a vector sized
by the largest position requested. It may not hold a structure that grows with
the file, and it may not hold a byte offset.

## api design

**No parameter the library could derive itself.** When a value has exactly one
correct source and passing the wrong one fails silently, the library computes
it. `Index::write(path, data)` stats the data file; it does not take a stamp.
Prefer the signature that cannot be called wrongly over the one that documents
how to call it rightly, and treat a correctness argument for the riskier shape
as a reason to find a third design.

**Rejected machinery comes out whole.** When something is judged not worth its
cost, remove it rather than keep it in a cheaper form. Raise a concern only if
removing it breaks something the task needs.

**Read order before memory.** Judge any combinator that selects or reorders
records first by the order of the reads it causes, then by what it allocates.

## vocabulary

The word "sidecar" is banned from this project: source, comments, tests,
markdown, commit messages and conversation alike. The thing it used to name is
the `.saidx` file that `Index::write` puts beside a data file. Call it the index
file, or name the extension.

## conventions

**Formatting.** Stock rustfmt, with no configuration. This workspace has no
`rustfmt.toml`. Do not add one: the code is written against rustfmt's defaults.

**Comments.** Follow the `rust-comments` skill, and load it before writing or
editing a comment rather than reconstructing the rules from memory. The default
is no comment.

**Tests.** `libsail` keeps every test inline in the file it covers; `sail` also
has `tests/cli.rs`, which drives the real binary so that argv parsing, exit
codes and stderr are pinned alongside the output. The golden files under
`sail/tests/golden/` are regenerated with `SAIL_BLESS=1`. Every backend of an
operation must write identical bytes: a streaming or indexed path that disagrees
with the in-memory one is invisible to any single-backend test.

Test names are full sentences describing the property, not `test_foo`. For
instance `a_delimiter_mid_line_is_not_a_record_boundary`,
`position_returns_the_first_of_a_repeated_name`, and
`an_index_saved_beside_a_file_is_refused_once_the_file_changes`.

**Fixtures.** `fixtures/` sits at the workspace root, reached through the
`fixture!` macro, which resolves `CARGO_MANIFEST_DIR/../fixtures` at compile
time because cargo runs a test binary from the package root rather than the
workspace root.

**Branches.** Work is committed on `dev`, in batches. A milestone or issue that
changes a public API or an on-disk format, and will not be finished before the
next release, goes on a branch named for it. That branch merges into `dev` once
the work is whole, and takes merges from `dev` when the two drift apart. Nothing
else branches off `dev` unless Jack asks. `main` takes a merge from `dev` only
when a crate is being released. Unreleased code of the other crate may come
along, and it reaches users only when its own crate is released.

**Releases.** Jack says when a release is, never Claude. Until he does, versions
in `Cargo.toml` stay where they are and changelog entries accumulate under
`## [Unreleased]` with `### Added` / `### Changed` / `### Removed` headings and
no date. The version heading, the date, the `Cargo.toml` bumps and the compare
link at the foot of the changelog are all part of the release commit. Release
tags are annotated and named `<crate>-v<version>`. Publish from `main` after
`dev` merges in, with a clean tree.

**Pushing `main`.** `scripts/main-gate.sh` runs before every push. It refuses
a push that moves `main` unless the pushed commit passes all of these:

- It carries an annotated `<crate>-v<version>` tag.
- That version matches the crate's `Cargo.toml`, and the crate's changelog has
  a heading for it.
- For a `libsail` tag, `cargo semver-checks` passes against the previous
  `libsail` tag.
- `cargo fmt --check`, `cargo clippy` with warnings denied, and `cargo test`
  pass.

The gate checks the pushed commit in a worktree under `target/main-gate`, so
the branch checked out does not matter. It refuses every push to `main` when
`cargo-semver-checks` is not installed. Git never copies hook settings into a
clone, so each clone registers the gate once:

```console
$ git config hook.main-gate.event pre-push
$ git config hook.main-gate.command scripts/main-gate.sh
```

**Versions.** Each crate is versioned by what changed in it, and the two are
never aligned for tidiness. A break in `libsail` moves `libsail`'s minor, under
cargo's rule that 0.x's minor is the breaking axis, and moves nothing in `sail`
unless `sail`'s own interface changed. For the CLI that interface is the command
line and what each subcommand writes, not the library underneath, and
`sail/tests/cli.rs` is what says whether it moved. Before a release,
`cargo semver-checks -p libsail` compares `libsail` with its last published
version and reports whether the bump is a patch or a minor.

**Formats.** A change to an on-disk layout bumps that format's version byte.
The reader keeps an arm for the old layout by default. A breaking release may
drop that arm, listed under `### Removed`.

**Benchmarks.** `benchmarks/` is a separate cargo project rather than a
workspace member. Its `setup.sh` fetches pinned builds of seqkit and the Easel
miniapps into `benchmarks/tools`, and nothing is installed outside that
directory: no package manager, no `~/.local/bin`, no PATH edits. The versions
are pinned because they are part of the measurement. Every measured run reads
from disk; `--warm` keeps the page cache, and the board's header says which it
was, since the two are not comparable.

## checks

```console
$ cargo fmt --all
$ cargo fmt --all --check
$ cargo clippy --workspace --all-targets --all-features
$ cargo test --workspace --all-features
```

# sail

A cargo workspace of two crates. `libsail` reads the file formats of profile-HMM
homology search and addresses the records inside them; `sail` is the
command-line tool over it. `fixtures/` at the workspace root holds a small file
for each format, shared by both crates' tests.

## the boundary between the crates

**The CLI never builds its own offset, index, or record data structures.** Those
belong in `libsail`. `sail`'s job is to choose among the backends `libsail`
already offers and to expose that choice as flags.

When the CLI needs something `libsail` lacks, the answer is to add it to
`libsail`, even when the CLI-side version would be shorter. Before changing
`libsail`'s public API, write the proposed change up in `API-CHANGES.md` and
wait for review rather than editing on the assumption it will be accepted.

The rule concerns state proportional to the *input*. It is not a rule about
memory in general. The CLI may hold what it is about to write or what the caller
asked for: a set of requested names, a ring buffer sized by `-n`, a vector sized
by the largest position requested. It may not hold a structure that grows with
the file, and it may not hold a byte offset.

## libsail

Organised by concern rather than by format. A format module owns its record
type, its `Delimit` impl and its `Parse` impl; everything downstream sees only
the shared traits.

| module | owns |
| --- | --- |
| `collection.rs` | `Indexable`/`Iterable`, the `View`/`Iter`/`Aggregate` adapters, every combinator |
| `error.rs` | `Error`, `ErrorKind`, `bail_parse!` |
| `format.rs` | `Format` and `detect`/`detect_path` |
| `index.rs` | `Delimit`, `scan`/`frame`/`Offset`, `Index` and the `.saidx` file, the streaming `Reader` |
| `lines.rs` | row framing for hit tables |
| `parse.rs` | the `Parse` trait |
| `pread.rs` | positioned reads, private |
| `seq/` | FASTA, Stockholm, HMMER3 p7hmm |
| `source.rs` | `Source`, `FileSource`, `MmapSource` |
| `tbl/` | `Hit`, `HitColumns`, `WriteRow`, `Table`/`IndexedTable`, and the BLAST, HMMER and nail layouts |

`Format` covers only the three sequence formats. Hit tables are framed by
`lines.rs` and parsed through `tbl/`, and are not a `Format` variant.

### the core abstractions

`Source` is bytes addressable at an offset, with no cursor: `range(&self, start,
len) -> Cow<[u8]>`. It takes `&self` because `Indexable::get` does, so a
`Read`/`Seek` bound would put a `RefCell` or a `Mutex` behind every disk-backed
collection. It returns `Cow` so a mapped source can borrow and copy nothing,
which is most of what it exists for. `FileSource` pays a syscall and an
allocation per range; `MmapSource` is behind the `mmap` feature.

`Delimit` is how a format marks a record boundary. An impl supplies `DELIM` and
`DELIM_TYPE`; `WIDTH`, `CARRY` and a const-evaluated validity check are derived,
so an impl can neither forget the newline anchor nor apply it twice. `scan`
finds every anchored delimiter across chunked reads, carrying `CARRY` bytes over
each chunk boundary so a split delimiter is never lost or doubled. `frame` turns
those positions into `Offset { start, n_bytes }`, and `build_index` composes the
two.

`lines.rs` looks like it should be a `Delimit` impl and deliberately is not: a
row's delimiter is the newline itself, which `Delimit` rejects, and skipping
comment lines needs a line's first byte, which `frame` never sees.

`Parse` takes the bytes of one framed record and returns a record. It has no
receiver, for the same reason `Source::range` takes `&self`.

`Indexable` is the one trait an implementor writes: `len` and `get(n) ->
Option<Yield<'_>>`. `Yield<'a>: Borrow<Record>` is the load-bearing piece, and
the part newcomers misread. It has nothing to do with mutability. It is what
lets an in-memory collection lend `&Record` while a disk-backed one yields an
owned `Record` freshly parsed, through one API. Every combinator is agnostic to
which.

`Indexable` carries a hard invariant: `get(i)` is `Some` for every `i < len()`.
Breaking it panics rather than returning an error, because an implementor that
lies is a bug and not a runtime condition to recover from.

`Iterable` is blanket-implemented over every `Indexable`. Nobody writes one by
hand.

### the three ways to read a file

- **In memory** (`Fasta`, `Stockholm`, `Hmm`) reads the file whole, frames it,
  and parses every record up front. Cheapest per access, most expensive in
  memory: the bytes and the parsed records are both resident.
- **Indexed** (`IndexedFasta<S>` and siblings, generic over `Source`) builds the
  offset table in one pass and parses nothing. `get(n)` reads and parses on
  every call. Memory is the offset table alone, and it supports random access.
  The `Mmap*` aliases are the same types with `S = MmapSource`.
- **Streaming** (`index::Reader<R>`) keeps no offset table at all, advancing one
  record at a time through a reused buffer. Memory is proportional to the
  largest record rather than to the file. Forward-only, one pass.

The indexed types are `cfg(any(unix, windows))`, and the `Mmap*` ones are
further behind the `mmap` feature.

### `Index` and the `.saidx` file

`Index` holds a format, a vector of offsets, and an optional names layer: a blob
of concatenated name bytes, a span per record, and a sorted permutation for
lookup by name. An absent name stays distinct from an empty one, which Stockholm
without `#=GF ID` needs.

On disk: magic `SLIX`, a version byte, a format byte, the source file's length
and mtime, the record count, a names flag, then the offsets and the names layer,
all little-endian. No length read from the file is allocated against until the
bytes back it, and a truncated file is always an error rather than a panic.

`save` stamps the source's length and mtime and writes through a temporary file
and a rename; `load` re-stats the source and refuses an index file that disagrees.
This catches ordinary edits, not all of them: a same-size change within one
mtime tick goes undetected, which is a known gap rather than an oversight.
`path_for` appends `.saidx` rather than substituting, so `.fa` and `.sto`
siblings do not collide.

`build` scans; `build_named` scans and names in the same pass, reading each name
out of the reader's own buffer instead of going back to the file.

### combinators

Every combinator is a default method on `Indexable`, and every one is eager: it
walks `0..len()` once, builds a vector of indices, and returns a `View<I>`
holding the inner collection and that permutation. `View` is itself `Indexable`,
so combinators chain. `&I` is `Indexable` too, which is what lets one collection
be adapted more than once without being consumed or reopened.

Sorting and grouping materialize the key once per record rather than per
comparison, because `get` may be a disk read and a parse. `group_by` returns an
`IndexMap` so grouping hits by query reports queries in file order.

`Aggregate<I>` addresses several collections of the same type as one range,
dispatching through a binary search over per-part starts, one `usize` per part
rather than one per record. It nests, and `locate` recovers which part a record
came from.

### errors

Two variants, `Parse { detail }` and `Io { source }`, held inline rather than
boxed, since no `Result` here sits on a per-line path. `ErrorKind` is
`#[non_exhaustive]`. `bail_parse!` is the sugar for returning a parse error, and
position is formatted into the message at the call site rather than carried in
structured fields.

## sail

`main.rs` parses and runs. `cli.rs` holds the clap tree, the dispatch match over
all nineteen subcommands, and the shared argument enums. `input.rs` locates,
identifies and opens inputs, and decides the backend. `output.rs` opens the
writer and writes records. `commands/` is one module per subcommand.

There is no command trait. Each module exports a clap `Args` struct with `pub fn
run(self) -> Result<()>`, and `cli.rs` calls it. Inside `run` the shape is: plan
the inputs, open a writer, then match on the settled backend, with the
`dispatch!` and `indexed!` macros in `input.rs` selecting the per-format libsail
type inside each arm.

### choosing a backend

Three things meet in `resolve()`, in `input.rs`:

- **`Mode`** is what the caller asked for: `auto`, `stream`, `memory`,
  `indexed`. `auto` is the only one that is not a backend, and `resolve` turns
  it into one.
- **`Needs`** is what the operation requires of its input: `Pass` for one
  forward scan, `Count` when the total must be known first, `Sparse` for a few
  records by position, `Whole` when every record is needed in an order the file
  does not hold.
- **The input itself**, since a pipe has no byte offsets.

`resolve` is one table rather than logic spread across nineteen `run` methods,
so "fastest by default" is a claim a reader can check in one place.

An explicit `--read` is an assertion and never falls back. `--read indexed` on a
pipe is an error, not a downgrade, because downgrading to memory would do the
one thing the caller was avoiding and downgrading to streaming would report a
number for a backend that never ran. `--format` behaves the same way: it
asserts, and fails when the bytes disagree.

Four commands take no `--read` because they need every record at once: `sort`,
`shuffle`, `sample`, `rename`. Three more take none for structural reasons:
`split` reindexes the whole collection, `validate` needs the disk-backed type to
name the record that failed, and `index` writes byte offsets and so refuses
stdin outright.

### input and output

The format is read from the first 64 bytes, since every format's magic is on the
first line. Stdin is refused at a terminal, so a bare `sail count` does not
hang, and is otherwise read whole, because a pipe cannot be reopened. That is
the only reason. Every other input is identified with a short read and handed on
unopened, so nothing is buffered twice.

`FormatArg` mirrors libsail's `Format` as a separate clap-facing type so that
libsail need not depend on clap. The `From` impl is the whole of the mapping.

Records reach the writer one of two ways. The in-memory path writes parsed
records through a per-format function pointer, which is a pointer rather than a
trait because the three formats take different extra arguments. The streaming
and indexed paths write a record's framed bytes straight through, on the grounds
that an operation which only selected a record has nothing to say about how it
was wrapped; `--rewrap` asks for the normalising behaviour instead, and matters
only for `grep`, `filter`, `dedup` and `fetch`.

Nothing is written through `Display`. A `Formatter` takes `&str` and a record
name need not be UTF-8, and a `Display` impl that fails panics out of `write!`
rather than returning an error.

## conventions

**Formatting.** Stock rustfmt, with no configuration. This workspace has no
`rustfmt.toml`. Do not add one: the code is written against rustfmt's defaults.

**Comments.** Follow the `rust-comments` skill, and load it before writing or
editing a comment rather than reconstructing the rules from memory. The default
is no comment.

**Tests.** `libsail` keeps every test inline in the file it covers; `sail` also
has `tests/cli.rs`, which drives the real binary so that argv parsing, exit
codes and stderr are pinned alongside the output. Roughly a hundred golden files
live under `sail/tests/golden/`, regenerated with `SAIL_BLESS=1`. Every backend
of an operation must write identical bytes: a streaming or indexed path that
disagrees with the in-memory one is invisible to any single-backend test.

Test names are full sentences describing the property, not `test_foo`. For
instance `a_delimiter_mid_line_is_not_a_record_boundary`,
`position_returns_the_first_of_a_repeated_name`, and
`an_index_saved_beside_a_file_is_refused_once_the_file_changes`.

**Fixtures.** `fixtures/` sits at the workspace root, reached through the
`fixture!` macro, which resolves `CARGO_MANIFEST_DIR/../fixtures` at compile
time because cargo runs a test binary from the package root rather than the
workspace root.

**Versions.** Each crate is versioned by what changed in it, and the two are
never aligned for tidiness. A break in `libsail` moves `libsail`'s minor, under
cargo's rule that 0.x's minor is the breaking axis, and moves nothing in `sail`
unless `sail`'s own interface changed. For the CLI that interface is the command
line and what each subcommand writes, not the library underneath, and
`sail/tests/cli.rs` is what says whether it moved.

**Working notes.** `API-CHANGES.md` and `CLI-CHANGES.md` are scratch. They hold
proposals being thought through, they are not tracked, and they are never
committed. Provenance lives in commit messages and in the code, so a proposal
here can be rewritten or thrown away without ceremony.

**Benchmarks.** `benchmarks/` is a separate cargo project, not a workspace
member, and is not tracked. Its `setup.sh` fetches pinned builds of seqkit and
the Easel miniapps into `benchmarks/tools`, and nothing is installed outside
that directory: no package manager, no `~/.local/bin`, no PATH edits. The
versions are pinned because they are part of the measurement: a number is
comparable only against another taken with the same peer.

## checks

```console
$ cargo fmt --all
$ cargo fmt --all --check
$ cargo clippy --workspace --all-targets --all-features
$ cargo test --workspace --all-features
```

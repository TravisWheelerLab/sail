# libsail API changes

Everything here is additive. No existing signature changes and nothing is
removed, so every current caller keeps compiling. All of it is implemented;
`cargo test --workspace` covers each item.

## Why

The CLI cannot walk a file without materializing all of it. `Fasta::new` parses
every record before it returns one, and `IndexedFasta` hands back only parsed
records — its `source` is private, so nothing outside libsail can reach a
record's framed bytes. The public `name_of(record: &[u8])` therefore has no
public way to be fed.

The cost is the same on every operation, because every operation reads the same
way. On a 1 GiB FASTA of 2.68M records, `sail count` and `sail sort` both spend
~1.3s of user time and 3.7 GiB of resident memory before doing their own work.
`sail index`, which uses the streaming scanner instead, crosses the same file in
0.21s user at 152 MiB.

Of the 19 subcommands, two need every record at once. Ten are a single forward
pass. None of the FASTA operations need a materialized residue `Vec` at all:
they need a name, a length, or the record's bytes to copy through.

## `index::Reader<R: Read>`

One record's bytes at a time, framed out of a stream, into a buffer reused
between calls.

```rust
pub struct Reader<R: Read> { /* private */ }

impl<R: Read> Reader<R> {
    pub fn new(bytes: R, format: Format) -> Self;
    pub fn detect(bytes: R) -> Result<Self>;

    /// Advance to the next record. `false` at the end of the stream.
    pub fn advance(&mut self) -> Result<bool>;

    /// The current record's bytes; empty before the first `advance`.
    pub fn record(&self) -> &[u8];

    pub fn offset(&self) -> Offset;
    pub fn format(&self) -> Format;
}
```

Advance-then-borrow rather than `Iterator`, because `Iterator::next` cannot lend
out of `self`, and an owned record per step is the allocation this exists to
avoid.

`Format` is a field rather than a `D: Delimit` type parameter. `Delimit` has
three impls, and it earns nothing in the inner loop: `scan` already builds its
needle at runtime and `memmem::Finder` already dispatches at runtime, so `D`
contributes one predicted branch per record — roughly 0.7 ms against a 280 ms
pass. As a type parameter it would instead cost four commands a three-arm macro
expansion each, of a body that is otherwise identical.

Making `scan` and `frame` public was the alternative. It does not work: `scan`
returns every delimiter position, 33.5 MiB on this fixture, so it is not the
O(1) walk it would be exported to provide.

```rust
let mut reader = Reader::new(File::open(path)?, Format::Fasta);
while reader.advance()? {
    out.write_all(reader.record())?;
}
```

## `seq::name_of`, `fasta::len_of`, `fasta::rewrap`

```rust
// libsail/src/seq/mod.rs
pub fn name_of(format: Format, record: &[u8]) -> Option<&[u8]>;

// libsail/src/seq/fasta.rs
pub fn len_of(record: &[u8]) -> usize;
pub fn rewrap(record: &[u8], w: &mut dyn Write, width: usize) -> Result<()>;
```

The format-dispatching `name_of` sits over the three that already exist. It is
what lets a name-matching command be one loop rather than three:

```rust
while reader.advance()? {
    let record = reader.record();
    if pattern.is_match(seq::name_of(format, record).unwrap_or_default()) {
        out.write_all(record)?;
    }
}
```

`len_of` reports residues without parsing, which is all `stats` and `filter`
need. FASTA only: parsing cost scales with record count, and FASTA is the only
one of the three formats that has millions of them in a file. Adding
`depth_of`, `width_of` and `leng_of` to match would be three functions no
measurement asks for.

`rewrap` writes a record's header and its residues at a given line width,
straight from the framed bytes. Without it `reformat` has to build a
`FastaRecord` to rewrap it, which is 0.37s of the budget on the fixture.

## Record bytes from the indexed backends

```rust
// on IndexedFasta<S>, IndexedStockholm<S>, IndexedHmm<S>
pub fn record(&self, n: usize) -> Result<Option<Cow<'_, [u8]>>>;
```

The parsed `try_get` is the only way in today, so an indexed read cannot emit
bytes without a round trip through the parser. This keeps the offset arithmetic
inside libsail, where `Source::range` already is.

## `Index::build` and `Index::build_named`

```rust
impl Index {
    pub fn build<R: Read>(bytes: R, format: Format) -> Result<Index>;
    pub fn build_named<R: Read>(bytes: R, format: Format) -> Result<Index>;
}
```

`name_index` currently names records through `Source::range`, which is a fresh
allocation and a `pread` each. Naming them during the scan, while the bytes are
already in the buffer, removes 2.68M reads and 5.4M allocations from `sail
index` and collapses the three-arm match in `commands/index.rs` to one line.

## Not an API change, but going in with it

- `FastaParser::parse` builds `seq` by collecting a filter whose `size_hint`
  lower bound is zero, so it doubles from 8 with a full regrowth copy each time.
  `Vec::with_capacity(residues.len())` overshoots only by the newline count.
  ~354 MiB and ~16M allocations on the fixture. Same in `Stockholm::new` and
  `Hmm::new`.
- `Fasta::new`'s `collect::<Result<Vec<_>>>()` has the same problem in the
  record spine: 302 MB of capacity for 193 MB of records.
- `Fasta::new` opens with `read_to_end`, which copies a slice the caller already
  owns when it is handed one. The CLI hands it one on every command.
- `Indexable::cloned` deep-copies a second time on the disk-backed collections,
  whose `Yield` is already owned. `Fasta` lends, so its clone is the necessary
  one; the three `Indexed*` impls should override.

## Also shipped, in `sail`

- `--read auto|stream|memory|indexed` on every converted operation, resolved in
  one function against what the operation needs. `--read indexed` on a pipe is
  an error, not a downgrade.
- `fetch --order file|asked`, file order by default.
- `--rewrap` on `grep`, `filter`, `dedup` and `fetch`, which otherwise copy a
  selected record's bytes through rather than re-wrapping them.

## What this does not solve

Peak memory becomes O(largest record) rather than O(1). For FASTA that is the
64 KiB buffer. A single large Stockholm alignment is still the whole alignment,
and `IndexedStockholm` remains the answer there.

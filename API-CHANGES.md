# libsail API changes

The first three items are additive and shipped in libsail 0.2.1;
`cargo test --workspace` covers each. The fourth changes `Index` and is still
proposed, for a 0.3.

## Why

Three things `nail-benchmarks` writes by hand against libsail 0.2, and one thing
`Index` gets wrong on its own.

`scores::scan` walks a hit table one row at a time, because libsail offers no way
to. A sequence file streams through `index::Reader`; a hit table has only the
in-memory `Table` or the random-access `IndexedTable`.

`util::split` deals records into `k` files weighted by residues or model length,
so that the files take about the same time to search rather than holding the same
number of records. libsail has `split(k)`, which balances counts.

`Index` reads every offset whether a caller wants one record or all of them, and
its file API asks for the data file in order to read an index that describes it.

## `lines::Rows<R: Read>`

One table row at a time, framed out of a stream, into a buffer reused between
calls. `index::Reader`'s shape, for the other framing.

```rust
pub struct Rows<R: Read> { /* private */ }

impl<R: Read> Rows<R> {
    /// `skip` sees a line's first byte, as `line_index`'s does.
    pub fn new(bytes: R, skip: fn(u8) -> bool) -> Self;

    /// Advance to the next row. `false` at the end of the stream.
    pub fn advance(&mut self) -> Result<bool>;

    /// The current row's bytes, without its newline.
    pub fn row(&self) -> &[u8];

    /// The line this row came from, counting from 1, for an error that has
    /// to name it.
    pub fn line(&self) -> u64;

    pub fn offset(&self) -> Offset;
}
```

`skip` is a parameter rather than a hardcoded `#`, so the two framings in
`lines` agree on what a row is. A line it rejects is skipped wherever it occurs,
not only at the head of the file: a table concatenated from shards repeats its
header partway down, and a reader that skipped only a leading block would parse
those as rows.

`Rows` frames a file exactly as `line_index` does, which a test asserts across
every hit-table fixture at five chunk sizes. That costs one non-obvious rule:
emptiness is judged before a trailing `\r` is dropped, because `line_index`
keeps a `\r\n` line and trimming first would silently skip it. `row()` drops the
line ending, `offset()` spans it.

## `tbl::fields`

```rust
pub fn fields<const K: usize>(row: &[u8], want: [usize; K]) -> Option<[&[u8]; K]>;
```

The columns at `want`, borrowed from the row. `None` when the row holds too few.

```rust
let mut rows = Rows::new(File::open(path)?);
while rows.advance()? {
    let Some([query, target, score]) =
        fields(rows.row(), [C::QUERY, C::TARGET, C::SCORE])
    else {
        continue;
    };
}
```

`K` is a const parameter against the counting rule, which says a type parameter
with fewer than five instantiations wants an enum or a field instead. It earns
its place on a different axis: the return is a fixed array rather than a `Vec`,
so the walk allocates nothing, and one pass reads every wanted column and stops
at the highest index rather than tokenising the rest of the row. The alternative
that avoids `K` is `field(row, n) -> Option<&[u8]>`, called once per column,
which is a pass per column.

Whitespace splitting only. There is no quoting and no escaping, because the three
table formats have neither.

## `Indexable::split_weighted`

```rust
fn split_weighted<F>(&self, k: usize, weight: F) -> Vec<View<&Self>>
where
    F: FnMut(&Self::Record) -> u64;
```

`k` parts of about equal weight, beside `split(k)`'s `k` parts of about equal
length. Longest-processing-time: the records sort by weight descending, and each
goes to the lightest part so far, ties to the lowest index, so one input always
deals the same way.

```rust
let parts = hmms.split_weighted(8, |m| m.header.leng as u64);
```

The weight is `u64` rather than a type parameter. Every weight in play is a
length: residues, model positions, alignment columns. A bound set of
`Ord + Default + Add` buys nothing for it.

The closure is also what makes this work for Stockholm, which the benchmark's
version does not handle. That one is a two-variant enum over FASTA and p7hmm, so
a third format means editing it.

Like `split`, this returns exactly `k` parts, some of them possibly empty. The
benchmark's version clamps to the record count instead. Agreeing with `split`
matters more than agreeing with the caller being replaced.

## `Index`: offsets that stay on disk, and a file API that reads one file

Two problems, one change.

`Index::load` reads every offset, 16 bytes a record, whether the caller wants one
record or all of them. There is no way to ask for fewer.

And the file API asks for the wrong thing. `load(source)` takes the *data* file
and infers the index path from it, because the staleness check is welded into
deserialisation. Reading a list of offsets and checking that those offsets still
describe something are two jobs, and one function does both. The names pay for
it: `source` also names a trait in `source.rs` meaning something else entirely,
`load_from`'s second path is called `at`, and `load_from` and `save_as` name one
axis two ways.

### The file API

```rust
// module level, because these are about files rather than values
pub fn path_for(data: &Path) -> PathBuf;
pub fn is_current(path: &Path, data: &Path) -> Result<bool>;
pub fn count_in(path: &Path) -> Result<usize>;

impl Index {
    /// Read the whole index into memory.
    pub fn read(path: &Path) -> Result<Index>;

    /// Attach to the index, leaving its offsets in the file.
    pub fn open(path: &Path) -> Result<Index>;

    /// Write this index to `path`, stamped with `data` as it is now.
    pub fn write(&self, path: &Path, data: &Path) -> Result<()>;
}
```

`read` and `open` take one path, the index, and never stat the data. `write`
takes two, because the stamp it records is a fact about the data file and cannot
be written without looking at one. That asymmetry is the whole of it: reading an
index needs nothing but the index, writing one needs the thing it describes.

`write` stats `data` itself rather than taking a prepared stamp. Handing the
caller a stamp to produce and pass back is a value with exactly one correct
source and no way to check it, so a caller who takes it from the wrong file, or
at the wrong moment, gets an index that looks valid and is not. The stamp type
stays private and `is_current` asks the only question anyone had for it.

There is no `load`/`load_from` or `save`/`save_as` pair. A caller wanting the
default location asks `path_for` for it.

```rust
let path = index::path_for(&data);

let index = Index::build_named(File::open(&data)?, Format::Fasta)?;
index.write(&path, &data)?;
```

```rust
if index::is_current(&path, &data)? {
    let index = Index::open(&path)?;
}
```

### Offsets that stay on disk

A private enum inside `Index`, which stays one concrete public type:

```rust
enum Offsets {
    Memory(Vec<Offset>),
    /// The index file. Offset `n` is the 16 bytes at `at + n * 16`.
    File { file: fs::File, at: u64, count: usize },
}
```

Not a type parameter: `IndexedFasta<S, O>` would put a second axis on every
collection and every error message, and could never be taken back. Not a boxed
trait either, which costs a vtable and an allocation to buy extensibility that a
private trait does not actually provide, since only this crate can implement
either. Because the public API is the same under all three, the choice stays
internal: if a foreign index format ever turns up, this becomes a public trait
without changing a caller.

The `File` arm reads through `pread::read_exact_at`, already in the crate.

Two signatures change, and `Offset` gains `Copy`:

```rust
pub fn offset(&self, n: usize) -> Result<Option<Offset>>;   // was Option<&Offset>
pub fn offsets(&self) -> Option<&[Offset]>;                 // was &[Offset]
```

`Offset` keeps `#[repr(Rust)]` and `n_bytes: usize`, and the file format does not
change.

### What it costs the callers

`offsets()` has no callers at all. `offset(n)` has three, one per format, each
already inside a function returning `Result`:

```diff
     pub fn record(&self, n: usize) -> Result<Option<Cow<'_, [u8]>>> {
-        let Some(o) = self.index.offset(n) else {
+        let Some(o) = self.index.offset(n)? else {
             return Ok(None);
         };
 
         Ok(Some(self.source.range(o.start, o.n_bytes)?))
     }
```

`Indexable::get` does not change, and nothing in the combinator layer does
either. It is already this:

```rust
fn get(&self, n: usize) -> Option<FastaRecord> {
    self.try_get(n).ok().flatten()
}
```

An indexed collection already reads and parses on every `get`, and already
collapses a failure into `None`. A file-backed index adds a second read to a call
that always had one. Had `get` been infallible, this design would have forced
`Result` through `Indexable` and every combinator with it.

No collection gains a constructor. `with_index` already takes whatever `Index` it
is handed and only checks the format:

```rust
IndexedFasta::with_index(&data, Index::read(&path)?)?   // offsets resident
IndexedFasta::with_index(&data, Index::open(&path)?)?   // offsets in the file
```

`Index` loses `Clone` and `PartialEq`. `fs::File` implements neither, so this is
a cost of holding an open file rather than of this design; a boxed trait pays it
too. `Clone` is unused. `PartialEq` has two call sites, both assertions in
`index.rs`'s own tests.

A memory-mapped arm was the alternative and costs more than it returns. It would
let a mapped region be borrowed as `&[Offset]`, keeping both signatures, but only
with `#[repr(C)]` and `n_bytes: u64` on a public struct, a little-endian target,
and five bytes of header padding so the offsets land 8-aligned, which is a v2
file format. It also buys nothing for sparse reads: a page fault costs about what
a `pread` costs, and the kernel reads ahead 128 KiB to use 16 bytes.

## Not settled

A file-backed index cannot lend names. `Names` holds `blob`, `spans` and
`by_name` as `Vec`s, and `name(n)` borrows out of `blob`. The simple answer is
for `open` to pull the name layer into memory anyway, since it is absent on the
indexes large enough to care. The alternative is for `open` to refuse a named
index outright.

`scan`'s buffer opens at 8 MiB and doubles to a 64 MiB cap, where
`index::Reader`'s opens at 64 KiB and grows only when a record will not fit. A
row is small and a record need not be, so the two want different numbers, and
nothing here has been measured.

The zero-allocation framing in `scan` describes the design rather than a result.
That repository holds no benchmark, no criterion dependency and no figure for
what the approach buys over `BufRead::read_line` and `split_whitespace`. The
measurement belongs before the buffer constants are chosen rather than after.

## What this does not solve

`Index::open` makes a sparse read cheap and leaves a full scan exactly where it
is. A scan should take `offsets()` once and walk the slice, which means it wants
`Index::read` and the whole allocation. That is the right trade for a scan
and the wrong one for a draw, and nothing chooses between them automatically.

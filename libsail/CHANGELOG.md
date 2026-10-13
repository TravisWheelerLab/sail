# Changelog

All notable changes to libsail are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html), with
cargo's reading of it below 1.0: the **minor** is the breaking axis, so 0.2 to
0.3 is a breaking release and 0.2.0 to 0.2.1 is not.

## [Unreleased]

## [0.4.1] - 2026-10-13

A streaming `Reader` reports the index of the record it holds and can advance
to a numbered one. `Readers` numbers the records of several files as one run.
A caller drawing a sample over a stream no longer writes its own counting and
position loops.

### Added

- `index::Framed<S>`, a file and the offset of every record in it, with no
  parser attached. `open(path, format)` scans, `detect(path)` reads the format
  from the first bytes, `with_index(path, index)` takes it from the index, and
  `map(path, format)` maps the file instead of reading it. `record(n)` hands
  back one record's framed bytes and `name_index` names them, so a caller with
  a parser of its own reads a file the way the disk-backed collections do.
- `index::count(bytes, format)`, the number of records in a stream, by one pass
  over it.
- `Reader::position()` and `Reader::advance_to(n)`. `position()` is the index of
  the record `record()` returns, and `None` before the first `advance`. An
  advance that returns `false` leaves the last record current, so
  `position() + 1` is then the record count. `advance_to(n)` advances until
  record `n` is current, and panics if record `n` is already behind the reader.
  `lines::Rows` gains the same two methods, counting rows.
- `index::Readers`, which reads several files as one run of records with one
  `position()`. Each file gets its own `Reader`, so a file without a final
  newline cannot swallow the first record of the next file as `Read::chain`
  does. `count()` keeps one count per file. Once it has run, `advance_to(n)`
  skips every file that ends before record `n` without opening it.

### Changed

- `Reader::record()` and `Rows::row()` still return the last record after
  `advance` returns `false`. On Stockholm and HMM, `Reader::record()` used to
  return bytes that the final refill had already moved.
- `IndexedFasta`, `IndexedStockholm` and `IndexedHmm` each hold a `Framed<S>`
  and forward `index`, `name_index` and `record` to it, where each used to
  carry its own copy of the source and the index. No signature moved.

## [0.4.0] - 2026-09-16

`permute` and `sample` compute their permutation from a seed instead of storing
one, which drops the `rand` dependency and makes a draw cost the same at any
length. The order records are read in is unchanged, and still the
permutation's. Allocation is the whole of the saving; the reads are as scattered
as they were.

### Added

- `select::Selection`, `m` of `0..n` drawn uniformly at random, ascending, in
  constant memory. Vitter's Method D while the gaps are long and Algorithm A
  once they are short, picking between them on the alpha crossover he reports.
  Reads nothing, so the choice is testable with no file behind it, and a caller
  crossing the gap between two selected records decides for itself whether to
  scan it or seek it.
- `Indexable::sample_in_order(m, seed)`, `m` records drawn uniformly at random
  and yielded in the order the collection holds them. A file-backed collection
  pays `m` ascending reads where `sample` pays `m` scattered ones, which is what
  decides the cost on a file too large to cache. Returns an iterator rather than
  a collection and allocates nothing, so it is additive: `sample` and `permute`
  are unchanged and the two are not interchangeable.
- `collection::Drawn<I>`, what `permute` and `sample` now return. Holds a
  `feisty::Permutation` over `0..len()` and a count of how many of its positions
  to keep. Allocates nothing per record. `permute_with` over 2.4 billion
  records built a 19.6 GB vector of indices before reading anything; `Drawn`
  builds a key.

### Changed

- `Indexable::permute(seed)` and `Indexable::sample(m, seed)` replace
  `permute`, `permute_with`, `sample` and `sample_with`. Both return
  `Drawn<Self>` where the old four returned `View<Self>`, and both take a `u64`
  where the `_with` pair took an `R: Rng`, so a draw reproduces from its
  recorded seed alone rather than from every prior use of a generator.
- A seed is required: the unseeded `permute()` and `sample(n)` are gone. With
  one variant of each combinator there is no `_with` suffix left to
  disambiguate.
- `sample(m, seed)` keeps the first `m` positions of `permute(seed)`, so two
  draws from one collection under one seed nest — the smaller is a subset of the
  larger. Two draws through `rand` were unrelated.
- A draw is now one of the 2^64 permutations a `u64` can key, where shuffling
  through a `ChaCha`-backed generator could reach far more of them. `sail`
  already funnelled `--seed` through `seed_from_u64`, so it never had more than
  2^64 either.

### Removed

- The `rand` dependency, and with it seven crates from every downstream build.
  `feisty` has none of its own. What a seed reproduces is now fixed by feisty
  rather than by `rand`'s version: `DEFAULT_ROUNDS` and the round function for
  `permute` and `sample`, and `splitmix64` for `sample_in_order`, which is the
  generator keying the other two. One seed therefore means one thing across all
  three. Needs feisty 0.1.1, where `splitmix64` became public.

## [0.3.0] - 2026-09-16

Records hold their text as bytes, an index can leave its offsets on disk, and a
hit table can be walked one row at a time.

### Added

- `lines::Rows<R: Read>`, one table row at a time out of a stream and into a
  reused buffer. Frames a file exactly as `line_index` does, which a test
  asserts over every hit-table fixture at five chunk sizes.
- `tbl::fields<const K>(row, want)`, K columns borrowed from a row in one pass,
  stopping at the highest column asked for.
- `Indexable::split_weighted(k, weight)`, `k` parts of near-equal weight beside
  `split`'s `k` parts of near-equal length. Longest-processing-time, and
  deterministic for a given input.
- `seq::p7hmm::leng_of`, a model's `LENG` read off its bytes without parsing it.
- `Index::open`, which reads an index file's header and leaves the offsets and
  names in it, answering `offset(n)` in one `pread` and `name(n)` in two.
- `Index::read`, `Index::write`, and module-level `index::path_for`,
  `index::is_current` and `index::count_in`.
- `StockholmRecord::id_str` and `HmmHeader::name_str`, for a caller that wants
  text and will take the `Result`.
- `StockholmRecord::gr_get(name, feature)`, which reads a `#=GR` annotation
  without building the key it is stored under.

### Changed

- Records hold bytes. `StockholmRecord`'s `gf`, `gs`, `gc`, `gr` and
  `names`, and `HmmHeader`'s `version`, `name`, `acc`, `desc`, `date`, `com`
  and `extra`, are `Vec<u8>` where they were `String`. FASTA already worked this
  way; the other two were stricter than either specification requires.
  A non-UTF-8 byte in an annotation no longer fails the record, which is what
  kept Pfam's SEED alignments from parsing at all. Numbers and enums still read
  as text and fail as the field they are.
- `Index::offset(n)` returns `Result<Option<Offset>>` where it returned
  `Option<&Offset>`, and `Index::offsets()` returns `Option<&[Offset]>`, `Some`
  only where there is a slice to lend.
- `Index::name(n)` returns `Result<Option<Cow<'_, [u8]>>>` and
  `Index::position(name)` returns `Result<Option<usize>>`.
- `StockholmRecord::get`, `id`, `gf_first`, `gf_all`, `ss_cons` and `row` take
  and answer in `&[u8]`.
- `Offset` derives `Copy`.

### Removed

- `Index::load`, `Index::load_from`, `Index::save` and `Index::save_as`,
  replaced by `read`/`open`/`write` plus `index::path_for`. Reading an index no
  longer needs the data file: `read` and `open` take one path, the index, and
  only `write` takes two, because the stamp it records is a fact about the data.
- `Index` no longer derives `Clone` or `PartialEq`. An index may hold an open
  file, and `std::fs::File` implements neither.

## [0.2.0] - 2026-09-16

The hit tables gained per-format row types, and every format gained a way to be
read one record at a time.

### Added

- `format::Format`, `format::detect` and `format::detect_path`.
- `index::Index`, an offset table with an optional name layer that saves and
  loads itself beside the file it describes.
- `index::Reader<R: Read>`, one record's framed bytes at a time.
- `seq::name_of` and the per-format `name_of` scans, `seq::fasta::len_of` and
  `seq::fasta::rewrap`.
- `record(n)` on the indexed collections, for a record's framed bytes.
- Per-format hit-table rows carrying every column: `tbl::blast::BlastRow`,
  `tbl::nail::NailRow`, `tbl::hmmer::HmmerRow` and `HmmerDomRow`, with the
  `tbl::WriteRow` trait and the `*Rows` aliases.

### Changed

- `tbl::Table` and `tbl::IndexedTable` take a `P: Parse` where they took a
  `C: HitColumns`, so `try_get` yields `P::Record`.
- `HitColumns::N_FIELDS` became `N_COLUMNS`.
- `collection::Indexed<I>` became `View<I>` and `IndexedIter` became `Iter`,
  which changes the return type of every combinator.

### Removed

- `tbl::Span`. `Hit`'s `query_span` and `target_span` became four
  `Option<usize>` fields.

## [0.1.0] - 2026-09-07

First release: FASTA, Stockholm, HMMER3 p7hmm and the BLAST, HMMER and nail hit
tables, with the `collection`, `parse` and `source` layers behind them.

[Unreleased]: https://github.com/TravisWheelerLab/sail/compare/libsail-v0.4.1...HEAD
[0.4.1]: https://github.com/TravisWheelerLab/sail/compare/libsail-v0.4.0...libsail-v0.4.1
[0.4.0]: https://github.com/TravisWheelerLab/sail/compare/libsail-v0.3.0...libsail-v0.4.0
[0.3.0]: https://github.com/TravisWheelerLab/sail/compare/libsail-v0.2.0...libsail-v0.3.0
[0.2.0]: https://github.com/TravisWheelerLab/sail/compare/libsail-v0.1.0...libsail-v0.2.0
[0.1.0]: https://github.com/TravisWheelerLab/sail/releases/tag/libsail-v0.1.0

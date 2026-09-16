# Changelog

All notable changes to the `sail` command line tool are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html), with
cargo's reading of it below 1.0: the **minor** is the breaking axis. What is
versioned is the command line and what each subcommand writes, not the library
underneath, so `sail` and `libsail` move independently.

## [Unreleased]

### Changed

- `shuffle` and `sample` write a different order for a given `--seed`. libsail
  computes the permutation from the seed now rather than shuffling a vector of
  indices, so a given seed picks different records. The flags are unchanged, a
  draw is still `-n` distinct records taken without replacement, and the memory
  it costs no longer grows with the input.
- `--seed` is the seed itself. Both commands used to build a generator from it
  and hand that over; they pass the number straight through now, and an
  unseeded run picks one `u64` for itself.

## [0.2.1] - 2026-09-16

### Changed

- `--read indexed` reads the index beside its input when one is there and still
  describes the file, and scans for one when it is not. `sail index` wrote a
  file nothing read back, so a second command on the same file paid the same
  scan as the first. The flags and the output are unchanged, which the golden
  tests pin byte for byte.
- `rename` can rename a record to a name that is not valid UTF-8, which `grep`
  and `fetch` already allowed.

## [0.2.0] - 2026-09-16

First release as a command line tool. What was on crates.io as 0.1.0 and 0.1.1
was a library stub with no binary.

### Added

- Nineteen subcommands over FASTA, Stockholm and HMMER3 profile HMMs: `cat`,
  `count`, `dedup`, `fetch`, `filter`, `get`, `grep`, `head`, `index`, `names`,
  `reformat`, `rename`, `sample`, `shuffle`, `sort`, `split`, `stats`, `tail`
  and `validate`.
- `--read auto|stream|memory|indexed` on most operations, choosing how the
  input is walked. An explicit mode is an assertion and never falls back, so
  `--read indexed` on a pipe is an error rather than a downgrade.
- `--format`, which asserts the input's format and fails on anything else.
- `fetch --order file|asked`, and `--rewrap` on `grep`, `filter`, `dedup` and
  `fetch`, which otherwise copy a selected record's bytes through untouched.
- `sail index`, which writes a record index beside a file.

### Fixed

- Writing no longer goes through `Display`, which panicked on a record whose
  name is not UTF-8 and ignored FASTA's line width.
- Reading from a terminal is refused rather than hanging.

## [0.1.1] - 2026-07-01

## [0.1.0] - 2025-11-05

Placeholders. Both were library stubs holding a few lines of Rust, and neither
shipped a command.

[Unreleased]: https://github.com/TravisWheelerLab/sail/compare/sail-v0.2.1...HEAD
[0.2.1]: https://github.com/TravisWheelerLab/sail/compare/sail-v0.2.0...sail-v0.2.1
[0.2.0]: https://github.com/TravisWheelerLab/sail/compare/sail-v0.1.1...sail-v0.2.0
[0.1.1]: https://github.com/TravisWheelerLab/sail/compare/sail-v0.1.0...sail-v0.1.1
[0.1.0]: https://github.com/TravisWheelerLab/sail/releases/tag/sail-v0.1.0

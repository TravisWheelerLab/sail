# libsail

A Rust library for the file formats of profile-HMM homology search: FASTA,
Stockholm, HMMER3 p7hmm, and the BLAST, HMMER, and nail hit tables.

## install

```console
$ cargo add libsail
```

## use

The crate is organised by concern rather than by format. `seq` and `tbl` hold
the record types, `format` tells the formats apart, `index` finds record
boundaries, `collection` addresses records by position, and `parse` reads text
back into records.

Walking a file one record at a time, without holding it in memory:

```rust
use libsail::{format::Format, index::Reader};
use std::{fs::File, io::Write};

let mut out = std::io::stdout().lock();
let mut reader = Reader::new(File::open("proteins.fa")?, Format::Fasta);
while reader.advance()? {
    out.write_all(reader.record())?;
}
```

Addressing records by position, reading each one off the disk on demand:

```rust
use libsail::seq::fasta::IndexedFasta;

let fasta = IndexedFasta::open("proteins.fa")?;
let record = fasta.try_get(0)?;
```

An `Index` saves and loads itself, so a file need not be scanned twice:

```rust
use libsail::index::Index;
use std::{fs::File, path::Path};

let index = Index::build_named(File::open("proteins.fa")?, Format::Fasta)?;
index.save(Path::new("proteins.fa"))?;
```

## features

`mmap` memory-maps a file instead of reading it through `pread`. Off by
default.

## license

BSD-3-Clause.

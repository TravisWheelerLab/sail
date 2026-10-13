# sail

> **<u>s</u>equence <u>a</u>nalysis <u>I</u>/O too<u>l</u>**

A command-line toolkit for biological sequence files. `sail` reads FASTA,
Stockholm, and HMMER3 profile HMMs, and gives all three the same set of
operations.

## install

```console
$ cargo install sail
```

## use

A command takes a file, or `-` for stdin, and writes to stdout, or to the path
given by `-o`. `sail` works out the format from the input; pass `--format` to
require a particular one and fail on anything else.

```console
$ sail count proteins.fa
$ sail names proteins.fa
$ sail grep -i kinase proteins.fa
$ sail fetch proteins.fa P01308 P01315
$ sail sort --key size -r proteins.fa -o sorted.fa
$ sail split -k 4 --prefix chunk proteins.fa
```

Twenty commands in all. `sail --help` lists them, and each takes its own
`--help`.

## reading large files

Most commands take `--read`, which sets how the command walks its input:

- `stream`, one record at a time through a reused buffer
- `memory`, every record parsed up front
- `indexed`, a byte offset per record, read on demand
- `auto`, the default, which lets the command decide

`sail index` writes a `.saidx` file beside its input, holding the offset and
name of every record, so a later run can address records without scanning the
file again.

## license

BSD-3-Clause.

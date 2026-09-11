use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use clap::Args;
use libsail::collection::{Indexable, Iterable};
use libsail::format::Format;

use crate::cli::FormatArg;
use crate::input::{Inputs, dispatch};
use crate::output::{emit, writer};

#[derive(Args)]
pub struct SplitArgs {
    /// file to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: PathBuf,

    /// split into exactly this many files
    #[arg(short = 'k', long, conflicts_with_all = ["size", "each"])]
    pub parts: Option<usize>,

    /// split into files of this many records
    #[arg(short = 'n', long, conflicts_with = "each")]
    pub size: Option<usize>,

    /// write one file per record
    #[arg(long)]
    pub each: bool,

    /// stem for the files written [default: the input's name]
    #[arg(long)]
    pub prefix: Option<PathBuf>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,
}

impl SplitArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::read(std::slice::from_ref(&self.input), self.format)?;
        let entry = &inputs.entries()[0];
        let prefix = self.prefix.clone().unwrap_or_else(|| self.stem());
        let extension = extension(inputs.format());

        let (parts, size, each) = (self.parts, self.size, self.each);

        dispatch!(inputs.format(), entry, |collection| {
            let groups = group_sizes(parts, size, each, collection.len())?;

            let mut at = 0;
            for (i, take) in groups.iter().enumerate() {
                let path = part_path(&prefix, i + 1, extension);
                let out = writer(Some(&path))?;

                emit((&collection).skip(at).take(*take).iter(), out)?;
                at += take;

                println!("{take}\t{}", path.display());
            }

            Ok(())
        })
    }

    /// The input's path with its extension dropped, or `part` for stdin.
    fn stem(&self) -> PathBuf {
        if self.input.as_os_str() == "-" {
            return PathBuf::from("part");
        }

        self.input.with_extension("")
    }
}

// ---

/// How many records go in each file written.
fn group_sizes(
    parts: Option<usize>,
    size: Option<usize>,
    each: bool,
    len: usize,
) -> Result<Vec<usize>> {
    if each {
        return Ok(vec![1; len]);
    }

    if let Some(k) = parts {
        if k == 0 {
            bail!("splitting into 0 files writes nothing")
        }

        // exactly k parts, empty ones included: a caller
        // sizing a thread pool by the count would otherwise
        // get fewer workers than it asked for
        let base = len / k;
        let over = len % k;

        return Ok((0..k).map(|i| base + usize::from(i < over)).collect());
    }

    if let Some(n) = size {
        if n == 0 {
            bail!("splitting into files of 0 records never finishes")
        }

        // by size, so an empty input yields no files at all
        return Ok(std::iter::repeat_n(n, len / n)
            .chain((!len.is_multiple_of(n)).then_some(len % n))
            .collect());
    }

    bail!("split needs -k, -n or --each")
}

fn part_path(prefix: &Path, n: usize, extension: &str) -> PathBuf {
    let mut name = prefix.as_os_str().to_os_string();
    name.push(format!(".{n:04}.{extension}"));

    PathBuf::from(name)
}

/// The extension each format's files are written with.
fn extension(format: Format) -> &'static str {
    match format {
        Format::Fasta => "fa",
        Format::Stockholm => "sto",
        Format::Hmm => "hmm",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splitting_into_k_parts_always_writes_exactly_k_files() {
        // even when that means an empty one: a caller sizing
        // a pool by the count must not get fewer workers than
        // it asked for
        assert_eq!(group_sizes(Some(3), None, false, 7).unwrap(), [3, 2, 2]);
        assert_eq!(group_sizes(Some(3), None, false, 2).unwrap(), [1, 1, 0]);
        assert_eq!(group_sizes(Some(2), None, false, 0).unwrap(), [0, 0]);
    }

    #[test]
    fn splitting_by_size_yields_no_file_for_an_empty_input() {
        // the other way round from -k: a chunk size says how
        // big a file is, not how many there are
        assert_eq!(group_sizes(None, Some(2), false, 7).unwrap(), [2, 2, 2, 1]);
        assert_eq!(
            group_sizes(None, Some(2), false, 0).unwrap(),
            [] as [usize; 0]
        );
        assert_eq!(group_sizes(None, Some(10), false, 7).unwrap(), [7]);
    }

    #[test]
    fn each_is_one_file_per_record() {
        assert_eq!(group_sizes(None, None, true, 3).unwrap(), [1, 1, 1]);
    }

    #[test]
    fn every_group_size_sums_back_to_the_record_count() {
        // the invariant that matters: split must not lose or
        // repeat a record, whichever way it was asked for
        for len in 0..20usize {
            for k in 1..6 {
                let by_parts: usize = group_sizes(Some(k), None, false, len).unwrap().iter().sum();
                let by_size: usize = group_sizes(None, Some(k), false, len).unwrap().iter().sum();

                assert_eq!(by_parts, len, "-k {k} over {len}");
                assert_eq!(by_size, len, "-n {k} over {len}");
            }
        }
    }

    #[test]
    fn a_zero_sized_split_is_refused_rather_than_looping() {
        assert!(group_sizes(Some(0), None, false, 5).is_err());
        assert!(group_sizes(None, Some(0), false, 5).is_err());
        assert!(group_sizes(None, None, false, 5).is_err());
    }

    #[test]
    fn a_part_is_numbered_so_the_files_sort_into_record_order() {
        // zero-padded, so part 10 does not sort before part 2
        assert_eq!(
            part_path(Path::new("out/proteins"), 2, "fa"),
            PathBuf::from("out/proteins.0002.fa")
        );
    }
}

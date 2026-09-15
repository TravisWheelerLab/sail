use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Args;
use std::io::Write;

use libsail::collection::Indexable;
use libsail::index::Reader;

use crate::cli::{FormatArg, ReadArgs};
use crate::input::{Backend, Inputs, Needs, dispatch, indexed, indexed_path};
use crate::output::{emit, write_framed, writer};

#[derive(Args)]
pub struct GetArgs {
    /// file to read, or - for stdin
    pub input: PathBuf,

    /// positions to extract, 1-based and inclusive: 3, or 2-4
    #[arg(required = true)]
    pub positions: Vec<String>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    #[command(flatten)]
    pub read: ReadArgs,
}

impl GetArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::plan(
            std::slice::from_ref(&self.input),
            self.format,
            self.read,
            Needs::Sparse,
        )?;
        let entry = &inputs.entries()[0];
        let out = writer(self.output.as_deref())?;

        let wanted = self
            .positions
            .iter()
            .map(|spec| parse_span(spec))
            .collect::<Result<Vec<_>>>()?;

        match inputs.backend() {
            Backend::Indexed => {
                let path = indexed_path(entry)?;
                let mut out = out;

                indexed!(inputs.format(), path, |collection| {
                    let held = collection.index().len();

                    for &(_first, last) in &wanted {
                        if last > held {
                            bail!(
                                "position {last} was asked for and the input holds {held} records"
                            )
                        }
                    }

                    // only the positions asked for are read:
                    // `get 2000000` reads one record rather
                    // than two million
                    for &(first, last) in &wanted {
                        for at in first..=last {
                            let record = collection.record(at - 1)?.expect("a counted record");

                            write_framed(inputs.format(), &record, &mut out)?;
                        }
                    }
                });
                out.flush()?;

                Ok(())
            }
            Backend::Stream => {
                // held rather than written as they are found,
                // so a position past the end is reported the
                // way the in-memory path reports it: nothing
                // written, then the error
                let last = wanted.iter().map(|&(_, last)| last).max().unwrap_or(0);
                let mut held: Vec<Option<Vec<u8>>> = vec![None; last];

                let mut reader = Reader::new(entry.reader()?, inputs.format());
                let mut n = 0;

                while n < last && reader.advance()? {
                    // only the positions asked for are kept,
                    // so this is O(what was asked) and not
                    // O(how far into the file it sits)
                    if wanted.iter().any(|&(f, l)| n + 1 >= f && n < l) {
                        held[n] = Some(reader.record().to_vec());
                    }

                    n += 1;
                }

                if n < last {
                    bail!("position {last} was asked for and the input holds {n} records")
                }

                let mut out = out;
                for &(first, last) in &wanted {
                    for at in first..=last {
                        let record = held[at - 1]
                            .as_ref()
                            .expect("a position inside a span was kept");

                        write_framed(inputs.format(), record, &mut out)?;
                    }
                }
                out.flush()?;

                Ok(())
            }
            Backend::Memory => {
                dispatch!(inputs.format(), entry, |collection, _size, _name, write| {
                    let picked = pick(&collection, &wanted)?;

                    emit(picked, write, out)
                })
            }
        }
    }
}

// ---

/// The records at every position in `wanted`, in the order asked for.
fn pick<C: Indexable>(collection: &C, wanted: &[(usize, usize)]) -> Result<Vec<C::Record>> {
    let mut picked = Vec::new();

    for &(first, last) in wanted {
        // named before the loop, because a range running past
        // the end must say so rather than come back short
        if last > collection.len() {
            bail!(
                "position {last} was asked for and the input holds {} records",
                collection.len()
            )
        }

        for n in first..=last {
            picked.push(
                collection
                    .cloned(n - 1)
                    .expect("a position below len() is a record"),
            );
        }
    }

    Ok(picked)
}

/// One position or one inclusive range, both 1-based, as a `(first, last)`
/// pair.
//
// 1-based because the domain's own tools count records from
// 1, and inclusive at both ends for the same reason -- a
// caller writing 2-4 means three records
fn parse_span(spec: &str) -> Result<(usize, usize)> {
    let number = |s: &str| -> Result<usize> {
        let n: usize = s
            .parse()
            .with_context(|| format!("{spec:?} is a position or a range like 2-4"))?;

        if n == 0 {
            bail!("positions count from 1, and {spec:?} names 0")
        }

        Ok(n)
    };

    let (first, last) = match spec.split_once('-') {
        None => {
            let n = number(spec)?;
            (n, n)
        }
        Some((a, b)) => (number(a)?, number(b)?),
    };

    if first > last {
        bail!("{spec:?} runs backwards: {first} is past {last}")
    }

    Ok((first, last))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_bare_number_is_the_one_record_at_that_position() {
        assert_eq!(parse_span("3").unwrap(), (3, 3));
    }

    #[test]
    fn a_range_counts_from_one_and_includes_both_ends() {
        // 2-4 is three records, not two: the domain's tools
        // number records from 1 and count the last one in
        let (first, last) = parse_span("2-4").unwrap();

        assert_eq!((first, last), (2, 4));
        assert_eq!(last - first + 1, 3);
    }

    #[test]
    fn position_zero_is_refused_rather_than_read_as_the_first_record() {
        // an off-by-one here hands back the wrong record in
        // silence, so it has to be an error
        assert!(parse_span("0").is_err());
        assert!(parse_span("0-2").is_err());
    }

    #[test]
    fn a_backwards_range_is_refused_rather_than_yielding_nothing() {
        assert!(parse_span("4-2").is_err());
    }

    #[test]
    fn something_that_is_not_a_position_names_itself_in_the_error() {
        let error = parse_span("two").unwrap_err().to_string();

        assert!(error.contains("two"), "{error}");
    }
}

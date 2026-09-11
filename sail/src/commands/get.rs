use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Args;
use libsail::collection::Indexable;

use crate::cli::FormatArg;
use crate::input::{Inputs, dispatch};
use crate::output::{emit, writer};

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
}

impl GetArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::read(std::slice::from_ref(&self.input), self.format)?;
        let entry = &inputs.entries()[0];
        let out = writer(self.output.as_deref())?;

        let wanted = self
            .positions
            .iter()
            .map(|spec| parse_span(spec))
            .collect::<Result<Vec<_>>>()?;

        dispatch!(inputs.format(), entry, |collection| {
            let picked = pick(&collection, &wanted)?;

            emit(picked, out)
        })
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

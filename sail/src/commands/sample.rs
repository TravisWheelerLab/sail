use std::path::PathBuf;

use anyhow::{Result, bail};
use clap::Args;
use libsail::collection::{Indexable, Iterable};

use crate::cli::FormatArg;
use crate::input::{Inputs, dispatch};
use crate::output::{emit, writer};

#[derive(Args)]
pub struct SampleArgs {
    /// file to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: PathBuf,

    /// how many records to draw
    #[arg(short = 'n', long, conflicts_with = "fraction")]
    pub count: Option<usize>,

    /// what share of the records to draw, between 0 and 1
    #[arg(short = 'p', long)]
    pub fraction: Option<f64>,

    /// seed the draw, so the same input gives the same sample
    #[arg(long)]
    pub seed: Option<u64>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,
}

impl SampleArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::read(std::slice::from_ref(&self.input), self.format)?;
        let entry = &inputs.entries()[0];
        let out = writer(self.output.as_deref())?;
        let (count, fraction) = (self.count, self.fraction);

        // an unseeded run picks a seed rather than taking a
        // second code path: a draw that ends up in a paper
        // has to reproduce, and libsail has one draw
        let seed = self.seed.unwrap_or_else(rand::random);

        dispatch!(inputs.format(), entry, |collection, _size, _name, write| {
            let n = draw_size(count, fraction, collection.len())?;

            emit(collection.sample(n, seed).iter(), write, out)
        })
    }
}

// ---

/// How many records to draw, from `-n` or from `-p`.
fn draw_size(count: Option<usize>, fraction: Option<f64>, len: usize) -> Result<usize> {
    match (count, fraction) {
        (Some(n), _) => {
            // sampling without replacement cannot return more
            // records than there are, and the rand routine
            // panics rather than clamping -- so the error is
            // raised here, naming both numbers
            if n > len {
                bail!("{n} records were asked for and the input holds {len}")
            }

            Ok(n)
        }

        (None, Some(p)) => {
            if !(0.0..=1.0).contains(&p) {
                bail!("a fraction runs from 0 to 1, and {p} does not")
            }

            // round rather than truncate, so -p 0.5 over 5
            // records draws 3 and not 2
            Ok((len as f64 * p).round() as usize)
        }

        (None, None) => bail!("sample needs -n or -p"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn drawing_more_records_than_there_are_is_an_error_and_not_a_clamp() {
        // the rand routine panics rather than clamping, so
        // this has to be caught before it is reached
        assert!(draw_size(Some(6), None, 5).is_err());
        assert_eq!(draw_size(Some(5), None, 5).unwrap(), 5);
    }

    #[test]
    fn a_fraction_rounds_rather_than_truncating() {
        // half of five is 2.5, and truncating would quietly
        // draw a smaller sample than was asked for
        assert_eq!(draw_size(None, Some(0.5), 5).unwrap(), 3);
        assert_eq!(draw_size(None, Some(1.0), 5).unwrap(), 5);
        assert_eq!(draw_size(None, Some(0.0), 5).unwrap(), 0);
    }

    #[test]
    fn a_fraction_outside_zero_to_one_is_refused() {
        assert!(draw_size(None, Some(1.5), 5).is_err());
        assert!(draw_size(None, Some(-0.5), 5).is_err());
    }

    #[test]
    fn neither_n_nor_p_is_an_error_rather_than_an_empty_sample() {
        assert!(draw_size(None, None, 5).is_err());
    }
}

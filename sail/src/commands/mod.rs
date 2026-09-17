pub mod cat;
pub mod count;
pub mod dedup;
pub mod fetch;
pub mod filter;
pub mod get;
pub mod grep;
pub mod head;
pub mod index;
pub mod names;
pub mod reformat;
pub mod rename;
pub mod sample;
pub mod shuffle;
pub mod skim;
pub mod sort;
pub mod split;
pub mod stats;
pub mod tail;
pub mod validate;

use anyhow::{Result, bail};

/// How many records to draw, from `-n` or from `-p`.
pub(crate) fn draw_size(count: Option<usize>, fraction: Option<f64>, len: usize) -> Result<usize> {
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

        (None, None) => bail!("how many records to draw: give -n or -p"),
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

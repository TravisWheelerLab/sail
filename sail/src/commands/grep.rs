use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Args;
use libsail::collection::Iterable;
use regex::bytes::Regex;

use crate::cli::FormatArg;
use crate::input::{Inputs, dispatch};
use crate::output::{emit, writer};

#[derive(Args)]
pub struct GrepArgs {
    /// pattern to match against each record's name
    pub pattern: String,

    /// files to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: Vec<PathBuf>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// keep the records that do not match
    #[arg(short = 'v', long)]
    pub invert: bool,

    /// match the whole name rather than any part of it
    #[arg(short = 'x', long)]
    pub exact: bool,

    /// ignore case
    #[arg(short = 'i', long)]
    pub ignore_case: bool,
}

impl GrepArgs {
    pub fn run(self) -> Result<()> {
        let inputs = Inputs::read(&self.input, self.format)?;
        let pattern = self.pattern()?;
        let invert = self.invert;
        let mut out = writer(self.output.as_deref())?;

        for entry in inputs.entries() {
            dispatch!(inputs.format(), entry, |collection, _size, name| {
                let kept = collection
                    .iter()
                    // a record with no name matches nothing,
                    // so inverting keeps it: it is not one of
                    // the records the pattern found
                    .filter(|r| name(r).is_some_and(|n| pattern.is_match(n)) != invert);

                emit(kept, &mut out)?;
            });
        }

        Ok(())
    }

    /// The pattern, on bytes because a FASTA name is not required to be utf-8.
    fn pattern(&self) -> Result<Regex> {
        let mut source = if self.exact {
            format!("^(?:{})$", self.pattern)
        } else {
            self.pattern.clone()
        };

        if self.ignore_case {
            source.insert_str(0, "(?i)");
        }

        Regex::new(&source).with_context(|| format!("the pattern {:?}", self.pattern))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(pattern: &str) -> GrepArgs {
        GrepArgs {
            pattern: pattern.to_string(),
            input: Vec::new(),
            format: None,
            output: None,
            invert: false,
            exact: false,
            ignore_case: false,
        }
    }

    #[test]
    fn a_pattern_matches_any_part_of_a_name_by_default() {
        assert!(args("HUMAN").pattern().unwrap().is_match(b"DLG4_HUMAN"));
    }

    #[test]
    fn exact_anchors_the_whole_name_rather_than_each_alternative() {
        // grouping matters: ^a|b$ anchors only the branch it
        // sits beside, so "a" would match inside "xay"
        let exact = GrepArgs {
            exact: true,
            ..args("PDZ|SH3")
        };
        let pattern = exact.pattern().unwrap();

        assert!(pattern.is_match(b"PDZ"));
        assert!(pattern.is_match(b"SH3"));
        assert!(!pattern.is_match(b"PDZ_1"));
        assert!(!pattern.is_match(b"xSH3"));
    }

    #[test]
    fn a_name_that_is_not_utf_8_is_still_matchable() {
        // a FASTA name carries whatever the file carried.
        // decoding it to a String here would leave a latin-1
        // name unmatchable by any pattern at all
        //
        // note: the pattern stays unicode-aware, so `.` does
        //       not match the stray byte -- a literal does,
        //       which is what a name pattern is in practice
        assert!(args("c").pattern().unwrap().is_match(b"a\xffc"));
        assert!(!args("a.c").pattern().unwrap().is_match(b"a\xffc"));
    }

    #[test]
    fn a_pattern_that_is_not_a_regex_names_itself_in_the_error() {
        let error = args("[unclosed").pattern().unwrap_err().to_string();

        assert!(error.contains("[unclosed"), "{error}");
    }

    #[test]
    fn ignore_case_folds_both_sides() {
        let folded = GrepArgs {
            ignore_case: true,
            ..args("pdz")
        };

        assert!(folded.pattern().unwrap().is_match(b"PDZ"));
    }
}

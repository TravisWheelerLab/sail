use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Args;
use libsail::collection::Indexable;
use libsail::format::Format;
use regex::bytes::Regex;

use crate::cli::FormatArg;
use crate::input::Inputs;
use crate::output::{emit, writer};

#[derive(Args)]
pub struct RenameArgs {
    /// file to read, or - for stdin
    #[arg(default_value = "-")]
    pub input: PathBuf,

    /// put this in front of every name
    #[arg(long)]
    pub prefix: Option<String>,

    /// put this after every name
    #[arg(long)]
    pub suffix: Option<String>,

    /// replace what this pattern matches, with --with
    #[arg(long, requires = "with")]
    pub replace: Option<String>,

    /// what --replace puts in place of a match; ${1} is the first group
    #[arg(long, requires = "replace")]
    pub with: Option<String>,

    /// number the records, replacing each name with this stem plus its position
    #[arg(long)]
    pub number: Option<String>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// where to write [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,
}

impl RenameArgs {
    pub fn run(self) -> Result<()> {
        let rules = self.rules()?;
        let inputs = Inputs::read(std::slice::from_ref(&self.input), self.format)?;
        let entry = &inputs.entries()[0];
        let out = writer(self.output.as_deref())?;

        // not through dispatch!: renaming writes to a record
        // rather than reading one, and where a name is stored
        // is the one thing the three formats do not share --
        // a Vec<u8> field, a #=GF line, and a header string
        match inputs.format() {
            Format::Fasta => {
                let collection = libsail::seq::fasta::Fasta::new(entry.reader()?)?;
                let renamed = (0..collection.len())
                    .map(|n| {
                        let mut record = collection.cloned(n).expect("a counted record");
                        record.name = rules.apply(&record.name, n);

                        record
                    })
                    .collect::<Vec<_>>();

                emit(renamed, crate::output::write::fasta, out)
            }

            Format::Stockholm => {
                let collection = libsail::seq::stockholm::Stockholm::new(entry.reader()?)?;
                let renamed = (0..collection.len())
                    .map(|n| {
                        let mut record = collection.cloned(n).expect("a counted record");
                        let was = record.id().unwrap_or_default().to_vec();
                        let now = rules.apply(&was, n);

                        set_id(&mut record, &now)?;

                        Ok(record)
                    })
                    .collect::<Result<Vec<_>>>()?;

                emit(renamed, crate::output::write::stockholm, out)
            }

            Format::Hmm => {
                let collection = libsail::seq::p7hmm::Hmm::new(entry.reader()?)?;
                let renamed = (0..collection.len())
                    .map(|n| {
                        let mut record = collection.cloned(n).expect("a counted record");
                        record.header.name = rules.apply(&record.header.name, n);

                        Ok(record)
                    })
                    .collect::<Result<Vec<_>>>()?;

                emit(renamed, crate::output::write::hmm, out)
            }
        }
    }

    fn rules(&self) -> Result<Rules> {
        if self.prefix.is_none()
            && self.suffix.is_none()
            && self.replace.is_none()
            && self.number.is_none()
        {
            bail!("rename needs --prefix, --suffix, --replace or --number")
        }

        let replace = match (&self.replace, &self.with) {
            (Some(pattern), Some(with)) => Some((
                Regex::new(pattern).with_context(|| format!("the pattern {pattern:?}"))?,
                with.as_bytes().to_vec(),
            )),
            _ => None,
        };

        Ok(Rules {
            prefix: self.prefix.clone().unwrap_or_default().into_bytes(),
            suffix: self.suffix.clone().unwrap_or_default().into_bytes(),
            replace,
            number: self.number.clone().map(String::into_bytes),
        })
    }
}

// ---

/// What each rewrite does to one name.
struct Rules {
    prefix: Vec<u8>,
    suffix: Vec<u8>,
    replace: Option<(Regex, Vec<u8>)>,
    number: Option<Vec<u8>>,
}

impl Rules {
    /// The new name for the record at position `n`.
    fn apply(&self, name: &[u8], n: usize) -> Vec<u8> {
        // the order is fixed and is what makes the flags
        // combinable: --number replaces the name outright,
        // then --replace rewrites it, then --prefix and
        // --suffix wrap what is left
        let mut now = match &self.number {
            // 1-based, and zero-padded to four so the names
            // sort into record order
            Some(stem) => {
                let mut numbered = stem.clone();
                numbered.extend_from_slice(format!("{:04}", n + 1).as_bytes());

                numbered
            }
            None => name.to_vec(),
        };

        if let Some((pattern, with)) = &self.replace {
            now = pattern.replace_all(&now, &with[..]).into_owned();
        }

        // after --replace, so a pattern cannot match the
        // prefix this same call just added
        let mut out = self.prefix.clone();
        out.extend_from_slice(&now);
        out.extend_from_slice(&self.suffix);

        out
    }
}

/// Rewrite an alignment's `#=GF ID`, adding one where it had none.
fn set_id(record: &mut libsail::seq::stockholm::StockholmRecord, name: &[u8]) -> Result<()> {
    let name = name.to_vec();

    // the record reads its id back out of gf rather than
    // storing it, so renaming means editing that line --
    // nothing else in the record would be written out
    match record.gf.iter_mut().find(|(feature, _)| feature == b"ID") {
        Some((_, value)) => *value = name,

        // in front, because Pfam writes ID first and a reader
        // scanning for it should not have to pass the rest
        None => record.gf.insert(0, (b"ID".to_vec(), name)),
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rules(args: RenameArgs) -> Rules {
        args.rules().unwrap()
    }

    fn args() -> RenameArgs {
        RenameArgs {
            input: PathBuf::from("-"),
            prefix: None,
            suffix: None,
            replace: None,
            with: None,
            number: None,
            format: None,
            output: None,
        }
    }

    #[test]
    fn a_prefix_and_a_suffix_wrap_the_name_they_are_given() {
        let r = rules(RenameArgs {
            prefix: Some("x_".into()),
            suffix: Some("_y".into()),
            ..args()
        });

        assert_eq!(r.apply(b"AB", 0), b"x_AB_y");
    }

    #[test]
    fn numbering_counts_from_one_and_pads_so_the_names_sort() {
        // seq10 sorting before seq2 is the whole reason for
        // the padding
        let r = rules(RenameArgs {
            number: Some("seq".into()),
            ..args()
        });

        assert_eq!(r.apply(b"anything", 0), b"seq0001");
        assert_eq!(r.apply(b"anything", 9), b"seq0010");
    }

    #[test]
    fn a_replacement_runs_before_the_prefix_it_would_otherwise_match() {
        // order is what makes the flags combinable: a pattern
        // must not match the prefix this same call just added
        let r = rules(RenameArgs {
            prefix: Some("HUMAN_".into()),
            replace: Some("HUMAN".into()),
            with: Some("H".into()),
            ..args()
        });

        assert_eq!(r.apply(b"DLG4_HUMAN", 0), b"HUMAN_DLG4_H");
    }

    #[test]
    fn a_replacement_reaches_a_captured_group_through_the_braced_form() {
        let braced = rules(RenameArgs {
            replace: Some("(.+)_(.+)".into()),
            with: Some("${2}_${1}".into()),
            ..args()
        });

        assert_eq!(braced.apply(b"DLG4_HUMAN", 0), b"HUMAN_DLG4");
    }

    #[test]
    fn a_bare_group_reference_swallows_an_underscore_after_it() {
        // **note: the trap this records is the regex crate's,
        //         not this command's: $2_ reads as a group
        //         *named* "2_", which does not exist and
        //         expands to nothing. the name comes out
        //         wrong in silence, so --with documents the
        //         braced form and this pins the behaviour
        //         against a future regex release changing it
        let bare = rules(RenameArgs {
            replace: Some("(.+)_(.+)".into()),
            with: Some("$2_$1".into()),
            ..args()
        });

        assert_eq!(bare.apply(b"DLG4_HUMAN", 0), b"DLG4");
    }

    #[test]
    fn renaming_with_no_rule_at_all_is_refused() {
        // silently emitting the input unchanged would look
        // like the rename worked
        assert!(args().rules().is_err());
    }
}

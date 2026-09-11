use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};
use libsail::format::Format;

use crate::commands::index::IndexArgs;

#[derive(Parser)]
#[command(name = "sail", version, about = "sequence analysis I/O tool")]
pub struct Cli {
    #[command(subcommand)]
    command: Command,
}

impl Cli {
    pub fn run(self) -> Result<()> {
        match self.command {
            Command::Index(args) => args.run(),
        }
    }
}

#[derive(Subcommand)]
enum Command {
    /// Write a record index beside a file
    Index(IndexArgs),
}

/// A format as named on the command line.
//
// libsail's Format with a ValueEnum impl, kept here
// because libsail must not depend on clap
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum FormatArg {
    Fasta,
    Stockholm,
    Hmm,
}

impl From<FormatArg> for Format {
    fn from(arg: FormatArg) -> Format {
        match arg {
            FormatArg::Fasta => Format::Fasta,
            FormatArg::Stockholm => Format::Stockholm,
            FormatArg::Hmm => Format::Hmm,
        }
    }
}

use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};
use libsail::format::Format;

use crate::commands::count::CountArgs;
use crate::commands::dedup::DedupArgs;
use crate::commands::fetch::FetchArgs;
use crate::commands::filter::FilterArgs;
use crate::commands::get::GetArgs;
use crate::commands::grep::GrepArgs;
use crate::commands::head::HeadArgs;
use crate::commands::index::IndexArgs;
use crate::commands::names::NamesArgs;
use crate::commands::stats::StatsArgs;
use crate::commands::tail::TailArgs;
use crate::commands::validate::ValidateArgs;

#[derive(Parser)]
#[command(name = "sail", version, about = "sequence analysis I/O tool")]
pub struct Cli {
    #[command(subcommand)]
    command: Command,
}

impl Cli {
    pub fn run(self) -> Result<()> {
        match self.command {
            Command::Count(args) => args.run(),
            Command::Dedup(args) => args.run(),
            Command::Fetch(args) => args.run(),
            Command::Filter(args) => args.run(),
            Command::Get(args) => args.run(),
            Command::Grep(args) => args.run(),
            Command::Head(args) => args.run(),
            Command::Index(args) => args.run(),
            Command::Names(args) => args.run(),
            Command::Stats(args) => args.run(),
            Command::Tail(args) => args.run(),
            Command::Validate(args) => args.run(),
        }
    }
}

#[derive(Subcommand)]
enum Command {
    /// Count the records in a file
    Count(CountArgs),

    /// Drop records with a repeated name, keeping the first
    Dedup(DedupArgs),

    /// Extract records by name
    Fetch(FetchArgs),

    /// Keep records by size
    Filter(FilterArgs),

    /// Extract records by position, 1-based
    Get(GetArgs),

    /// Keep records whose name matches a pattern
    Grep(GrepArgs),

    /// Keep the first n records
    Head(HeadArgs),

    /// Write a record index beside a file
    Index(IndexArgs),

    /// Print every record's name, one per line
    Names(NamesArgs),

    /// Report the record count and the distribution of their sizes
    Stats(StatsArgs),

    /// Keep the last n records
    Tail(TailArgs),

    /// Parse every record and report the first failure
    Validate(ValidateArgs),
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

/// Which of an alignment's two axes `--by` names.
//
// only Stockholm has two: a FASTA record's size is its
// residue count and a profile's is its LENG, and neither
// has a second axis to choose between
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Axis {
    /// number of aligned rows
    Depth,
    /// number of alignment columns
    Width,
}

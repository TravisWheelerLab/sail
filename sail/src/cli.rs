use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};
use libsail::format::Format;

use crate::commands::cat::CatArgs;
use crate::commands::count::CountArgs;
use crate::commands::dedup::DedupArgs;
use crate::commands::fetch::FetchArgs;
use crate::commands::filter::FilterArgs;
use crate::commands::get::GetArgs;
use crate::commands::grep::GrepArgs;
use crate::commands::head::HeadArgs;
use crate::commands::index::IndexArgs;
use crate::commands::names::NamesArgs;
use crate::commands::reformat::ReformatArgs;
use crate::commands::rename::RenameArgs;
use crate::commands::sample::SampleArgs;
use crate::commands::shuffle::ShuffleArgs;
use crate::commands::sort::SortArgs;
use crate::commands::split::SplitArgs;
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
            Command::Cat(args) => args.run(),
            Command::Count(args) => args.run(),
            Command::Dedup(args) => args.run(),
            Command::Fetch(args) => args.run(),
            Command::Filter(args) => args.run(),
            Command::Get(args) => args.run(),
            Command::Grep(args) => args.run(),
            Command::Head(args) => args.run(),
            Command::Index(args) => args.run(),
            Command::Names(args) => args.run(),
            Command::Reformat(args) => args.run(),
            Command::Rename(args) => args.run(),
            Command::Sample(args) => args.run(),
            Command::Shuffle(args) => args.run(),
            Command::Sort(args) => args.run(),
            Command::Split(args) => args.run(),
            Command::Stats(args) => args.run(),
            Command::Tail(args) => args.run(),
            Command::Validate(args) => args.run(),
        }
    }
}

#[derive(Subcommand)]
enum Command {
    /// Concatenate several files into one stream
    Cat(CatArgs),

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

    /// Re-emit every record through this crate's writer
    Reformat(ReformatArgs),

    /// Rewrite every record name
    Rename(RenameArgs),

    /// Draw records at random
    Sample(SampleArgs),

    /// Put the records in a random order
    Shuffle(ShuffleArgs),

    /// Order records by name or by size
    Sort(SortArgs),

    /// Write the records out across several files
    Split(SplitArgs),

    /// Report the record count and the distribution of their sizes
    Stats(StatsArgs),

    /// Keep the last n records
    Tail(TailArgs),

    /// Parse every record and report the first failure
    Validate(ValidateArgs),
}

/// A format as named on the command line.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum FormatArg {
    Fasta,
    Stockholm,
    Hmm,
}

impl From<FormatArg> for Format {
    fn from(arg: FormatArg) -> Format {
        // a mapping onto an identical enum because libsail
        // must not depend on clap: FormatArg is libsail's
        // Format with a ValueEnum impl
        match arg {
            FormatArg::Fasta => Format::Fasta,
            FormatArg::Stockholm => Format::Stockholm,
            FormatArg::Hmm => Format::Hmm,
        }
    }
}

/// How an operation reads its input.
#[derive(clap::Args, Clone, Copy, Debug, Default)]
#[command(next_help_heading = "Input")]
pub struct ReadArgs {
    /// stream the input, hold it in memory, or address it on disk
    //
    // one value enum rather than --stream/--indexed
    // booleans: booleans need a conflicts_with for every
    // pair, and their help never says which one is the
    // default
    #[arg(long = "read", value_enum, default_value_t, value_name = "MODE")]
    pub mode: Mode,
}

/// A read backend as named on the command line.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum Mode {
    /// let the operation pick: streaming unless it needs a record count or a
    /// few records by name, and never indexed on stdin
    //
    // the one Mode that is not an input::Backend, the same
    // split FormatArg makes against libsail's Format
    #[default]
    Auto,
    /// one forward pass, one record at a time, in a reused buffer
    Stream,
    /// parse every record into memory before doing anything
    Memory,
    /// keep a byte offset per record and read one on demand
    Indexed,
}

/// What order `fetch` writes the records it found.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum Order {
    /// the order the input holds them
    #[default]
    File,
    /// the order the names were asked for
    Asked,
}

/// Which of an alignment's two axes `--by` names.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Axis {
    /// number of aligned rows
    Depth,
    /// number of alignment columns
    Width,
}

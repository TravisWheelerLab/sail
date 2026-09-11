use std::process::ExitCode;

use clap::Parser;

mod cli;
mod commands;
mod input;
mod output;

fn main() -> ExitCode {
    if let Err(e) = cli::Cli::parse().run() {
        // {:#} so the context a command attached is
        // reported alongside what actually failed
        eprintln!("sail: {e:#}");

        return ExitCode::FAILURE;
    }

    ExitCode::SUCCESS
}

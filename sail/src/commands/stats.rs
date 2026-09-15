use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use libsail::collection::Iterable;
use libsail::index::Reader;

use crate::cli::{Axis, FormatArg, ReadArgs};
use crate::input::{
    Backend, Inputs, Needs, axis_for, dispatch, indexed, indexed_path, size_framed,
};

#[derive(Args)]
pub struct StatsArgs {
    /// files to summarise, or - for stdin
    #[arg(default_value = "-")]
    pub input: Vec<PathBuf>,

    /// assert the input is this format, and fail if it is not
    #[arg(long, value_enum)]
    pub format: Option<FormatArg>,

    /// for an alignment, measure rows or columns [default: depth]
    #[arg(long, value_enum)]
    pub by: Option<Axis>,

    #[command(flatten)]
    pub read: ReadArgs,
}

impl StatsArgs {
    pub fn run(self) -> Result<()> {
        for line in self.report()? {
            println!("{line}");
        }

        Ok(())
    }

    /// One `key<TAB>value` line per statistic, over every input together.
    fn report(&self) -> Result<Vec<String>> {
        let inputs = Inputs::plan(&self.input, self.format, self.read, Needs::Pass)?;
        let axis = axis_for(inputs.format(), self.by)?;

        // every size, so the order statistics can be read off
        // a sort rather than estimated. a size is a usize per
        // record against a file already held whole
        let mut sizes = Vec::new();
        for entry in inputs.entries() {
            match inputs.backend() {
                // a usize per record and never the record: a
                // FASTA length is counted off the framed
                // bytes, so the residues are never built
                Backend::Stream => {
                    let mut reader = Reader::new(entry.reader()?, inputs.format());

                    while reader.advance()? {
                        sizes.push(size_framed(inputs.format(), reader.record(), axis)?);
                    }
                }
                Backend::Indexed => {
                    let path = indexed_path(entry)?;

                    indexed!(inputs.format(), path, |collection| {
                        for n in 0..collection.index().len() {
                            let record = collection.record(n)?.expect("a counted record");

                            sizes.push(size_framed(inputs.format(), &record, axis)?);
                        }
                    });
                }
                Backend::Memory => {
                    dispatch!(
                        inputs.format(),
                        entry,
                        axis,
                        |collection, size, _name, _write| sizes.extend(collection.iter().map(size))
                    );
                }
            }
        }

        Ok(summarise(
            &mut sizes,
            inputs.format(),
            axis,
            inputs.entries().len(),
        ))
    }
}

// ---

/// The report for `sizes`, which is sorted in place to read the order
/// statistics off it.
fn summarise(
    sizes: &mut [usize],
    format: libsail::format::Format,
    axis: Axis,
    files: usize,
) -> Vec<String> {
    let mut lines = vec![format!("format\t{format}")];

    if files > 1 {
        lines.push(format!("files\t{files}"));
    }
    if format == libsail::format::Format::Stockholm {
        // an alignment's size is ambiguous, so the report
        // says which axis it measured rather than leaving the
        // reader to recall the default
        lines.push(format!("measured\t{}", axis_name(axis)));
    }

    lines.push(format!("records\t{}", sizes.len()));

    let total: usize = sizes.iter().sum();
    lines.push(format!("total\t{total}"));

    // an empty input has a count and a total and no
    // distribution at all; reporting min as 0 would read as a
    // record of size zero
    if sizes.is_empty() {
        return lines;
    }

    sizes.sort_unstable();

    lines.push(format!("min\t{}", sizes[0]));
    lines.push(format!("q1\t{}", quantile(sizes, 1, 4)));
    lines.push(format!("median\t{}", quantile(sizes, 1, 2)));
    lines.push(format!("q3\t{}", quantile(sizes, 3, 4)));
    lines.push(format!("max\t{}", sizes[sizes.len() - 1]));
    lines.push(format!("mean\t{:.1}", total as f64 / sizes.len() as f64));

    lines
}

/// The value at `num / den` of the way through `sorted`.
fn quantile(sorted: &[usize], num: usize, den: usize) -> usize {
    // the nearest-rank definition: the smallest value at or
    // above the p-th percentile, so the result is one of the
    // input values rather than an interpolation between two
    // of them -- a size counts residues or rows, and no
    // record has 376.5 of either
    let rank = (sorted.len() * num).div_ceil(den).max(1);

    sorted[rank - 1]
}

fn axis_name(axis: Axis) -> &'static str {
    match axis {
        Axis::Depth => "depth",
        Axis::Width => "width",
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use libsail::format::Format;

    use super::*;

    fn fixture(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../fixtures")
            .join(name)
    }

    fn stats(paths: &[PathBuf], by: Option<Axis>) -> Result<Vec<String>> {
        StatsArgs {
            read: ReadArgs::default(),
            input: paths.to_vec(),
            format: None,
            by,
        }
        .report()
    }

    fn value(lines: &[String], key: &str) -> String {
        lines
            .iter()
            .find_map(|l| l.strip_prefix(&format!("{key}\t")))
            .unwrap_or_else(|| panic!("no {key} line in {lines:?}"))
            .to_string()
    }

    #[test]
    fn a_fasta_record_is_measured_in_residues_and_not_in_bytes_of_text() {
        // the five proteins are wrapped over many lines, so a
        // byte count would be far larger than the residue one
        let lines = stats(&[fixture("proteins.fa")], None).unwrap();

        assert_eq!(value(&lines, "format"), "fasta");
        assert_eq!(value(&lines, "records"), "5");
        assert_eq!(value(&lines, "min"), "228");
        assert_eq!(
            value(&lines, "max"),
            "724",
            "PSD-95, the longest of the five"
        );
        assert_eq!(value(&lines, "total"), "2335");
        assert_eq!(value(&lines, "median"), "499");
    }

    #[test]
    fn a_profile_is_measured_by_its_leng_and_not_by_its_node_lines() {
        let lines = stats(&[fixture("models.hmm")], None).unwrap();

        assert_eq!(value(&lines, "format"), "hmm");
        assert_eq!(value(&lines, "records"), "2");
        assert_eq!(value(&lines, "min"), "46", "SH3_1");
        assert_eq!(value(&lines, "max"), "81", "PDZ");
        assert_eq!(value(&lines, "total"), "127");
    }

    #[test]
    fn an_alignment_measures_depth_by_default_and_width_when_asked() {
        // the one format with two answers, and the report has
        // to say which it gave
        let depth = stats(&[fixture("families.sto")], None).unwrap();
        let width = stats(&[fixture("families.sto")], Some(Axis::Width)).unwrap();

        assert_eq!(value(&depth, "measured"), "depth");
        assert_eq!(
            value(&depth, "min"),
            "6",
            "six sequences kept from each seed"
        );
        assert_eq!(value(&depth, "max"), "6");

        assert_eq!(value(&width, "measured"), "width");
        assert_eq!(value(&width, "min"), "56", "SH3_1");
        assert_eq!(value(&width, "max"), "110", "PDZ");
    }

    #[test]
    fn by_is_refused_for_a_format_with_only_one_size() {
        // accepted by clap and rejected here, because the
        // format is sniffed and so is not known at parse time
        assert!(stats(&[fixture("proteins.fa")], Some(Axis::Width)).is_err());
        assert!(stats(&[fixture("models.hmm")], Some(Axis::Depth)).is_err());
    }

    #[test]
    fn several_files_are_summarised_as_one_population() {
        // the same file twice: every count doubles and every
        // order statistic stays put
        let one = stats(&[fixture("proteins.fa")], None).unwrap();
        let two = stats(&[fixture("proteins.fa"), fixture("proteins.fa")], None).unwrap();

        assert_eq!(value(&two, "files"), "2");
        assert_eq!(value(&two, "records"), "10");
        assert_eq!(value(&two, "min"), value(&one, "min"));
        assert_eq!(value(&two, "max"), value(&one, "max"));
        assert_eq!(value(&two, "median"), value(&one, "median"));
    }

    #[test]
    fn the_quantiles_are_input_values_rather_than_interpolations() {
        // a size counts residues or rows, so no record has
        // 376.5 of either. nearest-rank keeps every reported
        // number one that some record actually has
        let sizes = [10, 20, 30, 40];

        assert_eq!(quantile(&sizes, 1, 4), 10);
        assert_eq!(quantile(&sizes, 1, 2), 20);
        assert_eq!(quantile(&sizes, 3, 4), 30);
    }

    #[test]
    fn an_empty_input_has_a_count_and_no_distribution() {
        // min over nothing is not 0, which would read as a
        // record of size zero
        let lines = summarise(&mut [], Format::Fasta, Axis::Depth, 1);

        assert_eq!(value(&lines, "records"), "0");
        assert_eq!(value(&lines, "total"), "0");
        assert!(!lines.iter().any(|l| l.starts_with("min\t")), "{lines:?}");
    }
}

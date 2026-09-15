//! HMMER3 profile HMM records.
//!
//! A record is a header of typed fields and a model body: per-node match and
//! insert emission scores plus state transitions, all as negative natural
//! logarithms of probabilities.

use std::borrow::Cow;
use std::fmt::Display;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use crate::Record;
use crate::collection::Indexable;
use crate::error::{Error, Result, bail_parse};
use crate::format::Format;
use crate::index::{Delimit, Delimiter, Index, build_index};
use crate::parse::Parse;
#[cfg(feature = "mmap")]
use crate::source::MmapSource;
use crate::source::{FileSource, Source};

/// Number of state transitions stored per node: m→m, m→i, m→d, i→m, i→i,
/// d→m, d→d.
pub const N_TRANSITIONS: usize = 7;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Alphabet {
    Amino,
    Dna,
    Rna,
    Other(String),
}

impl Alphabet {
    pub fn as_str(&self) -> &str {
        match self {
            Alphabet::Amino => "amino",
            Alphabet::Dna => "DNA",
            Alphabet::Rna => "RNA",
            Alphabet::Other(s) => s,
        }
    }
}

impl From<&str> for Alphabet {
    fn from(s: &str) -> Self {
        match s.to_ascii_lowercase().as_str() {
            "amino" => Alphabet::Amino,
            "dna" => Alphabet::Dna,
            "rna" => Alphabet::Rna,
            _ => Alphabet::Other(s.to_string()),
        }
    }
}

/// Parameters of a fitted score distribution, used to turn a bit score into a
/// P-value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Gumbel {
    /// Location parameter.
    pub tau: f64,
    /// Scale parameter.
    pub lambda: f64,
}

impl Gumbel {
    /// P-value of a bit score under this distribution.
    pub fn p_value(&self, score: f64) -> f64 {
        (-self.lambda * (score - self.tau)).exp()
    }
}

/// Everything above the model body.
#[derive(Debug, Clone, PartialEq)]
pub struct HmmHeader {
    /// The file magic, e.g. `HMMER3/f [3.3 | Nov 2019]`.
    pub version: String,
    pub name: String,
    pub acc: Option<String>,
    pub desc: Option<String>,
    /// Number of match states, the query axis of a search.
    pub leng: usize,
    pub maxl: Option<usize>,
    pub alph: Alphabet,
    pub rf: bool,
    pub mm: bool,
    pub cons: bool,
    pub cs: bool,
    pub map: bool,
    pub date: Option<String>,
    pub com: Vec<String>,
    pub nseq: Option<usize>,
    pub effn: Option<f64>,
    pub cksum: Option<u32>,
    /// Gathering threshold: the curated score above which a hit is a member.
    pub ga: Option<(f32, f32)>,
    /// Trusted cutoff.
    pub tc: Option<(f32, f32)>,
    /// Noise cutoff.
    pub nc: Option<(f32, f32)>,
    pub msv: Option<Gumbel>,
    pub viterbi: Option<Gumbel>,
    pub forward: Option<Gumbel>,
    /// Header lines this parser does not model, in order, so a model written
    /// back keeps them.
    pub extra: Vec<(String, String)>,
}

impl Default for HmmHeader {
    fn default() -> Self {
        HmmHeader {
            version: "HMMER3/f".to_string(),
            name: String::new(),
            acc: None,
            desc: None,
            leng: 0,
            maxl: None,
            alph: Alphabet::Amino,
            rf: false,
            mm: false,
            cons: false,
            cs: false,
            map: false,
            date: None,
            com: Vec::new(),
            nseq: None,
            effn: None,
            cksum: None,
            ga: None,
            tc: None,
            nc: None,
            msv: None,
            viterbi: None,
            forward: None,
            extra: Vec::new(),
        }
    }
}

/// The trailing tokens of a match line (MAP index, consensus residue, RF, MM,
/// CS), verbatim because which are present depends on the header flags.
pub type NodeAnnotation = Vec<String>;

/// The model body: emission and transition scores as negative natural log
/// probabilities, where `f32::INFINITY` is the file's `*`.
#[derive(Debug, Clone, PartialEq)]
pub struct HmmModel {
    /// Alphabet symbols from the `HMM` line.
    pub symbols: Vec<char>,
    /// The transition-header line, kept verbatim for round-tripping.
    pub transition_header: String,
    /// Average residue composition, if the model carries a `COMPO` line.
    pub compo: Option<Vec<f32>>,
    /// Insert emissions for the begin node.
    pub begin_insert: Vec<f32>,
    /// Transitions out of the begin node.
    pub begin_transitions: Vec<f32>,
    /// Match emissions, `leng` rows of `k` scores.
    pub match_emissions: Vec<f32>,
    /// Insert emissions, `leng` rows of `k` scores.
    pub insert_emissions: Vec<f32>,
    /// Transitions, `leng` rows of [`N_TRANSITIONS`] scores.
    pub transitions: Vec<f32>,
    /// Trailing annotations of each match line, `leng` entries.
    pub annotations: Vec<NodeAnnotation>,
}

impl HmmModel {
    /// Alphabet size.
    pub fn k(&self) -> usize {
        self.symbols.len()
    }

    /// Match emission scores for node `i`, counting from 1 as the file does.
    pub fn match_row(&self, i: usize) -> Option<&[f32]> {
        row(&self.match_emissions, i, self.k())
    }

    /// Insert emission scores for node `i`.
    pub fn insert_row(&self, i: usize) -> Option<&[f32]> {
        row(&self.insert_emissions, i, self.k())
    }

    /// Transition scores out of node `i`.
    pub fn transition_row(&self, i: usize) -> Option<&[f32]> {
        row(&self.transitions, i, N_TRANSITIONS)
    }
}

/// Row `i` of a flat table of `width`-wide rows, counting from 1.
fn row(flat: &[f32], i: usize, width: usize) -> Option<&[f32]> {
    // every step is checked: node 0 is the begin node and
    // lives in its own fields, and an absurd index must
    // return None rather than overflow into a
    // plausible-looking slice
    let i = i.checked_sub(1)?;
    let start = i.checked_mul(width)?;
    let end = start.checked_add(width)?;
    flat.get(start..end)
}

/// A profile HMM: header plus model body.
#[derive(Debug, Clone, PartialEq)]
pub struct HmmRecord {
    pub header: HmmHeader,
    pub model: HmmModel,
}

fn fmt_score(v: f32) -> String {
    if v.is_infinite() {
        "*".to_string()
    } else {
        format!("{v:.5}")
    }
}

impl HmmHeader {
    /// The gathering threshold as a P-value under the forward distribution.
    pub fn ga_p_value(&self) -> Option<f64> {
        let ga = self.ga?.0;

        // computed from GA and STATS LOCAL FORWARD rather
        // than stored beside them, which would let the three
        // drift apart
        Some(self.forward?.p_value(ga as f64))
    }
}

impl HmmRecord {
    pub fn write_to<W: Write>(&self, mut w: W) -> Result<()> {
        self.header.write_to(&mut w)?;

        let m = &self.model;
        write!(w, "HMM")?;
        for s in &m.symbols {
            write!(w, "  {s:>6}")?;
        }
        writeln!(w)?;
        writeln!(w, "{}", m.transition_header)?;

        let row = |w: &mut W, vals: &[f32]| -> Result<()> {
            for v in vals {
                write!(w, " {:>8}", fmt_score(*v))?;
            }
            Ok(())
        };

        if let Some(compo) = &m.compo {
            write!(w, "  COMPO")?;
            row(&mut w, compo)?;
            writeln!(w)?;
        }

        write!(w, "       ")?;
        row(&mut w, &m.begin_insert)?;
        writeln!(w)?;
        write!(w, "       ")?;
        row(&mut w, &m.begin_transitions)?;
        writeln!(w)?;

        for node in 1..=self.header.leng {
            write!(w, "{node:>7}")?;
            row(&mut w, m.match_row(node).unwrap_or_default())?;
            for a in &m.annotations[node - 1] {
                write!(w, " {a}")?;
            }
            writeln!(w)?;

            write!(w, "       ")?;
            row(&mut w, m.insert_row(node).unwrap_or_default())?;
            writeln!(w)?;

            write!(w, "       ")?;
            row(&mut w, m.transition_row(node).unwrap_or_default())?;
            writeln!(w)?;
        }

        writeln!(w, "//")?;
        Ok(())
    }
}

impl Display for HmmRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut buf = Vec::new();
        self.write_to(&mut buf).map_err(|_| std::fmt::Error)?;
        f.write_str(std::str::from_utf8(&buf).map_err(|_| std::fmt::Error)?)
    }
}

impl HmmHeader {
    pub fn write_to<W: Write>(&self, mut w: W) -> Result<()> {
        writeln!(w, "{}", self.version)?;
        writeln!(w, "NAME  {}", self.name)?;
        if let Some(v) = &self.acc {
            writeln!(w, "ACC   {v}")?;
        }
        if let Some(v) = &self.desc {
            writeln!(w, "DESC  {v}")?;
        }
        writeln!(w, "LENG  {}", self.leng)?;
        if let Some(v) = self.maxl {
            writeln!(w, "MAXL  {v}")?;
        }
        writeln!(w, "ALPH  {}", self.alph.as_str())?;
        for (key, flag) in [
            ("RF", self.rf),
            ("MM", self.mm),
            ("CONS", self.cons),
            ("CS", self.cs),
            ("MAP", self.map),
        ] {
            writeln!(
                w,
                "{key}{}{}",
                " ".repeat(6 - key.len()),
                if flag { "yes" } else { "no" }
            )?;
        }
        if let Some(v) = &self.date {
            writeln!(w, "DATE  {v}")?;
        }
        for v in &self.com {
            writeln!(w, "COM   {v}")?;
        }
        if let Some(v) = self.nseq {
            writeln!(w, "NSEQ  {v}")?;
        }
        if let Some(v) = self.effn {
            writeln!(w, "EFFN  {v:.6}")?;
        }
        if let Some(v) = self.cksum {
            writeln!(w, "CKSUM {v}")?;
        }
        for (key, cutoff) in [("GA", self.ga), ("TC", self.tc), ("NC", self.nc)] {
            if let Some((a, b)) = cutoff {
                writeln!(w, "{key}    {a:.2} {b:.2}")?;
            }
        }
        for (which, stat) in [
            ("MSV", self.msv),
            ("VITERBI", self.viterbi),
            ("FORWARD", self.forward),
        ] {
            if let Some(g) = stat {
                writeln!(
                    w,
                    "STATS LOCAL {which}{} {:.4} {:.5}",
                    " ".repeat(8usize.saturating_sub(which.len())),
                    g.tau,
                    g.lambda
                )?;
            }
        }
        for (key, value) in &self.extra {
            writeln!(w, "{key}  {value}")?;
        }
        Ok(())
    }
}

// ---

/// p7hmm's record delimiter: `//` at the start of a line, closing the model it
/// ends.
pub struct HmmDelimiter {}

impl Delimit for HmmDelimiter {
    const DELIM: &'static [u8] = b"//";
    const DELIM_TYPE: Delimiter = Delimiter::Terminating;
}

impl Record for HmmRecord {}

/// A parser for the bytes of one profile HMM.
pub struct HmmParser {}

impl Parse for HmmParser {
    type Record = HmmRecord;

    /// The inverse of [`HmmRecord::write_to`]: the magic, `KEY value` header
    /// lines, then three lines per node in the order the file writes them.
    fn parse(bytes: &[u8]) -> Result<HmmRecord> {
        // **note: unlike FASTA, this and every bail below it
        //         are reachable from a range build_index
        //         framed, because the body is positional.
        //         they are the None IndexedHmm::get can
        //         return and IndexedFasta::get cannot
        let Ok(text) = std::str::from_utf8(bytes) else {
            bail_parse!("a profile is utf-8 text")
        };

        // no blank line is meaningful anywhere in a model, and
        // dropping them here is what lets the body be read
        // positionally
        let mut lines = text
            .lines()
            .map(str::trim_end)
            .filter(|l| !l.trim().is_empty())
            .peekable();

        let mut header = read_header(&mut lines)?;
        let model = read_model(&mut lines, &mut header)?;

        Ok(HmmRecord { header, model })
    }
}

/// Everything from the magic through the line before `HMM`.
fn read_header<'a, I>(lines: &mut std::iter::Peekable<I>) -> Result<HmmHeader>
where
    I: Iterator<Item = &'a str>,
{
    let Some(version) = lines.next() else {
        bail_parse!("a profile starts with the HMMER3 magic, found nothing")
    };

    // the analogue of FASTA's '>': the one check that says
    // these bytes are a profile at all
    if !version.trim_start().starts_with("HMMER3") {
        bail_parse!(
            "a profile starts with the HMMER3 magic, found {:?}",
            version.chars().take(32).collect::<String>()
        )
    }

    let mut h = HmmHeader {
        version: version.trim().to_string(),
        ..Default::default()
    };

    while let Some(&line) = lines.peek() {
        let Some((key, value)) = split_first(line) else {
            break;
        };
        if key == "HMM" {
            break;
        }
        lines.next();

        match key {
            "NAME" => h.name = value.to_string(),
            "ACC" => h.acc = Some(value.to_string()),
            "DESC" => h.desc = Some(value.to_string()),
            "LENG" => h.leng = number(value, key)?,
            "MAXL" => h.maxl = Some(number(value, key)?),
            "ALPH" => h.alph = Alphabet::from(value),
            "RF" => h.rf = yes_no(value, key)?,
            "MM" => h.mm = yes_no(value, key)?,
            "CONS" => h.cons = yes_no(value, key)?,
            "CS" => h.cs = yes_no(value, key)?,
            "MAP" => h.map = yes_no(value, key)?,
            "DATE" => h.date = Some(value.to_string()),
            "COM" => h.com.push(value.to_string()),
            "NSEQ" => h.nseq = Some(number(value, key)?),
            "EFFN" => h.effn = Some(number(value, key)?),
            "CKSUM" => h.cksum = Some(number(value, key)?),
            "GA" => h.ga = Some(cutoff(value, key)?),
            "TC" => h.tc = Some(cutoff(value, key)?),
            "NC" => h.nc = Some(cutoff(value, key)?),
            "STATS" => {
                let fields: Vec<&str> = value.split_whitespace().collect();
                let [_config, which, tau, lambda] = fields[..] else {
                    bail_parse!(
                        "a STATS line is <config> <distribution> <tau> <lambda>, \
                         found {line:?}"
                    )
                };

                let g = Gumbel {
                    tau: number(tau, "a STATS tau")?,
                    lambda: number(lambda, "a STATS lambda")?,
                };

                match which {
                    "MSV" => h.msv = Some(g),
                    "VITERBI" => h.viterbi = Some(g),
                    "FORWARD" => h.forward = Some(g),

                    // reported rather than kept in extra: STATS
                    // is a modelled line, so an unknown
                    // distribution inside one is a surprise
                    _ => bail_parse!("unknown STATS distribution {which:?}"),
                }
            }
            _ => h.extra.push((key.to_string(), value.to_string())),
        }
    }

    Ok(h)
}

/// The `HMM` line onwards, through the last node. Takes the header because
/// `LENG` says how many nodes to read and the body carries no count of its own.
fn read_model<'a, I>(lines: &mut std::iter::Peekable<I>, header: &mut HmmHeader) -> Result<HmmModel>
where
    I: Iterator<Item = &'a str>,
{
    let Some(hmm_line) = lines.next() else {
        bail_parse!("the profile ends before its HMM line")
    };

    let symbols = hmm_line
        .split_whitespace()
        .skip(1)
        .map(|t| {
            let mut chars = t.chars();
            match (chars.next(), chars.next()) {
                (Some(c), None) => Ok(c),
                _ => bail_parse!("an alphabet symbol is one character, found {t:?}"),
            }
        })
        .collect::<Result<Vec<char>>>()?;

    if symbols.is_empty() {
        bail_parse!("the HMM line names no alphabet symbols: {hmm_line:?}")
    }
    let k = symbols.len();

    let Some(transition_header) = lines.next() else {
        bail_parse!("the profile ends before its transition header")
    };

    // optional, and the only line of the body that may or may
    // not be there -- everything after it is positional
    let compo = match lines.peek().and_then(|l| split_first(l)) {
        Some(("COMPO", rest)) => {
            lines.next();
            Some(scores(rest, k, "COMPO")?)
        }
        _ => None,
    };

    let Some(line) = lines.next() else {
        bail_parse!("the profile ends before the begin node's insert emissions")
    };
    let begin_insert = scores(line, k, "the begin node's insert emissions")?;

    let Some(line) = lines.next() else {
        bail_parse!("the profile ends before the begin node's transitions")
    };
    let begin_transitions = scores(line, N_TRANSITIONS, "the begin node's transitions")?;

    let leng = header.leng;
    let mut match_emissions = Vec::with_capacity(leng * k);
    let mut insert_emissions = Vec::with_capacity(leng * k);
    let mut transitions = Vec::with_capacity(leng * N_TRANSITIONS);
    let mut annotations: Vec<NodeAnnotation> = Vec::with_capacity(leng);

    for node in 1..=leng {
        let Some(line) = lines.next() else {
            bail_parse!("the profile ends at node {node} of {leng}")
        };

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 1 + k {
            bail_parse!(
                "node {node} has {} fields, and an index plus {k} match emissions is {}",
                fields.len(),
                1 + k
            )
        }

        // the file numbers its own nodes, so a mismatch says
        // the body and LENG have come apart -- which is
        // otherwise silent until a row reads as another's
        let index: usize = number(fields[0], "a node index")?;
        if index != node {
            bail_parse!("expected node {node}, found node {index}")
        }

        for field in &fields[1..1 + k] {
            match_emissions.push(score(field, &format!("node {node}'s match emissions"))?);
        }

        // whichever of MAP, CONS, RF, MM and CS the header
        // switched on, kept verbatim rather than decoded --
        // which is why the count is checked and not the
        // meaning
        let annotation: NodeAnnotation = fields[1 + k..].iter().map(|s| s.to_string()).collect();

        // a row one emission short slides an annotation column
        // into the scores, and the MAP index is a number that
        // parses. HMMER writes every node through the same
        // format string, so an odd row out is what catches it
        if let Some(first) = annotations.first()
            && first.len() != annotation.len()
        {
            bail_parse!(
                "node {node} carries {} annotation columns and node 1 carries {}",
                annotation.len(),
                first.len()
            )
        }
        annotations.push(annotation);

        let Some(line) = lines.next() else {
            bail_parse!("node {node} has no insert emissions")
        };
        insert_emissions.extend(scores(line, k, "insert emissions")?);

        let Some(line) = lines.next() else {
            bail_parse!("node {node} has no transitions")
        };
        transitions.extend(scores(line, N_TRANSITIONS, "transitions")?);
    }

    // a LENG below the real node count leaves whole nodes
    // unread, and nothing above says so: the terminator
    // is where that shows up
    match lines.next() {
        None => {}
        Some(line) if line.trim_start().starts_with("//") => {}
        Some(line) => bail_parse!("expected // after node {leng}, found {line:?}"),
    }

    Ok(HmmModel {
        symbols,
        transition_header: transition_header.to_string(),
        compo,
        begin_insert,
        begin_transitions,
        match_emissions,
        insert_emissions,
        transitions,
        annotations,
    })
}

/// The `NAME` header value in the bytes of one framed profile, read without
/// parsing it. `None` where the header carries none.
pub fn name_of(record: &[u8]) -> Option<&[u8]> {
    for line in record.split(|&b| b == b'\n') {
        let line = line.trim_ascii();
        if line.is_empty() {
            continue;
        }

        let end = line
            .iter()
            .position(|b| b.is_ascii_whitespace())
            .unwrap_or(line.len());
        let (key, value) = line.split_at(end);

        // the header ends at HMM and the body below it is
        // positional, so a NAME past here is not one
        if key == b"HMM" {
            return None;
        }

        if key == b"NAME" {
            return Some(value.trim_ascii());
        }
    }

    None
}

/// The first whitespace-separated token and the rest of the line, trimmed.
fn split_first(line: &str) -> Option<(&str, &str)> {
    let line = line.trim_start();
    if line.is_empty() {
        return None;
    }

    let end = line.find(char::is_whitespace).unwrap_or(line.len());

    Some((&line[..end], line[end..].trim()))
}

/// One score, where the file's `*` is [`f32::INFINITY`].
fn score(token: &str, what: &str) -> Result<f32> {
    if token == "*" {
        return Ok(f32::INFINITY);
    }

    match token.parse::<f32>() {
        Ok(v) => Ok(v),
        Err(e) => bail_parse!("{what}: {token:?} is not a score: {e}"),
    }
}

/// Exactly `n` scores, so a row that is short or long is caught where it sits.
fn scores(line: &str, n: usize, what: &str) -> Result<Vec<f32>> {
    let fields: Vec<&str> = line.split_whitespace().collect();

    if fields.len() != n {
        bail_parse!(
            "{what}: expected {n} scores, found {}: {line:?}",
            fields.len()
        )
    }

    fields.iter().map(|f| score(f, what)).collect()
}

fn number<T>(value: &str, what: &str) -> Result<T>
where
    T: std::str::FromStr,
    T::Err: Display,
{
    match value.parse() {
        Ok(v) => Ok(v),
        Err(e) => bail_parse!("{what}: {value:?} is not a number: {e}"),
    }
}

fn yes_no(value: &str, key: &str) -> Result<bool> {
    match value {
        "yes" => Ok(true),
        "no" => Ok(false),
        _ => bail_parse!("{key} is yes or no, found {value:?}"),
    }
}

/// A `GA`/`TC`/`NC` threshold pair.
fn cutoff(value: &str, key: &str) -> Result<(f32, f32)> {
    // Pfam ends the line with a ';' that HMMER's own writer
    // omits, so both spellings have to read
    let value = value.trim_end_matches(';');
    let fields: Vec<&str> = value.split_whitespace().collect();

    let [a, b] = fields[..] else {
        bail_parse!("{key} is two thresholds, found {value:?}")
    };

    Ok((number(a, key)?, number(b, key)?))
}

// ---

/// Every profile of a p7hmm file, held in memory. Peak memory is roughly the
/// bytes plus the records; [`IndexedHmm`] suits a larger file.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Hmm {
    records: Vec<HmmRecord>,
}

impl Hmm {
    /// Read and parse every profile from `bytes`.
    pub fn new<R: Read>(mut bytes: R) -> Result<Self> {
        let mut all = Vec::new();
        bytes.read_to_end(&mut all)?;

        let offsets = build_index::<_, HmmDelimiter>(&all[..])?;

        // a loop rather than collect::<Result<Vec<_>>>, whose
        // size_hint floor is 0: the spine would double from 1
        // and land on the next power of two, which is 302 MB
        // of capacity for 193 MB of records on a 1 GiB file
        let mut records = Vec::with_capacity(offsets.len());
        for (i, o) in offsets.iter().enumerate() {
            let range = &all[o.start as usize..][..o.n_bytes];

            // which profile failed matters as much as why:
            // the positional body makes this reachable, so
            // the message has to place the record
            records.push(
                HmmParser::parse(range)
                    .map_err(|e| Error::parse(format!("record {i} at byte {}: {e}", o.start)))?,
            );
        }

        Ok(Hmm { records })
    }

    /// Read and parse every profile in the file at `path`.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Hmm::new(File::open(path)?)
    }
}

impl Indexable for Hmm {
    type Record = HmmRecord;

    // the profiles are right here, so lend one
    type Yield<'a> = &'a HmmRecord;

    fn len(&self) -> usize {
        self.records.len()
    }

    fn get(&self, n: usize) -> Option<&HmmRecord> {
        self.records.get(n)
    }
}

// ---

/// Every profile of a p7hmm file, addressed on disk instead of held in memory.
/// Keeps 16 bytes per record and reads one on demand, so it suits sparse
/// access. Unix and Windows only; [`MmapHmm`] maps the file instead.
///
/// # Panics
///
/// Iterating panics on a record that does not parse, since [`Indexable`] reads
/// a `None` from `get` as a broken contract. [`validate`](Self::validate)
/// clears a file of unknown provenance up front.
#[cfg(any(unix, windows))]
pub struct IndexedHmm<S: Source = FileSource> {
    source: S,
    index: Index,
}

#[cfg(any(unix, windows))]
impl IndexedHmm<FileSource> {
    /// Index the file at `path`, parsing none of it.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let source = FileSource::open(path)?;
        let offsets = build_index::<_, HmmDelimiter>(source.file())?;

        Ok(IndexedHmm {
            source,
            index: Index::new(Format::Hmm, offsets),
        })
    }

    /// Address the file at `path` with `index` instead of scanning for it.
    ///
    /// # Errors
    ///
    /// An index of another format. Nothing here checks that it describes this
    /// file; [`Index::load_from`] does that when one comes off disk.
    pub fn with_index(path: impl AsRef<Path>, index: Index) -> Result<Self> {
        if index.format() != Format::Hmm {
            bail_parse!(
                "an index of {} cannot address {}, which is hmm",
                index.format(),
                path.as_ref().display()
            )
        }

        let source = FileSource::open(path)?;

        Ok(IndexedHmm { source, index })
    }
}

/// The same collection with the file mapped rather than read. See
/// [`MmapSource`] for what that costs.
#[cfg(feature = "mmap")]
pub type MmapHmm = IndexedHmm<MmapSource>;

#[cfg(feature = "mmap")]
impl IndexedHmm<MmapSource> {
    /// Map the file at `path` and index it, parsing none of it.
    ///
    /// # Safety
    ///
    /// The file must not be modified or truncated while the collection is
    /// alive; see [`MmapSource::open`].
    pub unsafe fn map(path: impl AsRef<Path>) -> Result<Self> {
        // SAFETY: the caller has promised, through this
        //         function's own contract, not to modify
        //         or truncate the file while the
        //         collection is alive
        let source = unsafe { MmapSource::open(path)? };
        let offsets = build_index::<_, HmmDelimiter>(source.bytes())?;

        Ok(IndexedHmm {
            source,
            index: Index::new(Format::Hmm, offsets),
        })
    }
}

#[cfg(any(unix, windows))]
impl<S: Source> IndexedHmm<S> {
    /// Where each profile sits, and what each is called once
    /// [`name_index`](Self::name_index) has run.
    pub fn index(&self) -> &Index {
        &self.index
    }

    /// Name every profile, by a byte scan of the range the index frames rather
    /// than a parse of it.
    pub fn name_index(&mut self) -> Result<()> {
        // the fields are borrowed apart because name_with
        // takes &mut self.index while the scan reads
        // self.source
        let source = &self.source;

        self.index.name_with(|o| {
            let bytes = source.range(o.start, o.n_bytes)?;

            Ok(name_of(&bytes).map(|name| Cow::Owned(name.to_vec())))
        })
    }

    /// Read every profile, returning the first that fails and its index.
    pub fn validate(&self) -> Result<()> {
        // a full pass over the file, which is why open does
        // not do it -- that would be Hmm with extra steps
        for n in 0..self.len() {
            self.try_get(n)
                .map_err(|e| Error::parse(format!("record {n}: {e}")))?;
        }

        Ok(())
    }

    /// The `n`th profile, or `None` past the end, keeping the failure reason
    /// that [`Indexable::get`] has to throw away.
    pub fn try_get(&self, n: usize) -> Result<Option<HmmRecord>> {
        self.record(n)?
            .map(|bytes| HmmParser::parse(&bytes))
            .transpose()
    }

    /// The `n`th profile's framed bytes, or `None` past the end.
    ///
    /// Borrowed from a mapped source and owned from a read one, exactly as
    /// [`Source::range`] hands them over.
    pub fn record(&self, n: usize) -> Result<Option<Cow<'_, [u8]>>> {
        let Some(o) = self.index.offset(n) else {
            return Ok(None);
        };

        Ok(Some(self.source.range(o.start, o.n_bytes)?))
    }
}

#[cfg(any(unix, windows))]
impl<S: Source> Indexable for IndexedHmm<S> {
    type Record = HmmRecord;

    // nothing to lend: the profile is parsed inside the call
    // and did not exist before it
    type Yield<'a>
        = HmmRecord
    where
        Self: 'a;

    fn len(&self) -> usize {
        self.index.len()
    }

    /// Reads and parses the record's bytes on every call. A failure comes back
    /// as `None`; [`try_get`](Self::try_get) keeps the reason.
    fn get(&self, n: usize) -> Option<HmmRecord> {
        self.try_get(n).ok().flatten()
    }

    // get parses a fresh profile and hands it over owned, so
    // the default -- get, borrow, then clone -- would deep
    // copy a record dropped on the next line
    fn cloned(&self, n: usize) -> Option<HmmRecord> {
        self.get(n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_name_scan_agrees_with_the_parsed_record_on_every_fixture_record() {
        // two independent readings of one file: the scan
        // reads bytes and the parser reads a profile, and
        // over real tool output agreeing is what catches
        // either of them drifting
        let path = fixture!("models.hmm");
        let bytes = std::fs::read(path).expect(path);
        let offsets = build_index::<_, HmmDelimiter>(&bytes[..]).unwrap();
        assert_eq!(offsets.len(), 2);

        for o in &offsets {
            let range = &bytes[o.start as usize..][..o.n_bytes];
            let record = HmmParser::parse(range).unwrap();

            assert_eq!(
                name_of(range),
                Some(record.header.name.as_bytes()),
                "the profile at byte {}",
                o.start
            );
        }
    }

    use crate::fixture;

    /// A two-node model over a three-symbol alphabet.
    fn model() -> HmmModel {
        HmmModel {
            symbols: vec!['A', 'C', 'G'],
            transition_header: "  m->m   m->i   m->d".into(),
            compo: Some(vec![1.0, 2.0, 3.0]),
            begin_insert: vec![0.1, 0.2, 0.3],
            begin_transitions: vec![0.0; N_TRANSITIONS],
            match_emissions: vec![11.0, 12.0, 13.0, 21.0, 22.0, 23.0],
            insert_emissions: vec![1.1, 1.2, 1.3, 2.1, 2.2, 2.3],
            transitions: vec![0.0; 2 * N_TRANSITIONS],
            annotations: vec![vec!["-".into()], vec!["-".into()]],
        }
    }

    #[test]
    fn nodes_are_numbered_from_one_as_the_file_numbers_them() {
        let m = model();
        assert_eq!(m.k(), 3);
        assert_eq!(m.match_row(1), Some(&[11.0, 12.0, 13.0][..]));
        assert_eq!(m.match_row(2), Some(&[21.0, 22.0, 23.0][..]));
        assert_eq!(m.insert_row(1), Some(&[1.1, 1.2, 1.3][..]));
        assert_eq!(m.transition_row(1).map(<[f32]>::len), Some(N_TRANSITIONS));
    }

    #[test]
    fn node_zero_and_past_the_end_are_none_not_a_panic() {
        // node 0 is the begin node and has its own
        // fields; asking the body for it must not wrap
        // around to the last row
        let m = model();
        assert_eq!(m.match_row(0), None);
        assert_eq!(m.insert_row(0), None);
        assert_eq!(m.transition_row(0), None);
        assert_eq!(m.match_row(3), None);
        assert_eq!(m.match_row(usize::MAX), None);
    }

    #[test]
    fn alphabets_are_recognised_by_name_and_anything_else_is_other() {
        assert!(matches!(Alphabet::from("amino"), Alphabet::Amino));
        assert!(matches!(Alphabet::from("DNA"), Alphabet::Dna));
        assert!(matches!(Alphabet::from("RNA"), Alphabet::Rna));
        assert_eq!(Alphabet::Amino.as_str(), "amino");

        let other = Alphabet::from("custom");
        assert!(matches!(other, Alphabet::Other(_)), "{other:?}");
        assert_eq!(
            other.as_str(),
            "custom",
            "an unknown name survives verbatim"
        );
    }

    #[test]
    fn a_gumbel_turns_a_bit_score_into_a_p_value() {
        // higher scores are rarer, which is the whole point
        let g = Gumbel {
            tau: 0.0,
            lambda: 0.7,
        };
        let low = g.p_value(5.0);
        let high = g.p_value(20.0);
        assert!(
            high < low,
            "a better score must be less probable: {high} vs {low}"
        );
        assert!((0.0..=1.0).contains(&low) && (0.0..=1.0).contains(&high));
    }

    #[test]
    fn a_header_without_stats_reports_no_p_value_for_its_threshold() {
        let mut h = HmmHeader::default();
        assert_eq!(h.ga_p_value(), None, "no GA and no fitted distribution");

        h.ga = Some((20.0, 20.0));
        assert_eq!(h.ga_p_value(), None, "GA alone is not enough");

        h.forward = Some(Gumbel {
            tau: 0.0,
            lambda: 0.7,
        });
        let p = h.ga_p_value().expect("GA and a fitted distribution");
        assert!((0.0..=1.0).contains(&p), "{p}");
    }

    #[test]
    fn a_record_writes_a_header_a_body_and_a_terminator() {
        let rec = HmmRecord {
            header: HmmHeader {
                version: "HMMER3/f".into(),
                name: "demo".into(),
                leng: 2,
                alph: Alphabet::Amino,
                ..Default::default()
            },
            model: model(),
        };

        let mut out = Vec::new();
        rec.write_to(&mut out).unwrap();
        let text = String::from_utf8(out).unwrap();

        assert!(text.starts_with("HMMER3/f"), "{text}");
        assert!(text.trim_end().ends_with("//"), "{text}");
        for needle in ["NAME  demo", "LENG  2", "ALPH  amino", "HMM  "] {
            assert!(text.contains(needle), "missing {needle:?} in:\n{text}");
        }
        assert_eq!(rec.to_string(), text, "Display and write_to must agree");
    }

    // ---

    use std::borrow::Borrow;

    use crate::collection::Iterable;

    /// The record `model()` belongs to, written out and read back by the
    /// round-trip tests.
    fn record() -> HmmRecord {
        HmmRecord {
            header: HmmHeader {
                version: "HMMER3/f [3.4 | Aug 2023]".into(),
                name: "demo".into(),
                acc: Some("DM00001".into()),
                desc: Some("a synthetic profile".into()),
                leng: 2,
                alph: Alphabet::Amino,
                cons: true,
                map: true,
                nseq: Some(24),
                effn: Some(3.5),
                cksum: Some(1234567890),
                ga: Some((22.0, 22.0)),
                forward: Some(Gumbel {
                    tau: -3.5,
                    lambda: 0.71,
                }),
                ..Default::default()
            },
            model: model(),
        }
    }

    fn written(rec: &HmmRecord) -> Vec<u8> {
        let mut out = Vec::new();
        rec.write_to(&mut out).unwrap();
        out
    }

    #[test]
    fn parsing_a_written_profile_returns_what_was_written() {
        // pinning parse against write_to checks both
        // directions with no hand-written expected value.
        // every optional header field is set here, since a
        // field silently dropped is what round-tripping catches
        let rec = record();
        assert_eq!(HmmParser::parse(&written(&rec)).unwrap(), rec);
    }

    #[test]
    fn a_profile_with_only_the_required_header_round_trips() {
        // the other end of the range: no ACC, no DESC, no
        // cutoffs, no fitted distributions
        let rec = HmmRecord {
            header: HmmHeader {
                name: "bare".into(),
                leng: 2,
                ..Default::default()
            },
            model: model(),
        };

        assert_eq!(HmmParser::parse(&written(&rec)).unwrap(), rec);
    }

    #[test]
    fn a_star_is_an_infinite_score_in_both_directions() {
        // the file writes '*' for a transition the model
        // forbids, and f32::INFINITY is what the body holds
        let mut rec = record();
        rec.model.transitions[2] = f32::INFINITY;

        let text = String::from_utf8(written(&rec)).unwrap();
        assert!(text.contains(" *"), "{text}");

        let back = HmmParser::parse(text.as_bytes()).unwrap();
        assert!(back.model.transition_row(1).unwrap()[2].is_infinite());
        assert_eq!(back, rec);
    }

    #[test]
    fn header_flags_read_as_the_yes_and_no_the_file_writes() {
        let rec = HmmParser::parse(&written(&record())).unwrap();
        assert!(rec.header.cons, "CONS  yes");
        assert!(rec.header.map, "MAP   yes");
        assert!(!rec.header.rf, "RF    no");

        let broken = String::from_utf8(written(&record()))
            .unwrap()
            .replace("RF    no", "RF    maybe");
        let err = HmmParser::parse(broken.as_bytes()).unwrap_err();
        assert!(err.to_string().contains("yes or no"), "{err}");
    }

    #[test]
    fn a_header_line_this_parser_does_not_model_survives_verbatim() {
        // what extra is for: an unknown key is kept in order
        // so a profile written back still carries it
        let text = String::from_utf8(written(&record()))
            .unwrap()
            .replace("HMM  ", "PARTICULAR  wat\nHMM  ");

        let rec = HmmParser::parse(text.as_bytes()).unwrap();
        assert_eq!(rec.header.extra, [("PARTICULAR".into(), "wat".into())]);
        assert_eq!(
            HmmParser::parse(&written(&rec)).unwrap(),
            rec,
            "and it survives a second trip"
        );
    }

    #[test]
    fn a_pfam_threshold_line_may_end_in_a_semicolon() {
        // HMMER's own writer omits it and Pfam's files carry
        // it, so both spellings have to read
        let text = String::from_utf8(written(&record()))
            .unwrap()
            .replace("GA    22.00 22.00", "GA    22.00 22.00;");

        assert_eq!(
            HmmParser::parse(text.as_bytes()).unwrap().header.ga,
            Some((22.0, 22.0))
        );
    }

    #[test]
    fn bytes_that_do_not_start_a_profile_say_so() {
        for bad in [&b">a\nACGT\n"[..], &b""[..], &b"# STOCKHOLM 1.0\n//\n"[..]] {
            let err = HmmParser::parse(bad).unwrap_err();
            assert!(matches!(err.kind(), crate::ErrorKind::Parse { .. }));
            assert!(
                err.to_string().contains("HMMER3"),
                "the message should name what was expected: {err}"
            );
        }
    }

    #[test]
    fn a_node_numbered_out_of_step_with_the_body_is_rejected() {
        // the file numbers its own nodes, and a mismatch is
        // the only sign that LENG and the body came apart
        // before a row is read as another node's
        let text = String::from_utf8(written(&record()))
            .unwrap()
            .replace("      2 ", "      7 ");

        let err = HmmParser::parse(text.as_bytes()).unwrap_err();
        assert!(err.to_string().contains("expected node 2"), "{err}");
    }

    #[test]
    fn a_length_smaller_than_the_body_is_rejected() {
        // otherwise the last node is read as nothing at all,
        // and the profile comes back quietly short
        let text = String::from_utf8(written(&record()))
            .unwrap()
            .replace("LENG  2", "LENG  1");

        let err = HmmParser::parse(text.as_bytes()).unwrap_err();
        assert!(err.to_string().contains("expected //"), "{err}");
    }

    #[test]
    fn a_truncated_body_names_the_node_it_stopped_at() {
        let text = String::from_utf8(written(&record())).unwrap();
        let cut = text.find("      2 ").unwrap();

        let err = HmmParser::parse(&text.as_bytes()[..cut]).unwrap_err();
        assert!(err.to_string().contains("node 2 of 2"), "{err}");
    }

    #[test]
    fn a_match_row_one_emission_short_is_caught_by_the_node_it_ruins() {
        // the annotation columns are kept verbatim, so a short
        // row does not run off the end -- it slides an
        // annotation into the scores instead. here that is
        // "-", which is not a number
        let text = String::from_utf8(written(&record()))
            .unwrap()
            .replace(" 11.00000 12.00000 13.00000", " 11.00000 12.00000");

        let err = HmmParser::parse(text.as_bytes()).unwrap_err();
        assert!(err.to_string().contains("node 1"), "{err}");
    }

    #[test]
    fn a_short_row_is_caught_even_when_the_column_that_slides_in_is_a_number() {
        // the case the score check cannot see: HMMER writes
        // the MAP index first, and a node one emission short
        // reads it as the last score without complaint. every
        // node is written through one format string, so the
        // node left with an annotation column too few is the
        // signal
        let mut rec = record();
        rec.model.annotations = vec![vec!["1".into(), "p".into()], vec!["2".into(), "g".into()]];

        let text = String::from_utf8(written(&rec))
            .unwrap()
            .replace(" 21.00000 22.00000 23.00000", " 21.00000 22.00000");

        let err = HmmParser::parse(text.as_bytes()).unwrap_err();
        assert!(
            err.to_string().contains("annotation columns"),
            "a numeric annotation slid into the scores unnoticed: {err}"
        );
    }

    // ---

    fn names<C: Indexable<Record = HmmRecord>>(c: &C) -> Vec<String> {
        c.iter().map(|r| r.borrow().header.name.clone()).collect()
    }

    fn two_profiles() -> Vec<u8> {
        let mut out = written(&record());
        let mut second = record();
        second.header.name = "other".into();
        second.header.leng = 2;
        out.extend_from_slice(&written(&second));
        out
    }

    #[test]
    fn a_file_becomes_its_profiles_in_file_order() {
        let hmm = Hmm::new(&two_profiles()[..]).unwrap();

        assert_eq!(hmm.len(), 2);
        assert!(!hmm.is_empty());
        assert_eq!(names(&hmm), ["demo", "other"]);
        assert!(hmm.get(2).is_none(), "past the end is None, not a panic");
    }

    #[test]
    fn an_input_with_no_profiles_is_an_empty_collection() {
        let hmm = Hmm::new(&b""[..]).unwrap();
        assert_eq!(hmm.len(), 0);
        assert!(hmm.is_empty());
        assert_eq!(hmm.iter().count(), 0);
    }

    #[test]
    fn a_malformed_profile_names_the_record_that_failed() {
        // the in-memory reader parses everything up front, so
        // it can refuse the file rather than hand back a
        // collection with a hole in it
        let mut bytes = written(&record());
        bytes.extend_from_slice(b"not a profile\n//\n");

        let err = Hmm::new(&bytes[..]).unwrap_err();
        assert!(err.to_string().contains("record 1"), "{err}");
        assert!(err.to_string().contains("HMMER3"), "{err}");
    }

    #[test]
    fn the_combinators_work_on_profiles_too() {
        let bytes = two_profiles();
        let named = Hmm::new(&bytes[..])
            .unwrap()
            .filter(|r| r.header.name == "other");
        assert_eq!(names(&named), ["other"]);

        let by_name = Hmm::new(&bytes[..])
            .unwrap()
            .sort_by_key(|r| std::cmp::Reverse(r.header.name.clone()));
        assert_eq!(names(&by_name), ["other", "demo"]);
        assert_eq!(Hmm::new(&bytes[..]).unwrap().sample(2).len(), 2);
    }

    #[test]
    fn every_profile_in_the_real_fixture_parses() {
        let path = fixture!("models.hmm");
        let hmm = Hmm::open(path).expect(path);

        // the two profiles hmmbuild makes from families.sto
        assert_eq!(hmm.len(), 2);
        assert_eq!(names(&hmm), ["SH3_1", "PDZ"]);

        let sh3 = hmm.get(0).unwrap();
        assert_eq!(sh3.header.leng, 46);
        assert_eq!(sh3.model.k(), 20);
        assert_eq!(sh3.header.acc.as_deref(), Some("PF00018.35"));

        // GA/TC/NC are Pfam's own curated cutoffs, carried through
        // hmmbuild from the SEED's #=GF GA rather than invented here
        assert_eq!(sh3.header.ga, Some((22.90, 22.90)));
        assert_eq!(sh3.header.nc, Some((22.80, 22.80)));
        assert!(sh3.model.compo.is_some(), "the fixture has a COMPO line");
        assert_eq!(sh3.model.match_row(46).map(<[f32]>::len), Some(20));
        assert_eq!(sh3.model.match_row(47), None, "LENG is 46");
        assert!(
            sh3.header.ga_p_value().unwrap() < 1e-3,
            "GA is a threshold a real hit should clear"
        );

        // MAP and CONS are both on, so every match line
        // carries trailing annotation the parser kept
        assert!(!sh3.model.annotations[0].is_empty());

        for n in 0..hmm.len() {
            let rec = hmm.get(n).unwrap();
            assert_eq!(
                HmmParser::parse(&written(rec)).unwrap(),
                *rec,
                "record {n} does not survive a round trip"
            );
        }
    }

    // ---

    /// Write `bytes` somewhere [`IndexedHmm`] can open them, since it needs a
    /// real file descriptor to read at an offset.
    fn temp_hmm(tag: &str, bytes: &[u8]) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("libsail-{tag}-{}.hmm", std::process::id()));
        std::fs::write(&p, bytes).unwrap();
        p
    }

    #[test]
    fn the_disk_backed_reader_agrees_with_the_in_memory_one() {
        // the claim of Indexable is that one API covers both
        // backends; this is what checks it beyond the types
        let bytes = two_profiles();
        let path = temp_hmm("agree", &bytes);
        let mem = Hmm::new(&bytes[..]).unwrap();
        let disk = IndexedHmm::open(&path).unwrap();

        assert_eq!(disk.len(), mem.len());
        for n in 0..mem.len() {
            assert_eq!(disk.cloned(n), mem.cloned(n), "record {n}");
        }

        assert!(disk.get(mem.len()).is_none(), "past the end is None");
        assert!(disk.try_get(mem.len()).unwrap().is_none());
        disk.validate().unwrap();
        assert_eq!(names(&disk), ["demo", "other"]);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn validate_finds_the_record_get_can_only_call_none() {
        // the gap FASTA does not have: parsing an indexed
        // range can fail here, so get() has a None that
        // Indexable treats as a broken contract
        let mut bytes = written(&record());
        bytes.extend_from_slice(b"not a profile\n//\n");

        let path = temp_hmm("bad", &bytes);
        let disk = IndexedHmm::open(&path).unwrap();

        assert_eq!(disk.len(), 2);
        assert!(disk.get(0).is_some());
        assert!(disk.get(1).is_none(), "get has nowhere to put the reason");

        let err = disk.validate().unwrap_err();
        assert!(err.to_string().contains("record 1"), "{err}");

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn a_file_that_shrinks_under_the_index_is_an_io_failure_not_a_parse_one() {
        // the other way a record fails here: the offsets
        // outlive the bytes they were taken from, and the
        // positioned read is what reports it
        //
        // **note: FileSource only, always. under MmapSource a
        //         truncation raises SIGBUS, so this test
        //         would abort the process rather than fail
        let bytes = two_profiles();
        let path = temp_hmm("truncated", &bytes);
        let disk = IndexedHmm::open(&path).unwrap();
        assert_eq!(disk.len(), 2);

        std::fs::write(&path, written(&record())).unwrap();

        assert!(disk.get(1).is_none(), "get has nowhere to put the reason");
        let err = disk.try_get(1).unwrap_err();
        assert!(matches!(err.kind(), crate::ErrorKind::Io { .. }), "{err}");
        assert!(err.to_string().contains("fill whole buffer"), "{err}");

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn the_two_backends_agree_on_the_real_fixture() {
        let path = fixture!("models.hmm");
        let (mem, disk) = (
            Hmm::open(path).expect(path),
            IndexedHmm::open(path).expect(path),
        );

        assert_eq!(disk.len(), mem.len());
        for n in 0..mem.len() {
            assert_eq!(disk.cloned(n), mem.cloned(n), "record {n}");
        }
    }

    #[test]
    fn an_index_handed_in_addresses_the_file_the_same_as_one_scanned_for() {
        // what writing an index down is worth: with_index
        // skips the scan, and nothing downstream may be able
        // to tell that it did
        let path = fixture!("models.hmm");
        let scanned = IndexedHmm::open(path).expect(path);

        let offsets = build_index::<_, HmmDelimiter>(File::open(path).unwrap()).unwrap();
        let handed = IndexedHmm::with_index(path, Index::new(Format::Hmm, offsets)).expect(path);

        assert_eq!(handed.len(), scanned.len());
        for n in 0..scanned.len() {
            assert_eq!(handed.cloned(n), scanned.cloned(n), "record {n}");
        }
    }

    #[test]
    fn an_index_of_another_format_is_refused() {
        // offsets alone cannot say what they frame, so the
        // format the index carries is the only thing between
        // a fasta index and a file it does not describe
        let path = fixture!("models.hmm");
        let offsets = build_index::<_, HmmDelimiter>(File::open(path).unwrap()).unwrap();

        assert!(IndexedHmm::with_index(path, Index::new(Format::Fasta, offsets)).is_err());
    }

    #[test]
    fn naming_the_index_finds_every_record_by_the_name_its_parser_reports() {
        let path = fixture!("models.hmm");
        let mut collection = IndexedHmm::open(path).expect(path);
        collection.name_index().unwrap();

        for n in 0..collection.len() {
            let Some(name) = collection.index().name(n) else {
                continue;
            };

            assert_eq!(collection.index().position(name), Some(n), "record {n}");
        }
    }

    #[cfg(feature = "mmap")]
    #[test]
    fn the_mapped_backend_agrees_with_the_read_one() {
        // a third backend under one API. Source is only
        // worth having if these agree profile for profile
        let path = temp_hmm("mapped", &two_profiles());
        let mem = Hmm::new(&two_profiles()[..]).unwrap();
        let read = IndexedHmm::open(&path).unwrap();

        // SAFETY: this test owns the file and does not touch
        //         it while the map is alive
        let mapped = unsafe { MmapHmm::map(&path) }.unwrap();

        assert_eq!(mapped.len(), mem.len());
        for n in 0..mem.len() {
            assert_eq!(
                mapped.cloned(n),
                mem.cloned(n),
                "profile {n} against memory"
            );
            assert_eq!(
                mapped.cloned(n),
                read.cloned(n),
                "profile {n} against pread"
            );
        }

        assert!(mapped.get(mem.len()).is_none(), "past the end is None");
        assert!(mapped.try_get(mem.len()).unwrap().is_none());
        mapped.validate().unwrap();
        assert_eq!(names(&mapped), ["demo", "other"]);

        // SAFETY: as above
        let sorted = unsafe { MmapHmm::map(&path) }.unwrap().sample(2);
        assert_eq!(sorted.len(), 2, "the combinators reach it too");

        std::fs::remove_file(path).ok();
    }
}

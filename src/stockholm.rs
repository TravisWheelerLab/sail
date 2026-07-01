use std::{
    fmt::Display,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::Path,
};

use anyhow::{anyhow, bail, Context};
use indexmap::IndexMap;

pub const HEADER: &str = "# STOCKHOLM 1.0";

#[derive(Default)]
pub struct StockholmRecord {
    pub id: String,
    pub gf_meta: IndexMap<String, String>,
    pub gs_meta: IndexMap<String, String>,
    pub sequences: IndexMap<String, String>,
}

impl StockholmRecord {
    pub fn parse(lines: Vec<String>) -> anyhow::Result<Self> {
        let mut rec = Self::default();
        let mut lines = lines.iter().filter(|line| !line.is_empty());

        match lines.next() {
            Some(line) => {
                if !line.starts_with("# STOCKHOLM") {
                    bail!("buffer doesn't begin with Stockholm header")
                }
            }
            None => bail!("empty buffer"),
        }

        for line in lines {
            if line.starts_with("#=GF") {
                let mut tokens = line.split_whitespace();
                let key = tokens.nth(1).ok_or(anyhow!("no GF key"))?.to_string();

                if key == "ID" {
                    rec.id = tokens.next().ok_or(anyhow!("no ID value"))?.to_string();
                }

                rec.gf_meta.insert(key, line.clone());
            } else if line.starts_with("#=GS") {
                let mut tokens = line.split_whitespace();
                let name = tokens.nth(1).ok_or(anyhow!("no GS name"))?.to_string();

                rec.gs_meta.insert(name.clone(), line.clone());
                rec.sequences.insert(name, "".to_string());
            } else if line.starts_with("#=GC") || line.starts_with("#=GR") {
                // ignore for now
            } else if !line.starts_with("#") {
                let tokens = line.split_whitespace().collect::<Vec<_>>();
                let name = tokens[0];
                let seq = tokens[1];

                rec.sequences
                    .entry(name.to_string())
                    .and_modify(|s| s.push_str(seq));
            } else {
                bail!("unexpected line format:\n{line}");
            }
        }

        if rec.id.is_empty() {
            bail!("failed to find Stockholm record ID");
        }

        Ok(rec)
    }

    pub fn get(&self, name: &str) -> Option<&String> {
        self.sequences.get(name)
    }
}

impl Display for StockholmRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{HEADER}")?;
        self.gf_meta.values().try_for_each(|s| writeln!(f, "{s}"))?;
        writeln!(f)?;
        self.gs_meta.values().try_for_each(|s| writeln!(f, "{s}"))?;
        writeln!(f)?;
        self.sequences
            .iter()
            .try_for_each(|(n, s)| writeln!(f, "{n} {s}"))?;
        writeln!(f, "//")
    }
}

#[derive(Default)]
pub struct Stockholm {
    pub records: IndexMap<String, StockholmRecord>,
}

impl Stockholm {
    pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let f = std::fs::File::open(path.as_ref())
            .with_context(|| format!("failed to open: {}", path.as_ref().to_string_lossy()))?;
        Self::parse(f)
    }

    pub fn parse<R: Read>(buf: R) -> anyhow::Result<Self> {
        let reader = BufReader::new(buf);

        let mut bufs: Vec<Vec<String>> = vec![];
        let mut buf: Vec<String> = vec![];

        for line in reader.lines() {
            // skip malformed lines for now
            let line = line.unwrap_or_default();

            if line.starts_with("//") {
                bufs.push(buf);
                buf = vec![];
            } else {
                buf.push(line);
            }
        }

        let mut sto = Self::default();

        for buf in bufs {
            let rec = StockholmRecord::parse(buf)?;
            sto.records.insert(rec.id.clone(), rec);
        }

        Ok(sto)
    }

    pub fn get(&self, name: &str) -> Option<&StockholmRecord> {
        self.records.get(name)
    }

    pub fn get_mut(&mut self, name: &str) -> Option<&mut StockholmRecord> {
        self.records.get_mut(name)
    }

    pub fn write<W: Write>(&self, buf: W) -> anyhow::Result<()> {
        let mut out = BufWriter::new(buf);
        self.records
            .values()
            .try_for_each(|r| writeln!(out, "{r}"))
            .context("Stockholm::write() failed")
    }

    pub fn seq_cnt(&self) -> usize {
        self.records.values().map(|r| r.sequences.len()).sum()
    }
}

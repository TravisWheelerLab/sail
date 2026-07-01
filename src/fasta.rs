use std::{
    fmt::Display,
    io::{BufRead, BufReader, Read, Seek},
    path::Path,
};

use anyhow::{anyhow, Context};
use indexmap::IndexMap;

#[derive(Default, Clone, PartialEq, Eq, Hash)]
pub struct FastaRecord {
    pub name: String,
    pub extra: String,
    pub seq: String,
}

impl FastaRecord {
    pub fn reverse(&mut self) {
        self.seq = self.seq.chars().rev().collect::<String>();
    }
}

impl Display for FastaRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, ">{} {}", self.name, self.extra)?;
        let mut chunks = self
            .seq
            .as_bytes()
            .chunks(60)
            .map(|c| std::str::from_utf8(c).unwrap());

        if let Some(last) = chunks.next_back() {
            chunks.try_for_each(|c| writeln!(f, "{c}"))?;
            write!(f, "{last}")
        } else {
            Ok(())
        }
    }
}

#[derive(Default)]
pub struct Fasta {
    pub records: IndexMap<String, FastaRecord>,
}

impl Fasta {
    pub fn parse<R: Read>(buf: R) -> anyhow::Result<Self> {
        let reader = BufReader::new(buf);

        let mut fasta = Self::default();
        let mut rec = FastaRecord::default();
        for line in reader.lines() {
            let line = line?;
            if let Some(line) = line.strip_prefix(">") {
                if !rec.name.is_empty() {
                    fasta.records.insert(rec.name.to_string(), rec);
                }
                rec = FastaRecord::default();

                let mut tokens = line.splitn(2, char::is_whitespace);
                rec.name = tokens.next().ok_or(anyhow!("no name"))?.to_string();
                rec.extra = tokens.next().unwrap_or_default().to_string();
            } else {
                rec.seq.push_str(&line)
            }
        }

        if rec != FastaRecord::default() {
            fasta.records.insert(rec.name.to_string(), rec);
        }

        Ok(fasta)
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let f = std::fs::File::open(path.as_ref())?;
        Self::parse(f).with_context(|| format!("failed to parse Fasta from: {:?}", path.as_ref()))
    }
}

pub struct ByteBlock<const L: usize> {
    pub byte_start: u64,
    pub byte_cnts: [u16; L],
}

impl<const L: usize> ByteBlock<L> {
    fn new(byte_start: u64) -> Self {
        Self {
            byte_start,
            byte_cnts: [0; L],
        }
    }
}

pub struct FastaByteIndex<R, const L: usize>
where
    R: Read + Seek,
{
    pub reader: BufReader<R>,
    pub blocks: Vec<ByteBlock<L>>,
    pub size: usize,
}

impl<R, const L: usize> FastaByteIndex<R, L>
where
    R: Read + Seek,
{
    pub fn new(buf: R) -> anyhow::Result<Self> {
        let mut reader = BufReader::new(buf);
        reader.seek(std::io::SeekFrom::Start(0))?;

        let mut total_record_cnt = 0usize;
        let mut block_record_cnt = 0usize;

        let mut record_byte_cnt = 0usize;
        let mut block_byte_start = 0usize;
        let mut byte_offset = 0usize;

        let mut blocks = vec![];
        let mut block: ByteBlock<L> = ByteBlock::new(0);

        for line in (&mut reader).lines() {
            let line = line?;

            if line.starts_with('>') {
                block_byte_start = byte_offset;
                block.byte_cnts[block_record_cnt] = record_byte_cnt as u16;
                block_record_cnt += 1;
                total_record_cnt += 1;
                record_byte_cnt = 0;
            }

            if block_record_cnt == L {
                blocks.push(block);
                block = ByteBlock::new(block_byte_start as u64);
                block_record_cnt = 0;
            }

            record_byte_cnt += line.len() + 1;
            byte_offset += line.len() + 1;
        }

        block.byte_cnts[block_record_cnt] = record_byte_cnt as u16;
        blocks.push(block);

        Ok(Self {
            reader,
            blocks,
            size: total_record_cnt,
        })
    }

    pub fn get_record(&mut self, index: usize) -> anyhow::Result<FastaRecord> {
        let s = self.get(index)?;
        let lines: Vec<&str> = s.lines().collect();

        let header = lines[0];
        let seq = lines
            .iter()
            .skip(1)
            .copied()
            .collect::<Vec<&str>>()
            .join("");

        let (name, extra) = header
            .split_once(char::is_whitespace)
            .unwrap_or((header, ""));

        Ok(FastaRecord {
            name: name[1..].to_string(),
            extra: extra.to_string(),
            seq,
        })
    }

    pub fn get(&mut self, index: usize) -> anyhow::Result<String> {
        let block = &self.blocks[index / L];
        let block_idx = index % L;

        let before = block
            .byte_cnts
            .iter()
            .take(block_idx)
            .map(|&l| l as u64)
            .sum::<u64>();

        let n_read = block.byte_cnts[block_idx] as usize;

        let mut buf = vec![0u8; n_read];
        self.reader
            .seek(std::io::SeekFrom::Start(block.byte_start + before))?;
        self.reader.read_exact(&mut buf)?;

        String::from_utf8(buf).context("failed to produce string")
    }
}

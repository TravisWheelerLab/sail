use std::{
    fs::File,
    io::{BufRead, BufReader, Read, Seek},
    path::{Path, PathBuf},
};

use anyhow::Context;
use regex::Regex;

#[derive(Clone, Copy)]
struct Descriptor {
    file_idx: usize,
    offset: u64,
    length: u64,
}

pub struct Db {
    file: File,
    paths: Vec<PathBuf>,
    descriptors: Vec<Descriptor>,
}

impl Db {
    pub fn from_path<P>(path: P) -> anyhow::Result<Self>
    where
        P: AsRef<Path>,
    {
        let path = path.as_ref();
        let dir = path.parent().context("path has no parent dir")?;
        let name = path
            .file_stem()
            .context("no file stem")?
            .to_str()
            .context("invalid utf8")?;

        // ---

        let index_path = path.with_extension("index");
        let index_reader =
            BufReader::new(File::open(&index_path).with_context(|| format!("{:?}", index_path))?);

        let mut offsets = vec![];
        let mut lengths = vec![];
        for line in index_reader.lines() {
            let line = line?;

            let tokens = line
                .split('\t')
                .map(|s| s.parse::<u64>().expect("failed to parse index line"))
                .collect::<Vec<_>>();

            offsets.push(tokens[1]);
            // -1 because the index lengths include a null byte
            lengths.push(tokens[2] - 1);
        }

        offsets.shrink_to_fit();
        lengths.shrink_to_fit();

        // ---

        let re = Regex::new(&format!(r"^{}\.\d+$", regex::escape(name)))?;
        let mut paths = std::fs::read_dir(dir)?
            .filter_map(Result::ok)
            .map(|e| e.path())
            .filter(|p| {
                p.file_name()
                    .and_then(|s| s.to_str())
                    .map(|s| re.is_match(s))
                    .unwrap_or(false)
            })
            .collect::<Vec<_>>();

        paths.sort_by_key(|p| {
            p.extension()
                .and_then(|e| e.to_str())
                .and_then(|s| s.parse::<usize>().ok())
                .expect("bad DB split suffix")
        });

        // if we have no paths here, that
        // means the prefilter is one file
        if paths.is_empty() {
            paths.push(path.to_path_buf())
        }

        let file_sizes: Vec<u64> = paths
            .iter()
            .map(|p| p.metadata().map(|m| m.len()))
            .collect::<Result<_, _>>()?;

        let prefix_sum = file_sizes
            .into_iter()
            .scan(0u64, |acc, x| {
                let cur = *acc;
                *acc += x;
                Some(cur)
            })
            .collect::<Vec<_>>();

        // --

        let mut descriptors = offsets
            .into_iter()
            .zip(lengths)
            .map(|(offset, length)| {
                let file_idx = match prefix_sum.binary_search(&offset) {
                    // this means the binary search found the exact offset
                    Ok(idx) => idx,
                    // this means the binary search landed between offsets
                    Err(idx) => idx.saturating_sub(1),
                };

                let relative_offset = offset - prefix_sum[file_idx];
                Descriptor {
                    file_idx,
                    offset: relative_offset,
                    length,
                }
            })
            .collect::<Vec<_>>();

        for desc in descriptors.iter_mut() {
            let mut file = File::open(&paths[desc.file_idx])?;
            file.seek(std::io::SeekFrom::Start(desc.offset + desc.length))?;

            let mut b = [0u8; 1];
            file.read_exact(&mut b)?;
            let byte = b[0];
            assert!(byte == 0);
        }

        let file = File::open(&paths[0])?;

        Ok(Self {
            file,
            paths,
            descriptors,
        })
    }

    pub fn len(&self) -> usize {
        self.descriptors.len()
    }

    pub fn open_file(&mut self, file_idx: usize) -> anyhow::Result<()> {
        self.file = File::open(&self.paths[file_idx])?;
        Ok(())
    }

    pub fn get(&mut self, idx: usize) -> anyhow::Result<String> {
        let desc = self.descriptors[idx];

        if desc.length == 0 {
            return Ok(String::new());
        }

        self.open_file(desc.file_idx)?;

        self.file.seek(std::io::SeekFrom::Start(desc.offset))?;
        let mut taken = (&mut self.file).take(desc.length);

        let mut s = String::new();
        taken.read_to_string(&mut s)?;

        Ok(s)
    }
}

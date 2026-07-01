use std::{
    collections::HashMap,
    fmt::Display,
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
};

use anyhow::Context;

pub mod nail {
    use std::{
        collections::HashMap,
        fs::File,
        io::{BufRead, BufReader, Read},
        path::Path,
    };

    use anyhow::Context;

    #[derive(Clone)]
    pub struct NailHit {
        pub query: String,
        pub target: String,
        pub query_start: usize,
        pub query_end: usize,
        pub target_start: usize,
        pub target_end: usize,
        pub score: f64,
        pub e_value: f64,
        pub cell_frac: f64,
    }

    pub struct NailTable {
        pub name: String,
        pub hits: Vec<NailHit>,
    }

    impl NailTable {
        pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
            let path = path.as_ref();
            let file = File::open(path)?;
            let name = path
                .file_stem()
                .and_then(|s| s.to_str())
                .context("invalid path")?;

            Self::parse::<File>(file, name)
        }

        pub fn parse<R: Read>(buf: R, name: &str) -> anyhow::Result<Self> {
            let reader = BufReader::new(buf);

            let hits = reader
                .lines()
                .filter_map(|line| {
                    let line = line.ok()?;
                    if line.starts_with('#') {
                        None
                    } else {
                        let tokens = line.split_whitespace().collect::<Vec<_>>();
                        Some(NailHit {
                            target: tokens[0].to_string(),
                            query: tokens[1].to_string(),
                            target_start: tokens[2].parse::<usize>().ok()?,
                            target_end: tokens[3].parse::<usize>().ok()?,
                            query_start: tokens[4].parse::<usize>().ok()?,
                            query_end: tokens[5].parse::<usize>().ok()?,
                            score: tokens[6].parse::<f64>().ok()?,
                            e_value: tokens[8].parse::<f64>().ok()?,
                            cell_frac: tokens[9].parse::<f64>().ok()?,
                        })
                    }
                })
                .collect();
            Ok(Self {
                name: name.to_string(),
                hits,
            })
        }

        pub fn to_map(self) -> HashMap<(String, String), NailHit> {
            self.hits
                .into_iter()
                .map(|h| ((h.query.clone(), h.target.clone()), h))
                .collect()
        }
    }
}

pub mod hmmer {
    use std::{
        collections::HashMap,
        fs::File,
        io::{BufRead, BufReader, Read},
        path::Path,
    };

    use anyhow::Context;

    pub struct HmmerHit {
        pub score: f32,
        pub e_value: f64,
        pub domains: Vec<Domain>,
    }

    pub struct Domain {
        pub query_start: usize,
        pub query_end: usize,
        pub target_start: usize,
        pub target_end: usize,
        pub score: f32,
        pub e_value: f64,
    }

    pub struct HmmerDomainTable {
        pub name: String,
        pub hits: HashMap<(String, String), HmmerHit>,
    }

    impl HmmerDomainTable {
        pub fn from_path<P, F>(path: P, filter: F) -> anyhow::Result<Self>
        where
            P: AsRef<Path>,
            F: Fn(&(String, String)) -> bool,
        {
            let path = path.as_ref();

            let file = File::open(path)?;
            let name = path
                .file_stem()
                .and_then(|s| s.to_str())
                .context("invalid path")?;

            Self::parse::<File, F>(file, name, filter)
        }

        pub fn parse<R, F>(buf: R, name: &str, filter: F) -> anyhow::Result<Self>
        where
            R: Read,
            F: Fn(&(String, String)) -> bool,
        {
            let reader = BufReader::new(buf);

            let mut hits: HashMap<(String, String), HmmerHit> = HashMap::new();
            for line in reader.lines() {
                let line = line.unwrap_or_default();

                if line.starts_with('#') {
                    continue;
                }

                let tokens = line.split_whitespace().collect::<Vec<_>>();

                let query = tokens[3].to_string();
                let target = tokens[0].to_string();
                let key = (query, target);

                if !filter(&key) {
                    continue;
                }

                let e_value = tokens[6].parse::<f64>()?;
                let score = tokens[7].parse::<f32>()?;

                let dom_e_value = tokens[12].parse::<f64>()?;
                let dom_score = tokens[13].parse::<f32>()?;

                let query_start = tokens[15].parse::<usize>()?;
                let query_end = tokens[16].parse::<usize>()?;

                let target_start = tokens[17].parse::<usize>()?;
                let target_end = tokens[18].parse::<usize>()?;

                let dom = Domain {
                    query_start,
                    query_end,
                    target_start,
                    target_end,
                    score: dom_score,
                    e_value: dom_e_value,
                };

                match hits.get_mut(&key) {
                    Some(hit) => {
                        hit.domains.push(dom);
                    }
                    None => {
                        hits.insert(
                            (key.0, key.1),
                            HmmerHit {
                                score,
                                e_value,
                                domains: vec![dom],
                            },
                        );
                    }
                }
            }

            Ok(Self {
                name: name.to_string(),
                hits,
            })
        }
    }
}

#[derive(Clone)]
pub struct Hit {
    pub query: String,
    pub target: String,
    // pub query_start: Option<usize>,
    // pub query_end: Option<usize>,
    // pub target_start: Option<usize>,
    // pub target_end: Option<usize>,
    pub score: f32,
    pub e_value: f64,
}

impl Display for Hit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {} {:.2e}", self.query, self.target, self.e_value)
    }
}

pub trait HitColumns {
    const QUERY: usize;
    const TARGET: usize;
    const Q_START: Option<usize>;
    const Q_END: Option<usize>;
    const T_START: Option<usize>;
    const T_END: Option<usize>;
    const SCORE: usize;
    const E_VALUE: usize;
}

// |      0      |    1     |     2      |     3    |    4   |  5   |  6  |    7    |   8  |  9  |  10 | 11| 12| 13| 14| 15| 16| 17|         18          |
// #                                                                             --- full sequence ---- --- best 1 domain ---- --- domain number estimation
// # target name   accession   query name  accession  E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
// #------------- ---------- ------------ ---------- -------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
pub struct HmmerTable {}
impl HitColumns for HmmerTable {
    const QUERY: usize = 2;
    const TARGET: usize = 0;
    const Q_START: Option<usize> = None;
    const Q_END: Option<usize> = None;
    const T_START: Option<usize> = None;
    const T_END: Option<usize> = None;
    const SCORE: usize = 5;
    const E_VALUE: usize = 4;
}

// |      0      |    1     |  2  |     3     |    4     |  5  |    6    |   7  |  8  | 9 | 10|   11    |    12   |  13  |  14 |  15 |  16 |  17 |  18 |  19 | 20  | 21 |         22          |
// #                                                             --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
// # target name  accession   tlen query name  accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
// #------------- ---------- ----- ----------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
pub struct HmmerDomTable {}
impl HitColumns for HmmerDomTable {
    const QUERY: usize = 3;
    const TARGET: usize = 0;
    const Q_START: Option<usize> = Some(15);
    const Q_END: Option<usize> = Some(16);
    const T_START: Option<usize> = Some(17);
    const T_END: Option<usize> = Some(18);
    const SCORE: usize = 13;
    const E_VALUE: usize = 12;
}

//  |   0   |   1  |  2   |   3  |  4  |  5  |  6  |  7 |   8  |  9  |
// #                target target query query       comp        cell
// # target  query  start  end    start end   score bias evalue frac
// # ------- ------ ------ ------ ----- ----- ----- ---- ------ -----
pub struct NailTable {}
impl HitColumns for NailTable {
    const QUERY: usize = 1;
    const TARGET: usize = 0;
    const Q_START: Option<usize> = Some(4);
    const Q_END: Option<usize> = Some(5);
    const T_START: Option<usize> = Some(2);
    const T_END: Option<usize> = Some(3);
    const SCORE: usize = 6;
    const E_VALUE: usize = 8;
}

//          |      0      |        1       |     2     |       3         |     4     |     5    |    6    |   7   |    8    |   9   |   10  |    11    |
// # Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
pub struct BlastTable {}
impl HitColumns for BlastTable {
    const QUERY: usize = 0;
    const TARGET: usize = 1;
    const Q_START: Option<usize> = Some(6);
    const Q_END: Option<usize> = Some(7);
    const T_START: Option<usize> = Some(8);
    const T_END: Option<usize> = Some(9);
    const SCORE: usize = 11;
    const E_VALUE: usize = 10;
}

pub struct HitTable {
    pub name: String,
    pub hits: Vec<Hit>,
}

impl HitTable {
    pub fn from_path<P: AsRef<Path>, C: HitColumns>(path: P) -> anyhow::Result<Self> {
        let path = path.as_ref();
        let file =
            File::open(path).with_context(|| format!("failed to open hit table: {path:?}"))?;
        let name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .context("invalid path")?;

        Self::parse::<File, C>(file, name)
    }

    pub fn parse<R: Read, C: HitColumns>(buf: R, name: &str) -> anyhow::Result<Self> {
        let reader = BufReader::new(buf);

        let mut hits = vec![];
        for line in reader.lines() {
            let line = line.unwrap_or_default();
            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let tokens = line.split_whitespace().collect::<Vec<_>>();

            hits.push(Hit {
                query: tokens[C::QUERY].to_string(),
                target: tokens[C::TARGET].to_string(),
                // query_start: C::Q_START.map(|i| tokens[i].parse::<usize>()).transpose()?,
                // query_end: C::Q_END.map(|i| tokens[i].parse::<usize>()).transpose()?,
                // target_start: C::T_START.map(|i| tokens[i].parse::<usize>()).transpose()?,
                // target_end: C::T_END.map(|i| tokens[i].parse::<usize>()).transpose()?,
                score: tokens[C::SCORE].parse()?,
                e_value: tokens[C::E_VALUE].parse()?,
            })
        }

        Ok(Self {
            name: name.to_string(),
            hits,
        })
    }

    pub fn from_path_filtered<P, C, F>(path: P, f: F) -> anyhow::Result<Self>
    where
        P: AsRef<Path>,
        C: HitColumns,
        F: Fn(&Hit) -> bool,
    {
        let path = path.as_ref();
        let file =
            File::open(path).with_context(|| format!("failed to open hit table: {path:?}"))?;
        let name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .context("invalid path")?;

        Self::parse_filtered::<File, C, F>(file, f, name)
    }

    pub fn parse_filtered<R, C, F>(buf: R, f: F, name: &str) -> anyhow::Result<Self>
    where
        R: Read,
        C: HitColumns,
        F: Fn(&Hit) -> bool,
    {
        let reader = BufReader::new(buf);

        let mut hits = vec![];
        for line in reader.lines() {
            let line = line.unwrap_or_default();
            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let tokens = line.split_whitespace().collect::<Vec<_>>();
            let hit = Hit {
                query: tokens[C::QUERY].to_string(),
                target: tokens[C::TARGET].to_string(),
                score: tokens[C::SCORE].parse()?,
                e_value: tokens[C::E_VALUE].parse()?,
            };

            if f(&hit) {
                hits.push(hit);
            }
        }

        Ok(Self {
            name: name.to_string(),
            hits,
        })
    }

    pub fn to_map(self) -> HashMap<(String, String), Hit> {
        self.hits
            .into_iter()
            .map(|h| ((h.query.clone(), h.target.clone()), h))
            .collect()
    }

    pub fn to_query_map(self) -> HashMap<String, Vec<Hit>> {
        let mut map: HashMap<String, Vec<Hit>> = HashMap::new();
        self.hits
            .into_iter()
            .for_each(|h| map.entry(h.query.clone()).or_default().push(h));
        map
    }

    pub fn to_target_map(self) -> HashMap<String, Vec<Hit>> {
        let mut map: HashMap<String, Vec<Hit>> = HashMap::new();
        self.hits
            .into_iter()
            .for_each(|h| map.entry(h.target.clone()).or_default().push(h));
        map
    }
}

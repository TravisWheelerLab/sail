//! Which records a sequential sample keeps, as arithmetic over counts.
//!
//! [`Selection`] answers "which `m` of `n`" in ascending order and constant
//! memory, and never reads anything. That is the same split as [`Delimit`],
//! which marks record boundaries without seeing bytes: the choice is testable
//! with no file behind it, and crossing the gap between one selected record and
//! the next is the caller's to make — a scan on a stream, a seek on an index.
//!
//! [`Delimit`]: crate::index::Delimit

// --- the generator

/// SplitMix64's golden-ratio increment.
const GAMMA: u64 = 0x9e37_79b9_7f4a_7c15;

/// Stafford's Mix13 multipliers, whose shift amounts came out of the same
/// search, so do not change one without the others.
const MIX_A: u64 = 0xbf58_476d_1ce4_e5b9;
const MIX_B: u64 = 0x94d0_49bb_1331_11eb;

/// Scales a 53-bit integer into the unit interval.
const UNIT: f64 = 1.0 / (1u64 << 53) as f64;

// TODO: feisty keys its permutation with this same
//       generator and is about to make it public. swap this
//       for the export once it lands, so the two combinators
//       draw from one implementation rather than two copies.
//       the constants are identical, so the swap changes no
//       draw and invalidates no recorded seed
fn step(state: &mut u64) -> u64 {
    *state = state.wrapping_add(GAMMA);

    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(MIX_A);
    z = (z ^ (z >> 27)).wrapping_mul(MIX_B);

    z ^ (z >> 31)
}

// ---

/// `m` of `0..n` drawn uniformly at random without replacement, ascending, in
/// constant memory.
///
/// Yields exactly `m` indices, every one below `n`, each larger than the last.
/// All `C(n, m)` subsets are equally likely.
///
/// # Panics
///
/// If `m > n`.
pub struct Selection {
    /// How many indices are still to be yielded.
    want: usize,

    /// How many of `0..n` have not yet been passed over or yielded.
    left: usize,

    /// Where the next skip is counted from.
    at: usize,

    state: u64,
}

impl Selection {
    /// `m` of `0..n`, drawn under `seed`.
    pub fn new(n: usize, m: usize, seed: u64) -> Self {
        assert!(
            m <= n,
            "Selection::new({n}, {m}) -- a draw without replacement cannot \
             return more indices than there are"
        );

        Selection {
            want: m,
            left: n,
            at: 0,
            // a caller handing over 0 is handing over a
            // number, not an absence, and splitmix64's step
            // is what separates it from any other seed
            state: seed,
        }
    }

    /// A uniform draw from `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (step(&mut self.state) >> 11) as f64 * UNIT
    }

    /// A uniform draw from `(0, 1]`, for the calls that take its logarithm.
    fn positive(&mut self) -> f64 {
        ((step(&mut self.state) >> 11) + 1) as f64 * UNIT
    }
}

impl Iterator for Selection {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        if self.want == 0 {
            return None;
        }

        // Vitter's own switch: D pays several transcendentals
        // per draw to skip in O(1), which only wins while the
        // gaps are long. alphainv = 13 is the crossover he
        // reports in section 4
        const ALPHA_INV: usize = 13;

        let skip = if self.want == 1 {
            // one record out of what is left, uniformly
            (self.left as f64 * self.unit()) as usize
        } else if ALPHA_INV * self.want < self.left {
            self.skip_d()
        } else {
            self.skip_a()
        };

        let n = self.at + skip;

        self.at = n + 1;
        self.left -= skip + 1;
        self.want -= 1;

        Some(n)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.want, Some(self.want))
    }
}

impl ExactSizeIterator for Selection {}

impl Selection {
    /// How many records to pass over before the next one to keep, by Algorithm
    /// A: one uniform draw and a walk down the remaining records.
    //
    // Vitter 1987, "An Efficient Algorithm for Sequential
    // Random Sampling", ACM TOMS 13(1), Algorithm A. the
    // inner walk costs one multiply per record passed over,
    // which is why D exists and why this only runs once the
    // gaps are short
    fn skip_a(&mut self) -> usize {
        let v = self.unit();

        let mut left = self.left as f64;
        let mut top = (self.left - self.want) as f64;
        let mut quot = top / left;
        let mut skip = 0;

        while quot > v {
            skip += 1;
            top -= 1.0;
            left -= 1.0;
            quot = quot * top / left;
        }

        skip
    }

    /// How many records to pass over before the next one to keep, by Method D:
    /// a skip drawn directly, accepted or redrawn.
    //
    // Vitter 1987, Algorithm D. one divergence: he carries
    // V' across calls, where this redraws it at the top of
    // every skip. the value he carries is distributed as
    // U^(1/n) for the next call's n, which is exactly what
    // redrawing produces, so the draws are the same ones --
    // it costs a log and an exp per selected record and
    // leaves next() with no state to get wrong
    fn skip_d(&mut self) -> usize {
        let n = self.want as f64;
        let left = self.left as f64;

        // the largest skip that still leaves want - 1
        // records to choose from
        let qu1 = (self.left - self.want + 1) as f64;
        let nmin1inv = 1.0 / (n - 1.0);

        let mut vprime = (self.positive().ln() / n).exp();

        loop {
            let (x, skip) = loop {
                let x = left * (1.0 - vprime);
                let skip = x.floor();

                if skip < qu1 {
                    break (x, skip);
                }

                vprime = (self.positive().ln() / n).exp();
            };

            let y1 = ((self.positive() * left / qu1).ln() * nmin1inv).exp();
            vprime = y1 * (1.0 - x / left) * (qu1 / (qu1 - skip));

            if vprime <= 1.0 {
                return skip as usize;
            }

            // the expensive test, reached only when the cheap
            // one was inconclusive: a falling factorial over
            // min(skip, want - 1) terms
            let mut y2 = 1.0;
            let mut top = left - 1.0;
            let (mut bottom, limit) = if n - 1.0 > skip {
                (left - n, left - skip)
            } else {
                (left - skip - 1.0, qu1)
            };

            let mut t = left - 1.0;
            while t >= limit {
                y2 = y2 * top / bottom;
                top -= 1.0;
                bottom -= 1.0;
                t -= 1.0;
            }

            if left / (left - x) >= y1 * (y2.ln() * nmin1inv).exp() {
                return skip as usize;
            }

            vprime = (self.positive().ln() / n).exp();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// Every (n, m) pair worth covering, small enough to check exhaustively.
    //
    // the pairs straddle the alpha crossover deliberately:
    // 13 * m < n takes Method D and anything else takes
    // Algorithm A, and a test that only ever hit one of them
    // would pass with the other unwritten
    const CASES: [(usize, usize); 10] = [
        (1, 0),
        (1, 1),
        (10, 1),
        (10, 5),
        (10, 10),
        (200, 5),
        (200, 100),
        (200, 200),
        (5000, 3),
        (5000, 4999),
    ];

    #[test]
    fn a_draw_yields_exactly_what_was_asked_for_ascending_and_in_range() {
        for (n, m) in CASES {
            for seed in 0..50u64 {
                let drawn: Vec<usize> = Selection::new(n, m, seed).collect();

                assert_eq!(drawn.len(), m, "n={n} m={m} seed={seed}");
                assert!(drawn.iter().all(|&i| i < n), "n={n} m={m} seed={seed}");
                assert!(
                    drawn.windows(2).all(|w| w[0] < w[1]),
                    "n={n} m={m} seed={seed}: {drawn:?}"
                );
            }
        }
    }

    #[test]
    fn asking_for_the_whole_population_yields_every_index_in_order() {
        for seed in 0..10u64 {
            let drawn: Vec<usize> = Selection::new(64, 64, seed).collect();
            assert_eq!(drawn, (0..64).collect::<Vec<_>>());
        }
    }

    #[test]
    fn an_empty_draw_and_an_empty_population_both_yield_nothing() {
        assert_eq!(Selection::new(100, 0, 1).count(), 0);
        assert_eq!(Selection::new(0, 0, 1).count(), 0);
    }

    #[test]
    fn the_same_seed_draws_the_same_indices_and_another_seed_does_not() {
        let drawn = |seed| Selection::new(1000, 20, seed).collect::<Vec<_>>();

        assert_eq!(drawn(7), drawn(7));
        assert_ne!(drawn(7), drawn(8));
    }

    #[test]
    fn the_length_is_known_before_the_draw_is_walked() {
        let s = Selection::new(1000, 20, 1);

        assert_eq!(s.len(), 20);
        assert_eq!(s.size_hint(), (20, Some(20)));
    }

    #[test]
    fn drawing_more_indices_than_there_are_panics() {
        assert!(std::panic::catch_unwind(|| Selection::new(3, 4, 1)).is_err());
    }

    #[test]
    fn every_subset_of_the_same_size_comes_up_about_as_often() {
        // 13 * 3 >= 8, so this is the Algorithm A path, and
        // C(8, 3) = 56 subsets is few enough to count them
        // all rather than settle for a marginal
        let (n, m, trials) = (8usize, 3usize, 56_000u64);

        let mut counts: HashMap<Vec<usize>, u32> = HashMap::new();
        for seed in 0..trials {
            *counts
                .entry(Selection::new(n, m, seed).collect())
                .or_default() += 1;
        }

        assert_eq!(counts.len(), 56, "every subset should be reachable");

        // 1000 expected each; a quarter either way is far
        // wider than sampling noise and far narrower than
        // any real bias would be
        let (lo, hi) = counts
            .values()
            .fold((u32::MAX, 0), |(lo, hi), &c| (lo.min(c), hi.max(c)));

        assert!(lo > 750, "least common subset came up {lo} times");
        assert!(hi < 1250, "most common subset came up {hi} times");
    }

    #[test]
    fn every_position_is_drawn_about_as_often_on_the_skipping_path() {
        // 13 * 5 < 200, so this is Method D. C(200, 5) is far
        // too large to enumerate, so this checks the marginal
        // instead: each position should turn up in m/n of the
        // draws
        let (n, m, trials) = (200usize, 5usize, 40_000u64);

        let mut counts = vec![0u32; n];
        for seed in 0..trials {
            for i in Selection::new(n, m, seed) {
                counts[i] += 1;
            }
        }

        // 40,000 * 5 / 200 = 1000 expected per position
        let lo = *counts.iter().min().unwrap();
        let hi = *counts.iter().max().unwrap();

        assert!(lo > 800, "least drawn position came up {lo} times");
        assert!(hi < 1200, "most drawn position came up {hi} times");
    }

    #[test]
    fn the_first_and_last_positions_are_no_less_likely_than_the_middle() {
        // a sequential sampler that gets its conditioning
        // wrong biases toward one end of the file, which a
        // whole-population average hides
        let (n, m, trials) = (1000usize, 10usize, 20_000u64);

        let mut counts = vec![0u32; n];
        for seed in 0..trials {
            for i in Selection::new(n, m, seed) {
                counts[i] += 1;
            }
        }

        let band = |r: std::ops::Range<usize>| counts[r].iter().sum::<u32>();
        let (front, middle, back) = (band(0..100), band(450..550), band(900..1000));

        // 20,000 * 10 / 1000 * 100 = 20,000 expected per band
        for (what, got) in [("front", front), ("middle", middle), ("back", back)] {
            assert!(
                (18_000..22_000).contains(&got),
                "the {what} of the file was drawn {got} times, expected ~20000"
            );
        }
    }
}

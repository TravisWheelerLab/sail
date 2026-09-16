//! Collections of records addressed by position.
//!
//! [`Indexable`] is the whole of what an implementor writes — a length and a
//! `get` — and the rest is built on it: [`Iterable`] by blanket impl, and
//! `filter`/`sort_by`/`take` through [`View`], which reorders a permutation of
//! indices and leaves the records themselves alone.
//!
//! [`Drawn`] is the other reordering. `permute` and `sample` compute their
//! permutation from a seed instead of storing one, so they allocate the same
//! amount at any length.
//!
//! [`Aggregate`] is the last adapter: several collections addressed as one.

use std::borrow::Borrow;
use std::collections::HashMap;
use std::hash::Hash;

use feisty::Permutation;
use indexmap::IndexMap;

use crate::Record;

// ---

/// A reordered or narrowed view of another collection.
pub struct View<I>
where
    I: Indexable,
{
    inner: I,

    // a permutation of indices into inner, so filter,
    // permute and the sorts cost one usize per surviving
    // record and never move a record
    indices: Vec<usize>,
}

impl<I> Indexable for View<I>
where
    I: Indexable,
{
    type Record = I::Record;
    type Yield<'a>
        = I::Yield<'a>
    where
        Self: 'a;

    fn len(&self) -> usize {
        self.indices.len()
    }

    fn get(&self, n: usize) -> Option<Self::Yield<'_>> {
        self.indices.get(n).and_then(|&x| self.inner.get(x))
    }
}

/// Another collection read through a permutation of `0..len()`, and how much of
/// that permutation this holds.
pub struct Drawn<I>
where
    I: Indexable,
{
    inner: I,

    /// The permutation, computed from a seed rather than stored.
    //
    // None only where inner is empty: feistel networks are
    // built over a non-empty domain, and an empty
    // collection has no position to map
    perm: Option<Permutation<usize>>,

    /// How many of the permuted positions this holds, counting from the front.
    //
    // a prefix, so two draws from one collection under one
    // seed nest: the smaller is a subset of the larger
    len: usize,
}

impl<I> Drawn<I>
where
    I: Indexable,
{
    /// The first `len` positions of `seed`'s permutation of `inner`.
    fn new(inner: I, len: usize, seed: u64) -> Self {
        let n = inner.len();

        Drawn {
            inner,
            perm: (n > 0).then(|| Permutation::new(n, seed)),
            len,
        }
    }
}

impl<I> Indexable for Drawn<I>
where
    I: Indexable,
{
    type Record = I::Record;
    type Yield<'a>
        = I::Yield<'a>
    where
        Self: 'a;

    fn len(&self) -> usize {
        self.len
    }

    fn get(&self, n: usize) -> Option<Self::Yield<'_>> {
        // Permutation::get asserts its argument is inside
        // the domain, so the bound is checked here rather
        // than left to panic out of the mapping
        if n >= self.len {
            return None;
        }

        self.inner.get(self.perm.as_ref()?.get(n))
    }
}

/// A shared reference to a collection is one too, so the combinators — which
/// all take `self` — can be pointed at a collection without consuming it.
impl<I> Indexable for &I
where
    I: Indexable + ?Sized,
{
    type Record = I::Record;
    type Yield<'a>
        = I::Yield<'a>
    where
        Self: 'a;

    fn len(&self) -> usize {
        (**self).len()
    }

    fn get(&self, n: usize) -> Option<Self::Yield<'_>> {
        (**self).get(n)
    }
}

/// Fetch record `n`, which [`Indexable::len`] has already counted.
fn expect_record<'a, I>(c: &'a I, n: usize) -> I::Yield<'a>
where
    I: Indexable + ?Sized,
{
    // a None here is a broken invariant and not a value:
    // skipping drops a record silently, and stopping
    // truncates the traversal where no caller can see it
    c.get(n).unwrap_or_else(|| {
        panic!(
            "Indexable::get({n}) returned None though len() is {} -- a \
             collection must be able to produce every record it counts",
            c.len()
        )
    })
}

// ---

/// A collection addressed by position. An implementor writes `len` and `get`,
/// and **`get(i)` must be `Some` for every `i < len()`** or the rest panics.
pub trait Indexable {
    /// The record itself, owned.
    type Record: Record;

    /// What [`get`](Self::get) hands back — `&Record` from a collection holding
    /// its records, `Record` from one building each on demand.
    //
    // the lifetime has to sit here and nowhere else. on the
    // trait instead, `filter(self) -> View<Self>` would
    // borrow from a collection it is moving, and no extra
    // associated type can introduce one
    //
    // Borrow is what keeps the two backends one API: std
    // gives both `impl Borrow<T> for &T` and
    // `impl Borrow<T> for T`, so a caller writes
    // `Indexable<Record = FastaRecord>` either way
    type Yield<'a>: Borrow<Self::Record>
    where
        Self: 'a;

    /// The number of records held.
    fn len(&self) -> usize;

    /// Whether the collection holds no records.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The `n`th record, or `None` past the end.
    //
    // an implementor that holds its records lends one; one
    // that parses from a file at a known offset has nothing
    // to lend and yields it owned. that second case is what
    // makes a call expensive, and why the sorts reach for a
    // record once per record rather than once per comparison
    fn get(&self, n: usize) -> Option<Self::Yield<'_>>;

    /// The `n`th record owned, cloning only where [`get`](Self::get) lends.
    fn cloned(&self, n: usize) -> Option<Self::Record> {
        self.get(n).map(|r| r.borrow().clone())
    }

    /// The index of the first record satisfying `predicate`.
    fn position<P>(&self, mut predicate: P) -> Option<usize>
    where
        P: FnMut(&Self::Record) -> bool,
    {
        (0..self.len()).find(|&n| predicate(expect_record(self, n).borrow()))
    }

    /// The first record satisfying `predicate`.
    fn find<P>(&self, mut predicate: P) -> Option<Self::Yield<'_>>
    where
        P: FnMut(&Self::Record) -> bool,
    {
        // one pass rather than position() then get(), which
        // would reach for the found record twice -- and on a
        // disk-backed collection that is a second read
        (0..self.len())
            .map(|n| expect_record(self, n))
            .find(|r| predicate(r.borrow()))
    }

    /// The record whose key is smallest, or `None` if there are none.
    fn min_by_key<K, F>(&self, mut f: F) -> Option<Self::Yield<'_>>
    where
        K: Ord,
        F: FnMut(&Self::Record) -> K,
    {
        let n = (0..self.len()).min_by_key(|&n| f(expect_record(self, n).borrow()))?;

        self.get(n)
    }

    /// The record whose key is largest, or `None` if there are none.
    fn max_by_key<K, F>(&self, mut f: F) -> Option<Self::Yield<'_>>
    where
        K: Ord,
        F: FnMut(&Self::Record) -> K,
    {
        let n = (0..self.len()).max_by_key(|&n| f(expect_record(self, n).borrow()))?;

        self.get(n)
    }

    /// A map from each record's key to its index, for lookup by name rather
    /// than by position.
    fn index_by<K, F>(&self, mut key: F) -> HashMap<K, usize>
    where
        K: Eq + Hash,
        F: FnMut(&Self::Record) -> K,
    {
        let mut map = HashMap::with_capacity(self.len());

        for n in 0..self.len() {
            // or_insert rather than insert, so a repeated key
            // keeps the first record's index: a duplicate
            // identifier means a malformed file, and
            // returning the last would answer a lookup with a
            // record the file does not lead with
            map.entry(key(expect_record(self, n).borrow())).or_insert(n);
        }

        map
    }

    /// `m` records drawn uniformly at random without replacement, in the order
    /// `seed`'s permutation puts them.
    ///
    /// # Panics
    ///
    /// If `m > len()`. A draw without replacement cannot return more records
    /// than there are, and asking for more is a bug rather than a request to
    /// clamp.
    fn sample(self, m: usize, seed: u64) -> Drawn<Self>
    where
        Self: Sized,
    {
        assert!(
            m <= self.len(),
            "sample({m}) from a collection of {} -- a draw without replacement \
             cannot return more records than there are",
            self.len()
        );

        Drawn::new(self, m, seed)
    }

    /// The records for which `predicate` holds, in their original order.
    fn filter<P>(self, mut predicate: P) -> View<Self>
    where
        Self: Sized,
        P: FnMut(&Self::Record) -> bool,
    {
        let indices = (0..self.len())
            .filter(|&i| predicate(expect_record(&self, i).borrow()))
            .collect::<Vec<_>>();

        View {
            inner: self,
            indices,
        }
    }

    /// The first `n` records, or all of them if there are fewer.
    fn take(self, n: usize) -> View<Self>
    where
        Self: Sized,
    {
        let indices = (0..n.min(self.len())).collect();

        View {
            inner: self,
            indices,
        }
    }

    /// Every record after the first `n`.
    fn skip(self, n: usize) -> View<Self>
    where
        Self: Sized,
    {
        let indices = (n.min(self.len())..self.len()).collect();

        View {
            inner: self,
            indices,
        }
    }

    /// The same records in the order `seed`'s permutation puts them, which is
    /// [`sample`](Self::sample) of all of them.
    fn permute(self, seed: u64) -> Drawn<Self>
    where
        Self: Sized,
    {
        let len = self.len();

        Drawn::new(self, len, seed)
    }

    /// The records ordered by `cmp`, materializing both at every comparison.
    /// Prefer [`sort_by_key`](Self::sort_by_key) for an order from one value.
    fn sort_by<F>(self, mut cmp: F) -> View<Self>
    where
        Self: Sized,
        F: FnMut(&Self::Record, &Self::Record) -> std::cmp::Ordering,
    {
        let mut indices: Vec<usize> = (0..self.len()).collect();

        let cmp2 = |a: &usize, b: &usize| {
            let x = expect_record(&self, *a);
            let y = expect_record(&self, *b);

            cmp(x.borrow(), y.borrow())
        };

        indices.sort_by(cmp2);

        View {
            inner: self,
            indices,
        }
    }

    /// The records satisfying `predicate` and those that do not, in one pass.
    fn partition<P>(&self, mut predicate: P) -> (View<&Self>, View<&Self>)
    where
        Self: Sized,
        P: FnMut(&Self::Record) -> bool,
    {
        // by reference, unlike the combinators above: two
        // views of one collection cannot each own it, and
        // &Self is itself Indexable so they borrow it
        let (yes, no) = (0..self.len()).partition(|&n| predicate(expect_record(self, n).borrow()));

        (
            View {
                inner: self,
                indices: yes,
            },
            View {
                inner: self,
                indices: no,
            },
        )
    }

    /// Views of `n` consecutive records each, the last shorter if it has to be.
    /// An empty collection yields none; [`split`](Self::split) always yields `k`.
    fn chunks(&self, n: usize) -> Vec<View<&Self>>
    where
        Self: Sized,
    {
        let n = n.max(1);

        (0..self.len())
            .step_by(n)
            .map(|start| View {
                inner: self,
                indices: (start..(start + n).min(self.len())).collect(),
            })
            .collect()
    }

    /// `k` views of near-equal size, for handing a collection to `k` workers.
    /// Always exactly `k`, even where that means empty ones.
    fn split(&self, k: usize) -> Vec<View<&Self>>
    where
        Self: Sized,
    {
        let k = k.max(1);

        // the remainder goes one record each to the parts in
        // front, so every part is within one of every other
        let (size, remainder) = (self.len() / k, self.len() % k);

        let mut parts = Vec::with_capacity(k);
        let mut start = 0;

        for part in 0..k {
            let n = size + usize::from(part < remainder);
            parts.push(View {
                inner: self,
                indices: (start..start + n).collect(),
            });
            start += n;
        }

        parts
    }

    /// `k` views of near-equal total weight, for `k` workers whose cost is not
    /// the record count. Always exactly `k`, as [`split`](Self::split) is.
    fn split_weighted<F>(&self, k: usize, mut weight: F) -> Vec<View<&Self>>
    where
        Self: Sized,
        F: FnMut(&Self::Record) -> u64,
    {
        let k = k.max(1);

        // materialized once per record rather than once per
        // comparison, because get may be a read and a parse
        // on a disk-backed collection
        let weights: Vec<u64> = (0..self.len())
            .map(|n| weight(expect_record(self, n).borrow()))
            .collect();

        // longest-processing-time: the heaviest record goes
        // first, onto the lightest part so far, which bounds
        // the worst part at 4/3 of optimal. the sort is
        // stable, so equal weights stay in file order and one
        // input always deals the same way
        let mut order: Vec<usize> = (0..self.len()).collect();
        order.sort_by(|&a, &b| weights[b].cmp(&weights[a]));

        let mut indices = vec![Vec::new(); k];
        let mut loads = vec![0u64; k];

        for n in order {
            // min_by_key reports the first minimum, so a tie
            // goes to the lowest-numbered part
            let lightest = (0..k).min_by_key(|&p| loads[p]).expect("k is at least 1");

            loads[lightest] += weights[n];
            indices[lightest].push(n);
        }

        indices
            .into_iter()
            .map(|mut part| {
                // the deal ordered these by weight, and a
                // caller writing a part out wants the order
                // the file holds them in
                part.sort_unstable();

                View {
                    inner: self,
                    indices: part,
                }
            })
            .collect()
    }

    /// The records gathered by key, each group in file order and the groups in
    /// the order their keys first appear.
    fn group_by<K, F>(&self, mut key: F) -> Vec<(K, View<&Self>)>
    where
        Self: Sized,
        K: Eq + Hash,
        F: FnMut(&Self::Record) -> K,
    {
        // an IndexMap and not a HashMap, so grouping hits
        // by query reports queries in the order the table
        // lists them rather than however they hash
        let mut groups: IndexMap<K, Vec<usize>> = IndexMap::new();

        for n in 0..self.len() {
            groups
                .entry(key(expect_record(self, n).borrow()))
                .or_default()
                .push(n);
        }

        groups
            .into_iter()
            .map(|(k, indices)| {
                (
                    k,
                    View {
                        inner: self,
                        indices,
                    },
                )
            })
            .collect()
    }

    /// The records ordered by a key extracted from each.
    fn sort_by_key<K, F>(self, mut f: F) -> View<Self>
    where
        Self: Sized,
        K: Ord,
        F: FnMut(&Self::Record) -> K,
    {
        let mut indices: Vec<usize> = (0..self.len()).collect();

        // every key up front, so get() is called once per
        // record rather than twice per comparison -- n
        // materializations instead of O(n log n)
        let keys = (0..self.len())
            .map(|i| f(expect_record(&self, i).borrow()))
            .collect::<Vec<_>>();

        let cmp = |a: &usize, b: &usize| keys[*a].cmp(&keys[*b]);

        indices.sort_by(cmp);

        View {
            inner: self,
            indices,
        }
    }
}

// ---

/// A collection that can be walked front to back.
pub trait Iterable: Indexable {
    /// Every record in order, yielded exactly as [`get`](Indexable::get) does.
    //
    // Indexable is a supertrait so this can name Yield. as a
    // standalone trait the item would only be known to
    // Borrow the record, and every call site would deref
    fn iter(&self) -> impl Iterator<Item = Self::Yield<'_>>;
}

/// An iterator over an [`Indexable`], by position.
pub struct Iter<'a, O>
where
    O: Indexable + ?Sized,
{
    inner: &'a O,
    pos: usize,
}

impl<'a, T> Iterator for Iter<'a, T>
where
    T: Indexable + ?Sized,
{
    type Item = T::Yield<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        // stop at len() rather than at the first None: an
        // impl that breaks the contract and keeps returning
        // Some then gives a wrong count, which is debuggable,
        // instead of looping forever, which is not
        if self.pos >= self.inner.len() {
            return None;
        }

        let n = self.pos;
        self.pos += 1;

        Some(expect_record(self.inner, n))
    }
}

impl<T> Iterable for T
where
    T: Indexable,
{
    fn iter(&self) -> impl Iterator<Item = Self::Yield<'_>> {
        Iter {
            inner: self,
            pos: 0,
        }
    }
}

// ---

/// Several collections addressed as one, in the order they were added.
pub struct Aggregate<I>
where
    I: Indexable,
{
    parts: Vec<I>,

    // global index of each part's first record, with the
    // total on the end. always one longer than parts, so
    // starts[i + 1] - starts[i] is part i's length
    //
    // one usize per *part* and not per record, binary
    // searched to dispatch a get: two 50-record files cost
    // 24 bytes of index between them
    //
    // cached rather than asked for on every get, which is
    // sound because a collection that counted differently
    // from one call to the next would already have broken
    // Indexable
    starts: Vec<usize>,
}

impl<I> Aggregate<I>
where
    I: Indexable,
{
    /// Address `parts` as one collection.
    pub fn new(parts: impl IntoIterator<Item = I>) -> Self {
        let mut agg = Aggregate::default();
        for part in parts {
            agg.push(part);
        }

        agg
    }

    /// Add a collection after the ones already held.
    pub fn push(&mut self, part: I) {
        let total = self.len() + part.len();
        self.parts.push(part);
        self.starts.push(total);
    }

    /// The collections held, in order.
    pub fn parts(&self) -> &[I] {
        &self.parts
    }

    /// Which part holds record `n`, and where in it — `None` past the end.
    /// The way to say *which file* a record came from.
    pub fn locate(&self, n: usize) -> Option<(usize, usize)> {
        if n >= self.len() {
            return None;
        }

        // a scan would do, but starts is sorted and this
        // stays O(log parts) for a directory of thousands.
        // it cannot underflow: starts[0] is 0 and n is
        // already known to be below the total
        //
        // an empty part repeats the previous start, so this
        // takes the *last* start at or before n and steps
        // over it
        let part = self.starts.partition_point(|&start| start <= n) - 1;

        Some((part, n - self.starts[part]))
    }
}

impl<I> Default for Aggregate<I>
where
    I: Indexable,
{
    fn default() -> Self {
        Aggregate {
            parts: Vec::new(),
            starts: vec![0],
        }
    }
}

impl<I> FromIterator<I> for Aggregate<I>
where
    I: Indexable,
{
    fn from_iter<T: IntoIterator<Item = I>>(parts: T) -> Self {
        Aggregate::new(parts)
    }
}

impl<I> Indexable for Aggregate<I>
where
    I: Indexable,
{
    type Record = I::Record;

    // whatever the parts yield: an aggregate of Fasta lends,
    // an aggregate of IndexedFasta owns
    type Yield<'a>
        = I::Yield<'a>
    where
        Self: 'a;

    fn len(&self) -> usize {
        self.starts.last().copied().unwrap_or(0)
    }

    fn get(&self, n: usize) -> Option<Self::Yield<'_>> {
        let (part, local) = self.locate(n)?;

        self.parts[part].get(local)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deliberately breaks the len/get contract: claims three records and hands
    /// one back for any index at all.
    struct Liar;

    #[derive(Clone, Debug, PartialEq)]
    struct Unit;

    impl Record for Unit {}

    impl Indexable for Liar {
        type Record = Unit;
        type Yield<'a> = Unit;

        fn len(&self) -> usize {
            3
        }

        fn get(&self, _: usize) -> Option<Unit> {
            Some(Unit)
        }
    }

    #[test]
    fn iteration_stops_at_len_even_if_get_never_says_no() {
        // without the bound in Iter::next this hangs
        // rather than fails
        assert_eq!(Liar.iter().count(), 3);
    }

    #[test]
    fn an_empty_collection_reports_itself_empty() {
        struct None_;
        impl Indexable for None_ {
            type Record = Unit;
            type Yield<'a> = Unit;
            fn len(&self) -> usize {
                0
            }
            fn get(&self, _: usize) -> Option<Unit> {
                None
            }
        }

        assert!(None_.is_empty());
        assert_eq!(None_.iter().count(), 0);
    }

    /// Counts three records and refuses to produce the middle one, the way a
    /// disk-backed collection behaves when a read fails.
    struct Gappy;

    impl Indexable for Gappy {
        type Record = Unit;
        type Yield<'a> = Unit;

        fn len(&self) -> usize {
            3
        }

        fn get(&self, n: usize) -> Option<Unit> {
            (n != 1).then_some(Unit)
        }
    }

    #[test]
    fn a_record_that_cannot_be_produced_is_loud_rather_than_missing() {
        // the panic has to name the index and the length;
        // silently ending at record 0 loses records 2..
        let err = std::panic::catch_unwind(|| Gappy.iter().count()).unwrap_err();
        let msg = err
            .downcast_ref::<String>()
            .expect("a formatted panic message");

        assert!(msg.contains("get(1)"), "{msg}");
        assert!(msg.contains("len() is 3"), "{msg}");
    }

    /// Holds its records, so it lends one.
    struct Lender(Vec<Unit>);

    impl Indexable for Lender {
        type Record = Unit;
        type Yield<'a> = &'a Unit;

        fn len(&self) -> usize {
            self.0.len()
        }

        fn get(&self, n: usize) -> Option<&Unit> {
            self.0.get(n)
        }
    }

    #[test]
    fn a_collection_holding_its_records_lends_one_and_a_parsing_one_owns_it() {
        // the whole point of Yield being separate from
        // Record. no other test here would catch a backend
        // changed back to cloning, because Borrow hides
        // which it is at every call site that uses one
        fn lends<'a, C: Indexable<Yield<'a> = &'a Unit>>(_: &'a C) {}
        fn owns<C: for<'a> Indexable<Yield<'a> = Unit>>(_: &C) {}

        lends(&Lender(vec![Unit]));
        owns(&Liar);

        // and the one bound both satisfy, which is what keeps
        // this a single API rather than two
        fn either<C: Indexable<Record = Unit>>(c: &C) -> usize {
            c.iter()
                .map(|r| {
                    let _: &Unit = r.borrow();
                    1
                })
                .sum()
        }

        assert_eq!(either(&Lender(vec![Unit, Unit])), 2);
        assert_eq!(either(&Liar), 3);
    }

    #[test]
    fn cloned_hands_back_an_owned_record_from_either_kind() {
        assert_eq!(Lender(vec![Unit]).cloned(0), Some(Unit));
        assert_eq!(Lender(vec![Unit]).cloned(9), None);
        assert_eq!(Liar.cloned(0), Some(Unit));
    }

    #[test]
    fn the_combinators_are_loud_about_it_too() {
        for (what, r) in [
            (
                "filter",
                std::panic::catch_unwind(|| Gappy.filter(|_| true).len()),
            ),
            // sample chooses indices and touches no record,
            // so the walk is what has to be loud
            (
                "sample",
                std::panic::catch_unwind(|| Gappy.sample(3, 0).iter().count()),
            ),
            (
                "sort_by_key",
                std::panic::catch_unwind(|| Gappy.sort_by_key(|_| 0).len()),
            ),
        ] {
            let Err(err) = r else {
                panic!("{what} should panic");
            };
            let msg = err.downcast_ref::<String>().expect("a formatted message");
            assert!(msg.contains("must be able to produce"), "{what}: {msg}");
        }
    }

    // ---

    /// Holds `n` records, numbered from `base` so a record says which part it
    /// came from.
    struct Part {
        base: u32,
        n: usize,
    }

    #[derive(Clone, Debug, PartialEq)]
    struct Num(u32);

    impl Record for Num {}

    impl Indexable for Part {
        type Record = Num;
        type Yield<'a> = Num;

        fn len(&self) -> usize {
            self.n
        }

        fn get(&self, n: usize) -> Option<Num> {
            (n < self.n).then(|| Num(self.base + n as u32))
        }
    }

    fn part(base: u32, n: usize) -> Part {
        Part { base, n }
    }

    #[test]
    fn an_aggregate_is_as_long_as_its_parts_together() {
        let agg = Aggregate::new([part(0, 50), part(100, 50)]);

        assert_eq!(agg.len(), 100);
        assert!(!agg.is_empty());
        assert_eq!(agg.parts().len(), 2);
    }

    #[test]
    fn an_index_lands_in_the_part_that_holds_it() {
        // the boundary is the whole of this: 49 is the last
        // of the first part and 50 is the first of the second
        let agg = Aggregate::new([part(0, 50), part(100, 50)]);

        assert_eq!(agg.get(0), Some(Num(0)));
        assert_eq!(agg.get(49), Some(Num(49)));
        assert_eq!(agg.get(50), Some(Num(100)), "the second part restarts");
        assert_eq!(agg.get(99), Some(Num(149)));
        assert_eq!(agg.get(100), None, "past the end is None, not a panic");
    }

    #[test]
    fn locate_says_which_part_a_record_came_from() {
        let agg = Aggregate::new([part(0, 50), part(100, 50)]);

        assert_eq!(agg.locate(0), Some((0, 0)));
        assert_eq!(agg.locate(49), Some((0, 49)));
        assert_eq!(agg.locate(50), Some((1, 0)));
        assert_eq!(agg.locate(99), Some((1, 49)));
        assert_eq!(agg.locate(100), None);
    }

    #[test]
    fn an_empty_part_is_stepped_over_rather_than_addressed() {
        // an empty part repeats the previous start, so the
        // search has to take the *last* start at or before n
        // or it dispatches into a part holding nothing
        let agg = Aggregate::new([part(0, 2), part(100, 0), part(200, 2)]);

        assert_eq!(agg.len(), 4);
        assert_eq!(agg.get(1), Some(Num(1)));
        assert_eq!(agg.get(2), Some(Num(200)), "not the empty part");
        assert_eq!(agg.locate(2), Some((2, 0)));
        assert_eq!(
            agg.iter().collect::<Vec<_>>(),
            [Num(0), Num(1), Num(200), Num(201)]
        );
    }

    #[test]
    fn empty_parts_at_either_end_are_harmless() {
        let agg = Aggregate::new([part(0, 0), part(100, 2), part(200, 0)]);

        assert_eq!(agg.len(), 2);
        assert_eq!(agg.get(0), Some(Num(100)));
        assert_eq!(agg.get(1), Some(Num(101)));
        assert_eq!(agg.get(2), None);
    }

    #[test]
    fn an_aggregate_of_nothing_is_an_empty_collection() {
        let agg: Aggregate<Part> = Aggregate::default();

        assert_eq!(agg.len(), 0);
        assert!(agg.is_empty());
        assert_eq!(agg.get(0), None);
        assert_eq!(agg.locate(0), None);
        assert_eq!(agg.iter().count(), 0);

        // and one holding only empty parts
        let empty = Aggregate::new([part(0, 0), part(1, 0)]);
        assert!(empty.is_empty());
        assert_eq!(empty.get(0), None);
    }

    #[test]
    fn pushing_extends_the_addressing_without_disturbing_it() {
        let mut agg = Aggregate::new([part(0, 2)]);
        assert_eq!(agg.len(), 2);

        agg.push(part(100, 3));
        assert_eq!(agg.len(), 5);
        assert_eq!(agg.get(1), Some(Num(1)), "the first part still addresses");
        assert_eq!(agg.get(2), Some(Num(100)));
        assert_eq!(agg.get(4), Some(Num(102)));
    }

    #[test]
    fn walking_an_aggregate_yields_every_part_in_order() {
        let agg = Aggregate::new([part(0, 2), part(100, 3)]);

        assert_eq!(
            agg.iter().collect::<Vec<_>>(),
            [Num(0), Num(1), Num(100), Num(101), Num(102)]
        );
    }

    #[test]
    fn the_combinators_work_over_an_aggregate() {
        // Aggregate is Indexable, so it composes like any
        // other collection rather than being a special case
        let agg = Aggregate::new([part(0, 3), part(100, 3)]);

        let high = agg.filter(|r| r.0 >= 100);
        assert_eq!(high.len(), 3);
        assert_eq!(high.get(0), Some(Num(100)));

        let sorted = Aggregate::new([part(100, 2), part(0, 2)]).sort_by_key(|r| r.0);
        assert_eq!(
            sorted.iter().collect::<Vec<_>>(),
            [Num(0), Num(1), Num(100), Num(101)]
        );

        assert_eq!(
            Aggregate::new([part(0, 3), part(100, 3)])
                .sample(4, 0)
                .len(),
            4
        );
    }

    #[test]
    fn an_aggregate_nests() {
        let inner = Aggregate::new([part(0, 2), part(100, 2)]);
        let outer = Aggregate::new([inner, Aggregate::new([part(200, 2)])]);

        assert_eq!(outer.len(), 6);
        assert_eq!(outer.get(3), Some(Num(101)));
        assert_eq!(outer.get(4), Some(Num(200)));
        assert_eq!(
            outer.locate(4),
            Some((1, 0)),
            "the outer part, then its own index"
        );
    }

    #[test]
    fn an_aggregate_collects_from_an_iterator() {
        let agg: Aggregate<Part> = (0..3).map(|i| part(i * 100, 2)).collect();

        assert_eq!(agg.len(), 6);
        assert_eq!(agg.get(2), Some(Num(100)));
    }

    #[test]
    fn an_aggregate_of_lenders_still_lends() {
        // Yield is forwarded, so aggregating collections that
        // hold their records does not force a clone back in
        fn lends<'a, C: Indexable<Yield<'a> = &'a Unit>>(_: &'a C) {}

        let agg = Aggregate::new([Lender(vec![Unit]), Lender(vec![Unit, Unit])]);
        assert_eq!(agg.len(), 3);
        lends(&agg);

        assert_eq!(agg.get(2), Some(&Unit));
    }

    #[test]
    fn a_part_that_breaks_the_contract_is_still_loud() {
        // Gappy counts three and refuses the middle one; the
        // aggregate forwards the None, so the panic names the
        // aggregate's index. locate() is how a caller turns
        // that back into a part
        let agg = Aggregate::new([Gappy, Gappy]);
        assert_eq!(agg.len(), 6);

        let err = std::panic::catch_unwind(|| agg.iter().count()).unwrap_err();
        let msg = err.downcast_ref::<String>().expect("a formatted message");
        assert!(msg.contains("get(1)"), "{msg}");
        assert!(msg.contains("len() is 6"), "{msg}");
    }

    // ---

    fn nums(from: u32, n: usize) -> Part {
        part(from, n)
    }

    fn seen<C>(c: &C) -> Vec<u32>
    where
        C: Indexable<Record = Num>,
    {
        c.iter().map(|r| r.borrow().0).collect()
    }

    #[test]
    fn a_shared_reference_is_a_collection_so_one_can_be_adapted_twice() {
        // without this every combinator consumes the
        // collection, and using it a second time means
        // rebuilding it -- or reopening the file
        let all = nums(0, 6);

        let evens = (&all).filter(|r| r.0 % 2 == 0);
        let odds = (&all).filter(|r| r.0 % 2 == 1);

        assert_eq!(seen(&evens), [0, 2, 4]);
        assert_eq!(seen(&odds), [1, 3, 5]);
        assert_eq!(all.len(), 6, "and the collection is still there");
    }

    #[test]
    fn a_borrowed_collection_reports_what_it_borrows() {
        let all = nums(10, 3);
        let borrowed = &all;

        assert_eq!(borrowed.len(), 3);
        assert_eq!(borrowed.get(1), Some(Num(11)));
        assert_eq!(borrowed.get(3), None);
        assert_eq!(seen(&borrowed), [10, 11, 12]);
    }

    #[test]
    fn take_and_skip_cut_a_contiguous_run() {
        assert_eq!(seen(&nums(0, 5).take(2)), [0, 1]);
        assert_eq!(seen(&nums(0, 5).skip(3)), [3, 4]);
        assert_eq!(seen(&(&nums(0, 5)).skip(1).take(2)), [1, 2], "they compose");
    }

    #[test]
    fn take_and_skip_clamp_rather_than_panic() {
        // asking for more than there is is the ordinary case
        // when a caller is paging through a file
        assert_eq!(seen(&nums(0, 3).take(99)), [0, 1, 2]);
        assert_eq!(nums(0, 3).skip(99).len(), 0);
        assert_eq!(nums(0, 0).take(5).len(), 0);
    }

    #[test]
    fn a_sample_reproduces_from_its_seed_and_a_different_seed_draws_differently() {
        let drawn = |seed: u64| seen(&(&nums(0, 100)).sample(10, seed));

        assert_eq!(drawn(7), drawn(7), "the same seed draws the same records");
        assert_eq!(drawn(7).len(), 10);
        assert_ne!(drawn(7), drawn(8), "a different seed draws differently");
    }

    #[test]
    fn a_permutation_reproduces_from_its_seed_and_holds_every_record_once() {
        let shuffled = |seed: u64| seen(&(&nums(0, 20)).permute(seed));

        assert_eq!(shuffled(3), shuffled(3));

        let mut sorted = shuffled(3);
        sorted.sort();
        assert_eq!(sorted, (0..20).collect::<Vec<_>>(), "every record, once");
    }

    #[test]
    fn a_sample_is_a_prefix_of_a_larger_one_drawn_under_the_same_seed() {
        // what a computed permutation gives that a shuffled
        // vector of indices does not: growing the draw keeps
        // everything the smaller one held
        let all = nums(0, 100);
        let small = seen(&(&all).sample(10, 7));
        let large = seen(&(&all).sample(40, 7));

        assert_eq!(small, large[..10]);
    }

    #[test]
    fn a_sample_of_every_record_is_the_permutation_itself() {
        let all = nums(0, 20);

        assert_eq!(seen(&(&all).sample(20, 5)), seen(&(&all).permute(5)));
    }

    #[test]
    fn drawing_from_an_empty_collection_is_empty_rather_than_a_panic() {
        // a feistel network is built over a non-empty
        // domain, so an empty collection has to be caught
        // before one is asked for
        assert_eq!(nums(0, 0).permute(1).len(), 0);
        assert_eq!(seen(&nums(0, 0).sample(0, 1)), []);
    }

    #[test]
    fn a_position_past_the_end_of_a_draw_is_none_rather_than_a_panic() {
        // Permutation::get asserts rather than answering,
        // so Drawn has to bound the call itself
        assert!(nums(0, 20).sample(5, 1).get(5).is_none());
        assert!(nums(0, 20).permute(1).get(20).is_none());
    }

    #[test]
    fn drawing_more_records_than_there_are_panics() {
        assert!(std::panic::catch_unwind(|| nums(0, 3).sample(4, 1)).is_err());
    }

    #[test]
    fn partition_splits_on_a_predicate_in_one_pass() {
        let all = nums(0, 6);
        let (small, large) = all.partition(|r| r.0 < 3);

        assert_eq!(seen(&small), [0, 1, 2]);
        assert_eq!(seen(&large), [3, 4, 5]);
        assert_eq!(small.len() + large.len(), all.len(), "nothing is lost");
    }

    #[test]
    fn chunks_are_consecutive_and_the_last_may_be_short() {
        let all = nums(0, 7);
        let chunks = all.chunks(3);

        assert_eq!(chunks.len(), 3);
        assert_eq!(seen(&chunks[0]), [0, 1, 2]);
        assert_eq!(seen(&chunks[2]), [6], "the remainder");
        assert_eq!(chunks.iter().map(|c| c.len()).sum::<usize>(), 7);
    }

    #[test]
    fn a_chunk_size_of_zero_does_not_hang() {
        let all = nums(0, 3);
        assert_eq!(all.chunks(0).len(), 3, "clamped to one record each");
        assert_eq!(nums(0, 0).chunks(4).len(), 0, "nothing to chunk");
    }

    #[test]
    fn split_always_yields_the_number_of_parts_asked_for() {
        // a caller sizing a thread pool by this would
        // otherwise get fewer workers than it asked for
        let all = nums(0, 7);
        let parts = all.split(3);

        assert_eq!(parts.len(), 3);
        assert_eq!(
            seen(&parts[0]),
            [0, 1, 2],
            "the remainder goes to the front"
        );
        assert_eq!(seen(&parts[1]), [3, 4]);
        assert_eq!(seen(&parts[2]), [5, 6]);

        let thin = nums(0, 2);
        let parts = thin.split(5);
        assert_eq!(parts.len(), 5, "even where two of them are empty");
        assert_eq!(parts.iter().map(|p| p.len()).sum::<usize>(), 2);
        assert_eq!(
            nums(0, 0).split(0).len(),
            1,
            "a split of zero is a split of one"
        );
    }

    #[test]
    fn split_weighted_balances_weight_where_split_balances_count() {
        // four cheap records and one expensive one, so an
        // even count is an uneven amount of work
        let all = nums(0, 5);
        let cost = |r: &Num| if r.0 == 4 { 10 } else { 1 };

        let parts = all.split_weighted(2, cost);
        assert_eq!(seen(&parts[0]), [4]);
        assert_eq!(
            seen(&parts[1]),
            [0, 1, 2, 3],
            "and a part holds its records in file order"
        );

        let by_count = all.split(2);
        assert_eq!(seen(&by_count[0]), [0, 1, 2]);
        assert_eq!(
            seen(&by_count[1]),
            [3, 4],
            "where split leaves the expensive record beside a cheap one"
        );
    }

    #[test]
    fn split_weighted_always_yields_the_number_of_parts_asked_for() {
        let thin = nums(0, 2);
        let parts = thin.split_weighted(5, |_| 1);

        assert_eq!(parts.len(), 5, "even where three of them are empty");
        assert_eq!(parts.iter().map(|p| p.len()).sum::<usize>(), 2);
        assert_eq!(
            nums(0, 0).split_weighted(0, |_| 1).len(),
            1,
            "a split of zero is a split of one"
        );
    }

    #[test]
    fn equal_weights_deal_the_same_way_every_time() {
        // a caller seeding a batch of jobs off the parts
        // gets the same batch from the same input
        let all = nums(0, 4);
        let once = all.split_weighted(2, |_| 1);
        let again = all.split_weighted(2, |_| 1);

        assert_eq!(seen(&once[0]), seen(&again[0]));
        assert_eq!(seen(&once[1]), seen(&again[1]));
        assert_eq!(seen(&once[0]), [0, 2]);
        assert_eq!(seen(&once[1]), [1, 3]);
    }

    #[test]
    fn group_by_keeps_file_order_inside_and_between_groups() {
        // grouping hits by query has to report queries in the
        // order the table lists them, not however they hash
        let all = nums(0, 6);
        let groups = all.group_by(|r| r.0 % 3);

        assert_eq!(groups.len(), 3);
        assert_eq!(groups[0].0, 0, "0 appears first");
        assert_eq!(seen(&groups[0].1), [0, 3]);
        assert_eq!(groups[1].0, 1);
        assert_eq!(seen(&groups[1].1), [1, 4]);
        assert!(all.group_by(|r| r.0).len() == 6, "all distinct keys");
        assert_eq!(nums(0, 0).group_by(|r| r.0).len(), 0);
    }

    #[test]
    fn find_and_position_stop_at_the_first_match() {
        let all = nums(10, 5);

        assert_eq!(all.position(|r| r.0 >= 12), Some(2));
        assert_eq!(all.find(|r| r.0 >= 12), Some(Num(12)));
        assert_eq!(all.position(|r| r.0 > 99), None);
        assert_eq!(all.find(|r| r.0 > 99), None);
    }

    #[test]
    fn min_and_max_by_key_pick_a_record_not_a_key() {
        let all = nums(0, 5);

        assert_eq!(all.min_by_key(|r| r.0), Some(Num(0)));
        assert_eq!(all.max_by_key(|r| r.0), Some(Num(4)));

        // the best hit is the smallest e-value, so reversing
        // the key is the ordinary use
        assert_eq!(all.min_by_key(|r| std::cmp::Reverse(r.0)), Some(Num(4)));
        assert_eq!(nums(0, 0).min_by_key(|r| r.0), None);
        assert_eq!(nums(0, 0).max_by_key(|r| r.0), None);
    }

    #[test]
    fn index_by_maps_a_key_to_a_position() {
        let all = nums(100, 4);
        let map = all.index_by(|r| r.0);

        assert_eq!(map.len(), 4);
        assert_eq!(map[&102], 2);
        assert_eq!(all.get(map[&103]), Some(Num(103)));
        assert!(!map.contains_key(&999));
    }

    #[test]
    fn a_repeated_key_keeps_the_first_record_that_had_it() {
        // duplicate identifiers are a malformed file, and
        // answering a lookup with the last one would hand back
        // a record the file does not lead with
        struct Same(usize);
        impl Indexable for Same {
            type Record = Num;
            type Yield<'a> = Num;
            fn len(&self) -> usize {
                self.0
            }
            fn get(&self, n: usize) -> Option<Num> {
                (n < self.0).then_some(Num(7))
            }
        }

        let map = Same(3).index_by(|r| r.0);
        assert_eq!(map.len(), 1);
        assert_eq!(map[&7], 0);
    }

    #[test]
    fn the_borrowing_adapters_work_over_an_aggregate_too() {
        // Aggregate is Indexable, so &Aggregate is as well,
        // and splitting a set of files across workers is the
        // pairing this exists for
        let agg = Aggregate::new([part(0, 3), part(100, 3)]);
        let parts = agg.split(2);

        assert_eq!(parts.len(), 2);
        assert_eq!(seen(&parts[0]), [0, 1, 2]);
        assert_eq!(seen(&parts[1]), [100, 101, 102]);

        let groups = agg.group_by(|r| r.0 >= 100);
        assert_eq!(groups.len(), 2);
    }
}

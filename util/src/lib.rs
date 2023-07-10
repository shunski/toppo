#![allow(unused)]

use std::collections::binary_heap::Iter;
use std::thread::current;

pub mod bitwise {
    pub fn count_ones(mut data: u64) ->  usize {
        let mut count = 0;
        while data != 0 {
            count+=1;
            data &= data-1;
        }
        return count;
    }
}



#[derive(Clone, Eq, PartialEq)]
pub struct BinaryPool {
    // "len" is the number of 0's and 1's
    len: usize,

    // "n_elements" is the number of 1's
    n_elements: usize,

    data: Vec<u64>
}

impl BinaryPool {
    pub fn new(len: usize) -> BinaryPool {
        BinaryPool {
            len: len,
            n_elements: 0,
            data: vec![0; len/64+1],
        }
    }

    pub fn from(len: usize, indeces: Vec<usize>) -> BinaryPool {
        if indeces.is_empty() { return BinaryPool::new(len) }
        if( len < indeces.iter().max().unwrap()/64+1) { panic!("specified length is too small."); }

        let data = vec![0; len/64+1];

        let mut pool = BinaryPool { 
            len: len, 
            n_elements: 0, 
            data: data,
        };

        indeces.iter().for_each(|&index| pool.add(index));
        pool
    }

    pub fn resize(&mut self, new_size: usize) {
        if new_size < self.len { return }

        while new_size/64 >= self.data.len() {
            self.data.push(0);
        }
        self.len = new_size;
    }

    pub fn add(&mut self, index: usize) {
        self.check_bounds(index);

        if !self.contains(index) {
            self.n_elements += 1;
        }

        self.data[index/64] |= 1 << (index % 64);
    }

    pub fn remove(&mut self, index: usize) {
        self.check_bounds(index);

        if self.contains(index) {
            self.n_elements -= 1;
        }

        self.data[index/64] &= !(1 << (index % 64));
    }

    pub fn union(mut self, mut other: Self ) -> Self {
        self.check_len(&other);
        self.data.iter_mut().zip(other.data.iter()).for_each(|(x, &y)| *x |= y);
        self.n_elements = self.data.iter().map(|&s| bitwise::count_ones(s)).sum();
        self
    }

    pub fn intersection(mut self, mut other: Self ) -> Self {
        self.check_len(&other);
        self.data.iter_mut().zip( other.data.iter() ).for_each(|(x, &y)| *x &= y);
        self.n_elements = self.data.iter().map(|&s| bitwise::count_ones(s)).sum();
        self
    }

    pub fn fill(&mut self) {
        (0..self.len()).for_each( |i| self.add(i) );
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline]
    pub fn n_elements(&self) -> usize {
        self.n_elements
    }

    #[inline]
    pub fn contains(&self, index: usize) -> bool {
        self.check_bounds(index);
        self.data[index/64] & (1 << (index%64)) == 1 << (index%64)
    } 

    #[inline]
    fn check_bounds(&self, index: usize) {
        if index >= self.len() { panic!("index out of bounds: len={} but index={}", self.len(), index); };
    }

    pub fn is_subset_of(&self, other: &Self) -> bool {
        self.check_len(other);
        self.data.iter().zip( other.data.iter() ).all(|(&x, &y)| x & y == x)
    }

    pub fn check_len(&self, other: &Self) {
        if self.len != other.len { panic!("two pools are not comparable. One has len {} but the other has len {}", self.len, other.len ) }
    }

    pub fn is_full(&self) -> bool {
        self.n_elements == self.len
    }
}


pub struct BinaryPoolIter<'a> {
    data: &'a BinaryPool,
    curr: std::ops::Range<usize>,
}

impl<'a> BinaryPool {
    pub fn iter(&'a self) -> BinaryPoolIter {
        BinaryPoolIter { 
            data: self, 
            curr: (0..self.len), 
        }
    }
}

impl<'a> Iterator for BinaryPoolIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        self.curr.find(|&i| self.data.contains(i))
    }
}

use std::fmt::{self, Binary};
impl fmt::Debug for BinaryPool {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "len = {}, ", self.len())?;
        write!(f, "n_elements = {}, ", self.n_elements())?;
        write!(f, "data = [")?;
        for i in 0..self.len() {
            if self.contains(i) {
                write!(f, "1")?;
            } else {
                write!(f, "0")?;
            }
            if i!=self.len()-1 {
                write!(f, ", ")?;
            }
        }
        write!(f, "]")?;
        write!(f, "")
    }
}


#[cfg(test)]
mod binary_pool_test {
    use crate::BinaryPool;

    #[test]
    fn binary_test() {
        let x = BinaryPool::from(10, Vec::new());
        assert_eq!(x.len(), 10);

        let y = BinaryPool::from(10, (0..10).collect());
        assert_eq!(y.len(), 10);

        let mut z = BinaryPool::from(10, vec![1, 3, 5, 7, 9]);
        let mut w = BinaryPool::from(10, vec![0, 2, 4, 6, 8]);

        let a = z.clone().intersection(w.clone());
        let b = z.clone().union(w.clone());

        assert_eq!(x.n_elements(), 0);
        assert_eq!(y.n_elements(), 10);
        assert_eq!(z.n_elements(), 5);
        assert_eq!(w.n_elements(), 5);

        assert_eq!(a, x);
        assert_eq!(b, y);

        assert_eq!(x.is_subset_of(&x), true);
        assert_eq!(x.is_subset_of(&y), true);
        assert_eq!(y.is_subset_of(&x), false);
        assert_eq!(z.is_subset_of(&w), false);
        assert_eq!(w.is_subset_of(&z), false);
        assert_eq!(x.is_subset_of(&w), true);
        assert_eq!(x.is_subset_of(&z), true);
        assert_eq!(z.is_subset_of(&y), true);
        assert_eq!(w.is_subset_of(&y), true);

        z.add(2);
        w.add(1);

        let a = z.clone().intersection(w.clone());
        let b = z.clone().union(w.clone());

        assert_eq!(z.n_elements(), 6);
        assert_eq!(w.n_elements(), 6);
        assert_eq!(a.n_elements(), 2);
        assert_eq!(b.n_elements(), 10);

        assert_eq!(a, BinaryPool::from(10, vec![1,2]));

        assert_eq!(x.is_subset_of(&y), true);
        assert_eq!(y.is_subset_of(&x), false);
        assert_eq!(z.is_subset_of(&w), false);
        assert_eq!(w.is_subset_of(&z), false);
        assert_eq!(x.is_subset_of(&w), true);
        assert_eq!(x.is_subset_of(&z), true);
        assert_eq!(z.is_subset_of(&y), true);
        assert_eq!(w.is_subset_of(&y), true);

        z.resize(100);
        let x = BinaryPool::from(100, (0..100).collect());
        assert_eq!(z.is_subset_of(&x), true);
        assert_eq!(x.is_subset_of(&z), false);
    }

    #[test] 
    fn iter_test() {
        let v = vec![1, 4, 5, 7, 8, 9];
        let x = BinaryPool::from(10, v.clone());
        assert_eq!( x.iter().collect::<Vec<_>>(), v);
    }
}



#[derive(Clone, PartialEq, Eq, Debug)]
pub struct LabeledSet<T> 
    where T: PartialEq + Eq + std::fmt::Debug
{
    set: BinaryPool,
    labels: Vec<T>,
}

impl <T> LabeledSet<T>
    where T: PartialEq + Eq + Clone + std::fmt::Debug
{
    pub fn new(labels: Vec<T>) -> LabeledSet<T> {
        LabeledSet { 
            set: BinaryPool::new( labels.len() ), 
            labels: labels,
        }
    }

    pub fn add(&mut self, label: T) {
        match self.labels.iter().position(|l| l == &label) {
            Some(p) => self.set.add(p),
            None => panic!("The label '{:?}' is not registered", label),
        }
    }

    pub fn remove(&mut self, label: T) {
        match self.labels.iter().position(|l| l == &label) {
            Some(p) => self.set.remove(p),
            None => panic!("The label '{:?}' is not registered", label),
        }
    }

    pub fn n_elements(&self) -> usize {
        self.set.n_elements()
    }
}

pub struct LabeledSetIter<'a, T> 
    where T: PartialEq + Eq + Clone + std::fmt::Debug
{
    data: &'a LabeledSet<T>,
    curr: BinaryPoolIter<'a>
}

impl<'a, T> Iterator for LabeledSetIter<'a, T>
    where T: PartialEq + Eq + Clone + std::fmt::Debug
{
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(i) = self.curr.next() {
            Some(&self.data.labels[i])
        } else {
            None
        }
    }
}

impl<'a, T> LabeledSet<T> 
    where T: PartialEq + Eq + Clone + std::fmt::Debug
{
    pub fn iter(&'a self) -> LabeledSetIter<'a, T> {
        LabeledSetIter { 
            data: self, 
            curr: self.set.iter()
        }
    }
}


impl<T> fmt::Display for LabeledSet<T> 
    where T: PartialEq + Eq + Clone + std::fmt::Debug + std::fmt::Display
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{{")?;
        let mut iter = self.iter();
        if let Some(elem) = iter.next() {
            write!(f, "{}", elem);
            while let Some(elem) = iter.next() { write!(f, ", {}", elem); }
        }
        write!(f, "}}")?;
        write!(f, "")
    }
}


#[cfg(test)]
mod labled_set_test {
    use crate::LabeledSet;

    #[test]
    fn init_test() {
        LabeledSet::new( vec!["1", "2", "3"] );
    }

    #[test]
    fn deref_test() {
        let mut set1 = LabeledSet::new( vec!["1", "2", "3"] );
        let mut set2 = set1.clone();

        set1.add("1"); set1.add("2"); set1.add("3"); set1.remove("1"); set1.remove("1"); 
        set2.add("2"); set2.add("3");

        assert_eq!(set1, set2);
        assert_eq!(set1.n_elements(), 2);
    }

    #[test]
    fn iter_test() {
        let mut set = LabeledSet::new( vec!["1", "2", "3"] );
        assert_eq!( set.iter().map(|&x| x).collect::<Vec<_>>(), Vec::<&str>::new());
        set.add("1"); 
        assert_eq!( set.iter().map(|&x| x).collect::<Vec<_>>(), vec!["1"]);
        set.add("2"); set.add("3"); set.remove("1"); 
        assert_eq!( set.iter().map(|&x| x).collect::<Vec<_>>(), vec!["2", "3"]);
    }
}


pub struct TupleIterator <'a, T> {
    data: &'a[T],
    indeces: Vec<usize>,
    curr_idx: usize,
    is_done: bool,
}

impl<'a, T> TupleIterator<'a, T> {
    fn new(data: &'a[T], size: usize) -> TupleIterator<'a, T> {
        let indeces = (0..size).collect();

        let is_done = data.len() < size;

        TupleIterator{
            data: data,
            indeces: indeces,
            curr_idx: size-1,
            is_done: is_done,
        }
    }

    fn is_blocked(&self) -> bool {
        if self.curr_idx >= self.indeces.len()-1 {
            if self.indeces[self.curr_idx] >= self.data.len()-1 {
                true
            } else {
                false
            }
        } else {
            if self.indeces[self.curr_idx]+1 == self.indeces[self.curr_idx+1] {
                true
            } else {
                false
            }
        }
    }
}

impl<'a, T: std::fmt::Debug> Iterator for TupleIterator<'a, T> {
    type Item = Vec<&'a T>;

    fn next(&mut self) -> Option<Self::Item> {
        let out = if self.is_done {
            return None;
        } else {
            Some(self.indeces.iter().map(|&idx| &self.data[idx]).collect())
        };

        if self.is_blocked() {
            while self.curr_idx > 0 && self.indeces[self.curr_idx]-self.indeces[self.curr_idx-1]==1 {
                self.curr_idx-=1;
            }

            if self.curr_idx == 0 {
                self.is_done = true;
            } else {
                self.indeces[self.curr_idx-1] += 1;
                let diff = self.indeces[self.curr_idx] - self.indeces[self.curr_idx-1] - 1;
                for idx in self.curr_idx..self.indeces.len() {
                    self.indeces[idx] -= diff;
                }

                self.curr_idx = self.indeces.len()-1;
            }
        }

        else{self.indeces[self.curr_idx] += 1;}
        out
    }
}

pub trait TupleIterable <'a, T> {
    fn tuple_iter(&'a self, size: usize) -> TupleIterator<'a, T>;
}

impl<'a, T> TupleIterable<'a, T> for Vec<T> {
    fn tuple_iter(&'a self, size: usize) -> TupleIterator<'a, T> {
        TupleIterator::new(self, size)
    }
}


#[cfg(test)]
mod tuple_iterator_test {
    use crate::TupleIterable;

    #[test]
    fn init() {
        let v:Vec<usize> = (0..5).collect();
        let mut tuple_iter = v.tuple_iter(3);
        assert_eq!(tuple_iter.next(), Some(vec![&0, &1, &2]));
        assert_eq!(tuple_iter.next(), Some(vec![&0, &1, &3]));
        assert_eq!(tuple_iter.next(), Some(vec![&0, &1, &4]));
        assert_eq!(tuple_iter.next(), Some(vec![&0, &2, &3]));
        assert_eq!(tuple_iter.next(), Some(vec![&0, &2, &4]));
        assert_eq!(tuple_iter.next(), Some(vec![&0, &3, &4]));
        assert_eq!(tuple_iter.next(), Some(vec![&1, &2, &3]));
        assert_eq!(tuple_iter.next(), Some(vec![&1, &2, &4]));
        assert_eq!(tuple_iter.next(), Some(vec![&1, &3, &4]));
        assert_eq!(tuple_iter.next(), Some(vec![&2, &3, &4]));
        assert_eq!(tuple_iter.next(), None);


        let v: Vec<usize> = Vec::new();
        for i in 1..10 {
            let mut tuple_iter = v.tuple_iter(i);
            assert_eq!(tuple_iter.next(), None);
        }

        let v:Vec<usize> = vec![0; 20];
        let sum: usize = v.tuple_iter(10).map(|_| 1).sum();
        assert_eq!(sum, choose(20, 10));
    }

    fn factorial(n: usize) -> usize {
        if n==0 { return 1; };
        n * factorial(n-1)
    }

    fn choose(n: usize, k: usize) -> usize {
        let mut out=1;
        for i in (k+1..=n) {
            out *= i;
            out /= i-k;
        }
        out
    }
}
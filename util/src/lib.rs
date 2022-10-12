#![allow(unused)]

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


#[derive(Clone, PartialEq)]
pub struct BinaryPool {
    len: usize,
    data: Vec<u64>
}

impl BinaryPool {
    pub fn new(len: usize) -> BinaryPool {
        BinaryPool {
            len: len,
            data: vec![0; len/64+1],
        }
    }

    pub fn from(len: usize, indeces: Vec<usize>) -> BinaryPool {
        if indeces.is_empty() { return BinaryPool::new(len) }
        if( len < indeces.iter().max().unwrap()/64+1) { panic!("specified length is too small."); }

        let data = vec![0; len/64+1];

        let mut pool = BinaryPool { 
            len: len, 
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
        self.data[index/64] |= 1 << (index % 64);
    }

    pub fn subtract(&mut self, index: usize) {
        self.check_bounds(index);
        self.data[index/64] &= !(1 << (index % 64));
    }

    pub fn union(mut self, mut other: Self ) -> Self {
        self.check_len(&other);
        self.data.iter_mut().zip(other.data.iter()).for_each(|(x, &y)| *x |= y);
        self
    }

    pub fn intersection(mut self, mut other: Self ) -> Self {
        self.check_len(&other);
        self.data.iter_mut().zip( other.data.iter() ).for_each(|(x, &y)| *x &= y);
        self
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn contains(&self, index: usize) -> bool {
        self.check_bounds(index);
        self.data[index/64] & (1 << (index%64)) == 1 << (index%64)
    } 

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
}

use std::fmt;
impl fmt::Debug for BinaryPool {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "len = {}, ", self.len())?;
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
        assert_eq!(z.clone().intersection(w.clone()), x);
        assert_eq!(z.clone().union(w.clone()), y);

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

        assert_eq!(z.clone().intersection(w.clone()), BinaryPool::from(10, vec![1,2]));

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
}


pub struct TupleIterator <'a, T: Copy> {
    data: &'a[T],
    indeces: Vec<usize>,
    curr_idx: usize,
    is_done: bool,
}

impl<'a, T: Copy> TupleIterator<'a, T> {
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

impl<'a, T: Copy + std::fmt::Debug> Iterator for TupleIterator<'a, T> {
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Self::Item> {
        let out;
        if self.is_done {
            return None;
        } else {
            out = Some(self.indeces.iter().map(|&idx| self.data[idx]).collect());
        }

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

pub trait TupleIterable <'a, T: Copy> {
    fn tuple_iter(&'a self, size: usize) -> TupleIterator<'a, T>;
}

impl<'a, T: Copy> TupleIterable<'a, T> for Vec<T> {
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
        assert_eq!(tuple_iter.next(), Some(vec![0,1,2]));
        assert_eq!(tuple_iter.next(), Some(vec![0,1,3]));
        assert_eq!(tuple_iter.next(), Some(vec![0,1,4]));
        assert_eq!(tuple_iter.next(), Some(vec![0,2,3]));
        assert_eq!(tuple_iter.next(), Some(vec![0,2,4]));
        assert_eq!(tuple_iter.next(), Some(vec![0,3,4]));
        assert_eq!(tuple_iter.next(), Some(vec![1,2,3]));
        assert_eq!(tuple_iter.next(), Some(vec![1,2,4]));
        assert_eq!(tuple_iter.next(), Some(vec![1,3,4]));
        assert_eq!(tuple_iter.next(), Some(vec![2,3,4]));
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
use std::iter::Product;
use std::ops::Mul;
use std::ops::BitXor;
use std::cmp::Eq;
use std::fmt::Display;
use crate::non_commutative::Group;
use crate::non_commutative::Gset;


#[derive(PartialEq, Eq, Debug, Clone)]
pub struct Cycle {
    data: Vec<usize>,
}

#[allow(unused)]
impl Cycle {
    fn is_valid(&self) -> bool {
        let mut checker = Vec::new();
        for &elem in self.iter() {
            if checker.contains(&elem) {
                return false;
            } else {
                checker.push(elem);
            }
        }
        true
    }

    fn cyclic_from(r: std::ops::Range<usize>) ->  Cycle {
        let data: Vec<_> = r.collect();
        if data.len() == 1 { return Cycle::identity(); };

        Cycle {
            data: data,
        }
        // this constructor does not need 'bring_smallest_frist()' since the smallest number is already the first
    }

    fn identity() -> Cycle {
        Cycle{ data: Vec::new() }
        // this constructor does not need 'bring_smallest_frist()' since this cycle does not have any element
    }

    pub fn from(data: Vec<usize>) -> Cycle {
        if data.len() == 1 { return Cycle::identity(); }

        let mut cycle = Cycle {
            data: data,
        };
        if !cycle.is_valid() { panic!("Input cycle is invalid.") }
        cycle = cycle.bring_smallest_first();
        cycle
    }

    fn bring_smallest_first(mut self) -> Self {
        // THIS FUNCTION MUST BE CALLED IN EVERY CONSTRUCTOR unless otherwise stated to be unnecessary.
        let (min_idx, &min) = self.iter().enumerate().min_by_key(|&(i, &e)| e).unwrap();

        // if the smallest index is already in the front, simply return.
        if min_idx == 0 { return self };

        // else, we bring the smallest element in the front 
        self.data = self.iter()
            .skip(min_idx)
            .map(|&elem| elem)
            .chain( self.iter().map(|&elem| elem).take_while(|&elem| elem!=min) )
            .collect::<Vec<_>>();
        self
    }

    fn map(&self, n: usize) -> usize {
        let mut pos  = match self.iter().position(|&k| k == n) {
            Some(pos) => pos,
            None => return n,
        };
        if pos==self.len()-1 {
            pos = 0;
        } else {
            pos += 1;
        }
        self[pos]
    }
}

impl std::ops::Deref for Cycle {
    type Target = Vec<usize>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

use std::cmp::Ordering;

impl Ord for Cycle {
    fn cmp(&self, other: &Self) -> Ordering {
        // we assume that the two input cycles are disjoint.

        // two cycles are equal
        if *self == *other {
            return Ordering::Equal;
        };

        // if one them is an empty cycle, the empty one is always smaller
        if self.is_empty() {
            return Ordering::Less;
        } else if other.is_empty() {
            return Ordering::Greater;
        }

        // Since the cycles are disjoint and nonempty, now the order is copmletely determined by their first element
        self.first().unwrap().cmp(other.first().unwrap())
    }
}

impl PartialOrd for Cycle {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some( self.cmp(other) )
    }
}

impl Display for Cycle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_empty() {
            return write!(f, "(1)");
        }

        write!(f, "(")?;
        for num in 0..self.len()-1 {
            write!(f, "{}, ", self[num])?;
        }
        // printing the last element
        write!(f, "{})", self.last().unwrap())
    }
}


#[cfg(test)]
mod cycle_test {
    use crate::non_commutative::permutation::Cycle;

    #[test]
    fn init_test() {
        // test 1
        let c1 = Cycle::from(vec![0,1,2]);
        let c2 = Cycle::from(vec![2,0,1]);
        let c3 = Cycle::from(vec![1,2,0]);
        let c4 = Cycle::cyclic_from(0..3);
        assert_eq!(c1, c2);
        assert_eq!(c2, c3);
        assert_eq!(c3, c4);

        // test 2
        let c1 = Cycle::identity();
        let c2 = Cycle::from(vec![1]);
        let c3 = Cycle::from(vec![2]);
        assert_eq!(c1, c2);
        assert_eq!(c2, c3);
    }

    #[test]
    #[should_panic]
    fn construction_panic() {
        Cycle::from(vec![0,1,2,2]);
    }
}

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct Permutation {
    // The variable 'cycles' is a collection of 'Cycle's that satisfies the following properties:
    // (1) if it has a non-empty 'Cycle,' then it contains no empty cycle.
    // (2) if it contains no non-empty cycle, then it contains exactly one empty cycle.
    // (3) 'Cycle's in 'cycles' must be disjoint
    // These condition must be satisfied at the construction of the object.
    cycles: Vec<Cycle>,
}


#[allow(unused)]
impl Permutation {
    fn is_valid(&self) -> bool {
        let mut checker = Vec::new();
        for elem in self.unfixed_indeces_iter() {
                if checker.contains(&elem) {
                    return false;
                } else {
                    checker.push(elem);
                }
        }
        true
    }

    pub fn cyclic_from(r: std::ops::Range<usize>) -> Self {
        Permutation{ cycles: vec![Cycle::cyclic_from(r)] }
    }

    pub fn swap(r1: std::ops::Range<usize>, r2: std::ops::Range<usize>) -> Self {
        // trivial cases:
        if r1.is_empty() || r2.is_empty() { Self::identity(); }

        // Two input ranges must be disjoint
        if !(r1.end <= r2.start || r2.end <= r1.start) { panic!("Two input ranges must be disjoint, but r1={:?} and r2={:?}", r1, r2); }

        // Now r1.end <= r2.start || r2.end <= r1.start, i.e. two ranges are disjoint
        if r2.end < r1.start { return Self::swap(r2, r1) };

        // Now r1.end <= r2.start
        let r1_size = r1.end - r1.start;
        let r2_size = r2.end - r2.start;

        return if r1_size < r2_size {
            Self::swap(r1.start..r2.start, r2.start..r2.end-r1_size) * Self::swap_ranges_of_the_same_length(r1, r2.end-r1_size..r2.end)
        } else if r1_size > r2_size {
            Self::swap(r1.start+r2_size..r1.end, r1.end..r2.end) * Self::swap_ranges_of_the_same_length(r1.start..r1.start+r2_size, r2)
        } else {
            Self::swap_ranges_of_the_same_length(r1, r2)
        }
    }

    fn swap_ranges_of_the_same_length(r1: std::ops::Range<usize>, r2: std::ops::Range<usize>) -> Self {
        let r1_size = r1.end - r1.start;
        let r2_size = r2.end - r2.start;

        debug_assert!( r1.end <= r2.start, "The condition 'r1.end <= r2.start' is required, but r1={:?} and r2={:?}", r1, r2 );
        debug_assert!( r1_size == r2_size , "The condition 'r1_size == r2_size' is required, but r1={:?} and r2={:?}", r1, r2 );

        let cycles = r1.zip( r2 ).map(|(i, j)| Cycle::from(vec![ i,j] )).collect();
        Permutation{ cycles: cycles }
    }

    pub fn from_disjoint_cycles(cycles: Vec<Cycle>) -> Permutation {
        // remove all the empty cycles
        let mut cycles: Vec<_> = cycles.into_iter().filter(|cycle| !cycle.is_empty()).collect();

        // if the resulting cycles are empty, then it is an identity permutation.
        if cycles.is_empty() { cycles.push( Cycle{ data: Vec::new() } ) }

        // order the cycles
        cycles.sort();

        // construct the permutation and check wheher it is valid
        let p = Permutation { cycles: cycles };
        if !p.is_valid() {
            panic!("Cannot construct permutation because the cycles are invalid. ");
        };

        // return the result
        p
    }

    pub fn from_cycle(cycle: Cycle) -> Permutation {
        return Permutation{ cycles: vec![ cycle ] };
    }

    pub fn map(&self, n: usize) -> usize {
        match self.iter().map(|cycle| cycle.map(n)).find(|&k| k != n) {
            Some(k) => k,
            None => n,
        }
    }
}


impl std::ops::Deref for Permutation {
    type Target = Vec<Cycle>;

    fn deref(&self) -> &Self::Target {
        &self.cycles
    }
}

pub struct UnfixedIndexIter<'a> {
    cycle_iter: std::slice::Iter<'a, Cycle>,
    unfixed_indeces_iter: std::slice::Iter<'a, usize>,
}

impl<'a> UnfixedIndexIter<'a> {
    fn new(p: &'a Permutation) -> Self {
        let mut cycle_iter = p.cycles.iter();
        let unfixed_indeces_iter = if let Some(x) = cycle_iter.next() {
            x.data.iter()
        } else {
            panic!("Permutation must contain at least one cycle.")
        };
        UnfixedIndexIter {
            cycle_iter: cycle_iter,
            unfixed_indeces_iter: unfixed_indeces_iter,
        }
    }
}


impl<'a> Iterator for UnfixedIndexIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(&elem) = self.unfixed_indeces_iter.next() {
            return Some(elem);
        };

        // else we go to the next cycle
        if let Some(cycle) = self.cycle_iter.next() {
            self.unfixed_indeces_iter = cycle.iter();
            match self.unfixed_indeces_iter.next() {
                Some(&elem) => Some(elem),
                None => None
            }
        } else {
            // if there is no more cycle, the iteration is done, so return None
            None
        }
    }
}

impl<'a> Permutation {
    pub fn unfixed_indeces_iter(&'a self) -> UnfixedIndexIter<'a> {
        UnfixedIndexIter::new(self)
    }
}


impl Display for Permutation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for cycle in self.cycles.iter() {
            write!(f, "{}", cycle)?;
        }
        write!(f, "")
    }
}


impl Mul for Permutation {
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        let mut product = Vec::new();
        let mut unfixed_indeces: Vec<_> = other.unfixed_indeces_iter().map(|x| x).collect();
        self.unfixed_indeces_iter().for_each(|x| unfixed_indeces.push(x) );
        while !unfixed_indeces.is_empty() {
            let first_elem = *unfixed_indeces.first().unwrap();

            let mut cycle = vec![first_elem];
            let mut elem = first_elem;
            while self.map( other.map(elem) ) != first_elem {
                elem = self.map( other.map(elem) );
                cycle.push(elem);
            }

            // update the 'unfixed_indeces' for the next loop:
            // that is, remove 'cycle' from 'unfixed_indeces'
            unfixed_indeces = unfixed_indeces.into_iter().filter( |x| !cycle.contains(x) ).collect();

            // add the cycle
            product.push( Cycle::from(cycle) );
        };
        Permutation::from_disjoint_cycles(product)
    }
}


#[macro_export]
macro_rules!  permutation {
    ($(($($x:expr), +)) +) => {
        {
            let mut cycles = Vec::new();
            $(
                let mut v = Vec::new();
                $(
                    v.push($x);
                )+
                cycles.push( Cycle::from(v) );
            )+
            Permutation::from_disjoint_cycles( cycles )
        }
    };
}

// implementation of 'Product' for 'Permutation' to be a Group
impl Product for Permutation {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::identity(), |accum, p| p * accum )
    }
}

// implementation of 'BitXor' for 'Permutation' to be a Group
impl BitXor<i64> for Permutation {
    type Output = Permutation;
    fn bitxor(self, power: i64) -> Self::Output {
        if power == 0 {
            Self::identity()
        } else if power > 0 {
            (1..=power).map(|_| self.clone()).product()
        } else { // if 'power < 0'
            let inv = self.inverse();
            (1..=power).map(|_| inv.clone()).product()
        }
    }
}

// we define that the set of generators for symmetric groups is the set of transpositions (even though it can be smaller)
impl Group for Permutation {
    fn identity() -> Self {
        Permutation { cycles: vec![ Cycle::identity() ] }
    }

    fn inverse_of_gen(self) -> Self {
        // this inplemeentation is in fact trivial
        self
    }

    fn in_terms_of_gens(self) -> Vec<Self> {
        // in the case that self is the identity element
        if self == Self::identity() {
            return vec![self];
        }

        let mut gens = Vec::new();
        for cycle in self.cycles.iter() {
            let first_elem = *cycle.first().unwrap();
            cycle.iter().filter(|&&elem| first_elem!=elem).for_each(|&elem| gens.push( permutation!{ (first_elem, elem) } ));
        };

        debug_assert!(
            self==gens.iter().map(|gen| gen.clone()).product(), 
            "self is : {}, but the product of gens is: {}", self, gens.iter().map(|gen| gen.clone()).product::<Permutation>() 
        );

        gens
    }
}



#[cfg(test)]
mod permutation_test {
    use crate::non_commutative::permutation::*;

    #[test]
    fn init_test() {
        assert_eq!( permutation!{ (0, 1, 2, 3, 5) }.to_string(), "(0, 1, 2, 3, 5)" );
        assert_eq!( permutation!{ (0, 1, 2, 3)(5, 6, 4) }.to_string(), "(0, 1, 2, 3)(4, 5, 6)" );
    }

    #[test]
    #[should_panic]
    fn validness_test() {
        permutation!{ (0, 1, 2, 3, 5)(4, 5) };
    }

    #[test]
    fn map_test() {
        let p = permutation!((0, 2, 3, 6));
        assert_eq!( p.map(0), 2);
        assert_eq!( p.map(1), 1);
        assert_eq!( p.map(2), 3);
        assert_eq!( p.map(3), 6);
        assert_eq!( p.map(4), 4);
        assert_eq!( p.map(5), 5);
        assert_eq!( p.map(6), 0);
    }

    #[test]
    fn swap_test() {
        assert_eq!( 
            Permutation::swap(0..7, 7..9), 
            permutation!{(0,2,4,6,8,1,3,5,7)}
        );

        assert_eq!(
            Permutation::swap(2..8, 8..10), 
            permutation!{ (2, 4, 6, 8)(3, 5, 7, 9) }
        )
    }

    #[test]
    fn unfixed_indeces_iter_test() {
        // test 1
        let p = permutation!((0, 2, 3, 6));
        let unfixed_indeces: Vec<_> = p.unfixed_indeces_iter().collect();
        assert_eq!(unfixed_indeces, vec![0,2,3,6]);

        // test 2
        let p = permutation!((0, 2, 3, 6)(1,4,8));
        let unfixed_indeces: Vec<_> = p.unfixed_indeces_iter().collect();
        assert_eq!(unfixed_indeces, vec![0,2,3,6,1,4,8]);
    }

    #[test]
    fn multuplication_test() {
        let p = permutation!{ (0, 5, 4, 6, 1, 2) };
        let q = permutation!{ (0, 1, 2, 3, 5) };
        let r = permutation!{ (0, 2, 3, 4, 6, 1) };
        assert_eq!( p*q, r);
    }

    #[test]
    fn identity_test() {
        assert_eq!(Permutation::identity(), permutation!{ (0) });
        assert_eq!(Permutation::identity(), permutation!{ (1) });
        assert_eq!(Permutation::identity(), permutation!{ (100) });

        let id = Permutation::identity();
        let p = permutation!{ (0, 5, 4, 6, 1, 2) };
        let q = permutation!{ (0, 1, 2, 3, 5) };
        let r = permutation!{ (0, 2, 3)(1, 6)(4, 8, 5) };
        let s = Permutation::cyclic_from(0..50);

        assert_eq!(p.clone() * id.clone(), p);
        assert_eq!(id.clone() * p.clone(), p);
        assert_eq!(q.clone() * id.clone(), q);
        assert_eq!(id.clone() * q.clone(), q);
        assert_eq!(r.clone() * id.clone(), r);
        assert_eq!(id.clone() * r.clone(), r);
        assert_eq!(s.clone() * id.clone(), s);
        assert_eq!(id.clone() * s.clone(), s);
        assert_eq!(id.clone() * id.clone(), id);
    }

    #[test]
    fn inverse_test() {
        let id = Permutation::identity();
        let p = permutation!{ (0, 5, 4, 6, 1, 2) };
        let q = permutation!{ (0, 1, 2, 3, 5) };
        let r = permutation!{ (0, 2, 3)(1, 6)(4, 8, 5) };
        let s = Permutation::cyclic_from(0..50);

        assert_eq!(id.clone().inverse(), id);

        assert_eq!(p.clone().inverse(), permutation!{ (0, 2, 1, 6, 4, 5) });
        assert_eq!(p.clone() * p.clone().inverse(), id);

        assert_eq!(q.clone().inverse(), permutation!{ (0, 5, 3, 2, 1) });
        assert_eq!(q.clone() * q.clone().inverse(), id);

        assert_eq!(r.clone().inverse(), permutation!{ (0, 3, 2)(1, 6)(4, 5, 8) });
        assert_eq!(r.clone() * r.clone().inverse(), id);
        
        assert_eq!(s.clone() * s.clone().inverse(), id);
    }

    #[test]
    fn power_test() {
        let p = permutation!{ (0, 5, 4, 6, 1, 2) };
        let q = permutation!{ (0, 1, 2, 3, 5) };
        let r = permutation!{ (0, 2, 3)(1, 6)(4, 8, 5) };

        assert_eq!( (p.clone())^2, permutation!{ (0,4,1)(5,6,2) } );
        assert_eq!( p^6, Permutation::identity() );
        assert_eq!( q^5, Permutation::identity() );
        assert_eq!( r^6, Permutation::identity() );
    }
}


// symmetric action
impl<T: Clone> Gset<Permutation> for Vec<T> {
    fn gen_action_by(&mut self, elem: Permutation) {
        // First, process the case that 'elem' is the identity
        if elem == Permutation::identity() { return; }

        // again, self is a generator, it must be a 2-cycle in our convension.
        let unfixed_indeces: Vec<_> = elem.unfixed_indeces_iter().collect();

        debug_assert!(unfixed_indeces.len() == 2, "gen_act must receive a 2-cycle, but it receives: {}", unfixed_indeces.len());

        let (i, j) = (unfixed_indeces[0], unfixed_indeces[1]);
        self.swap(i, j);
    }
}

impl<T: Clone> std::ops::Shl<Permutation> for Vec<T> {
    type Output = Vec<T>;

    fn shl(mut self, elem: Permutation) -> Self::Output {
        self.action_by(elem);
        self
    }
}

impl<T: Clone> std::ops::ShlAssign<Permutation> for Vec<T> {
    fn shl_assign(&mut self, elem: Permutation) {
        self.action_by(elem);
    }
}


#[cfg(test)]
mod symmetric_action_test {
    use crate::non_commutative::permutation::*;
    #[test]
    fn action_test() {
        let mut v = vec![0,1,2,3,4,5];
        let p = permutation!{ (0,1,2)(3,4) };
        v <<= p;
        assert_eq!( v, vec![2,0,1,4,3,5] );

        let p = Permutation::cyclic_from(0..3);
        v <<= p;
        assert_eq!( v, vec![1,2,0,4,3,5] );

        let p = Permutation::swap(0..2, 4..5);
        v <<= p;
        assert_eq!( v, vec![3,0,4,1,2,5] );

        let mut v:Vec<_> = (0..10).collect();

        let p = Permutation::swap(2..8, 8..10);
        v <<= p;
        assert_eq!( v, vec![0,1,8,9,2,3,4,5,6,7] );
    }
}
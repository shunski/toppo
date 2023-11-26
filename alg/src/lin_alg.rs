use crate::commutative::{PID, Field};
use std::cmp::Ordering;
use std::marker::PhantomData;
use std::ptr::NonNull;
use std::rc::Rc;

pub struct Matrix<T: PID> {
    ptr: NonNull<T>,
    _marker: PhantomData<T>,
    
    size: (usize, usize),
    offsets: (usize, usize),
}

pub struct SubMatrix<T: PID>([T]);

#[derive(Clone, PartialEq, Debug)]
pub struct SparseMatrix<T: PID> {
    rows: Vec<SparseNode<T>>,
    cols: Vec<SparseNode<T>>,
    size: (usize, usize),
    transposed: bool,
}

#[derive(PartialEq, Clone, Debug)]
enum SparseNode<T: PID> {
    Start(usize, Rc<SparseNode<T>>),
    End,
    Entry(
        /* index */ (usize, usize),
        /* value */ T,
        /* (next_row_node, next col node) */ (Rc<SparseNode<T>>, Rc<SparseNode<T>>)
    )
}

#[derive(Clone, PartialEq)]
pub struct ConstMatrix<T: PID, const N_ROW: usize, const N_COL: usize>([[T; N_COL]; N_ROW]);

#[derive(Clone, PartialEq, Debug)] // ToDo: remove the Debug and implement its own Debug
pub struct ConstVector<T: PID, const N: usize> ([T; N]);

pub mod matrix_impl;
mod matrix_test;

mod sparse_matrix_impl;

mod const_matrix_impl;

#[macro_export]
macro_rules! matrix {
    ( $coeff: ty; [$([$( $elem:expr ), *]), *]) => {
        {
            let mut v: Vec<Vec<$coeff>> = Vec::new();
            $(
                v.push(vec!($( $elem ), *));
            )*
            Matrix::<$coeff>::from_array(v)
        }
    };

    ( [$([$( $elem:expr ), *]), *]) => {
        {
            let mut v: Vec<Vec<_>> = Vec::new();
            $(
                v.push(vec!($( $elem ), *));
            )*
            Matrix::from_array(v)
        }
    };
}

pub struct Module<Coeff: PID> {
    free_rank: usize,
    elementary_divisors: Vec<(Coeff, usize)>
}


impl<Coeff: PID> Module<Coeff> {
    pub fn new(free_rank: usize, elementary_divisors: Vec<(Coeff, usize)>) -> Self {
        // TODO: check whether the elementary_divisors are indeed elementary divisors
        let elementary_divisors: Vec<_> = elementary_divisors.into_iter().filter(|(div, _)| *div != Coeff::one()).collect();
        Module {
            free_rank: free_rank,
            elementary_divisors: elementary_divisors,
        }
    }

    pub fn trivial() -> Self {
        Module {
            free_rank: 0,
            elementary_divisors: Vec::new(),
        }
    }

    pub fn is_isomorphic_to (self, rhs: Self) -> bool {
        self.free_rank == rhs.free_rank && self.elementary_divisors == rhs.elementary_divisors
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        return self.free_rank==0 && self.elementary_divisors.is_empty();
    }
}

use std::fmt::Display;
impl<Coeff: PID+Display> Display for Module<Coeff> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // in case it is the zero module
        if self.is_zero() {
            return write!(f, "0");
        }

        // general case
        // format free part first
        if self.free_rank == 1 {
            write!(f, "{}", Coeff::symbol())?;
            if !self.elementary_divisors.is_empty() { write!(f, " + ")? };
        } else if self.free_rank > 1 {
            write!(f, "{}^{}", Coeff::symbol(), self.free_rank)?;
            if !self.elementary_divisors.is_empty() { write!(f, " + ")? };
        }

        
        // then format the torsion part
        if let Some((div, mul)) = self.elementary_divisors.first() {
            if *mul==1 {
                write!(f, "{}_{}", Coeff::symbol(), div)?;
            } else {
                write!(f, "({}_{})^{}", Coeff::symbol(), div, mul)?;
            }
        } else {
            return write!(f, "")
        }
        for (divisor, multiplicity) in self.elementary_divisors.iter().skip(1) {
            if *multiplicity==1 {
                write!(f, " + {}_{}", Coeff::symbol(), divisor)?;
            } else {
                write!(f, " + ({}_{})^{}", Coeff::symbol(), divisor, multiplicity)?;
            }
        }
        write!(f, "")
    }
}


#[cfg(test)]
mod module_fmt_test {
    use crate::lin_alg::Module;
    #[test]
    fn fmt_test() {
        let m = Module::<i64>::new(0,Vec::new());
        assert_eq!(m.to_string(), "0");

        let m = Module::<i64>::new(1,Vec::new());
        assert_eq!(m.to_string(), "Z");

        let m = Module::<i64>::new(0,vec![(2,1)]);
        assert_eq!(m.to_string(), "Z_2");

        let m = Module::<i64>::new(1,vec![(2,2)]);
        assert_eq!(m.to_string(), "Z + (Z_2)^2");

        let m = Module::<i64>::new(2,vec![(2,2), (8,1), (24, 2)]);
        assert_eq!(m.to_string(), "Z^2 + (Z_2)^2 + Z_8 + (Z_24)^2");

        let m = Module::<i64>::new(2,vec![(1,3), (3,2), (24, 2)]);
        assert_eq!(m.to_string(), "Z^2 + (Z_3)^2 + (Z_24)^2");
    }
}


#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Vector<Coeff: PID, Basis: Clone + PartialEq + PartialOrd> (
    Vec<(Coeff, Basis)>,
);

impl<Basis: Clone + PartialEq + PartialOrd, Coeff: PID> Vector<Coeff, Basis> {
    pub fn zero() -> Self {
        Vector(Vec::new())
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn remove(&mut self, i: usize) -> (Coeff, Basis) {
        self.0.remove(i)
    }
}

impl<const N: usize, Basis: Clone + PartialEq + PartialOrd, Coeff: PID> std::convert::From<[Basis; N]> for Vector<Coeff, Basis> {
    fn from(mut summands: [Basis; N]) -> Self {
        summands.sort_by(|x,y| x.partial_cmp(y).unwrap());
        Vector(summands.into_iter().map(|b| (Coeff::one(), b)).collect::<Vec<_>>())
    }
}

impl<Basis: Clone + PartialEq + PartialOrd, Coeff: PID> std::convert::From<Vec<Basis>> for Vector<Coeff, Basis> {
    fn from(mut summands: Vec<Basis>) -> Self {
        summands.sort_by(|x,y| x.partial_cmp(y).unwrap());
        Vector(summands.into_iter().map(|b| (Coeff::one(), b)).collect::<Vec<_>>())
    }
}

impl<Basis: Clone + PartialEq + PartialOrd, Coeff: PID> std::convert::From<Basis> for Vector<Coeff, Basis> {
    fn from(summand: Basis) -> Self {
        Vector(vec![(Coeff::one(), summand)])
    }
}

impl<Basis: Clone + PartialEq + PartialOrd, Coeff: PID> std::iter::FromIterator<Basis> for Vector<Coeff, Basis> {
    fn from_iter<T: IntoIterator<Item = Basis>>(iter: T) -> Self {
        let summands: Vec<_> = iter.into_iter().collect();
        Self::from(summands)
    }
}

impl<'a, Basis: Clone + PartialEq + PartialOrd, Coeff: PID> Vector<Coeff, Basis> {
    pub fn iter(&'a self) -> impl Iterator<Item=&'a(Coeff, Basis)> {
        self.0.iter()
    }

    pub fn iter_mut(&'a mut self) -> impl Iterator<Item=&'a mut (Coeff, Basis)> {
        self.0.iter_mut()
    }
}

impl<Basis: Clone + PartialEq + PartialOrd, Coeff: PID> Vector<Coeff, Basis> {
    pub fn into_iter(self) -> impl Iterator<Item=(Coeff, Basis)> {
        self.0.into_iter()
    }
}

impl<Basis: Clone + PartialEq + PartialOrd, Coeff: PID> std::ops::Add for Vector<Coeff, Basis> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        [self, rhs].into_iter().sum()
    }
}

impl<Basis: Clone + PartialEq + PartialOrd, Coeff: PID> std::ops::Sub for Vector<Coeff, Basis> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<Basis: Clone + PartialEq + PartialOrd, Coeff: PID> std::ops::Neg for Vector<Coeff, Basis> {
    type Output = Self;
    fn neg(mut self) -> Self::Output {
        self.iter_mut().for_each(|(s, _)| *s *= -Coeff::one());
        self
    }
}

impl<Basis: Clone + PartialEq + PartialOrd, Coeff: PID> std::ops::Mul<Coeff> for Vector<Coeff, Basis> {
    type Output = Vector<Coeff, Basis>;
    fn mul(mut self, rhs: Coeff) -> Self::Output {
        if rhs == Coeff::zero() { return Vector::zero(); };
        self.iter_mut().for_each(|(s, _)| *s *= rhs.clone() );
        self
    }
}

impl<Basis: Clone + PartialEq + PartialOrd, Coeff: Field> std::ops::Div<Coeff> for Vector<Coeff, Basis> {
    type Output = Vector<Coeff, Basis>;
    fn div(mut self, rhs: Coeff) -> Self::Output {
        if rhs == Coeff::zero() { panic!("cannot devide a vector by 0!") };
        self.iter_mut().for_each(|(s, _)| *s /= rhs.clone() );
        self
    }
}

impl<Basis: Clone + PartialEq + PartialOrd, Coeff: PID> std::ops::Add<Basis> for Vector<Coeff, Basis> {
    type Output = Self;
    fn add(self, rhs: Basis) -> Self::Output {
        self + Vector::from( rhs )
    }
}

impl<Basis: Clone + PartialEq + PartialOrd, Coeff: PID> std::ops::Sub<Basis> for Vector<Coeff, Basis> {
    type Output = Self;
    fn sub(self, rhs: Basis) -> Self::Output {
        self - Vector::from( rhs )
    }
}

struct Merge <T, I, J, F>
    where 
        I: Iterator<Item = T>,
        J: Iterator<Item = T>,
        F: Fn(&T, &T)->Ordering
{
    curr: Option<T>,
    i: I,
    j: J,
    f: F
}

impl<T, I, J, F> Iterator for Merge<T, I, J, F>
    where 
        I: Iterator<Item = T>,
        J: Iterator<Item = T>,
        F: Fn(&T, &T)->Ordering
{
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        let mut self_curr = None;
        std::mem::swap(&mut self_curr, &mut self.curr);
        if let Some(mut curr) = self_curr {
            if let Some(mut next) = self.i.next() {
                if (self.f)(&curr, &next) == Ordering::Less {
                    std::mem::swap( &mut next, &mut curr );
                }
                self.curr = Some(curr);
                Some(next)
            } else if let Some(next) = self.j.next() { 
                self.curr = Some(next);
                Some(curr)
            } else {
                self.curr = None;
                Some(curr)
            }
        } else {
            None
        }
    }
}

fn merge_by<T, I: Iterator<Item=T>, J: Iterator<Item=T>, F: Fn(&T, &T)->Ordering>(mut i: I, mut j: J, f: F) -> impl Iterator<Item = T> {
    if let Some(next) = j.next() {
        Merge{
            curr: Some(next), 
            i, j, f
        }
    } else {
        Merge{
            curr: i.next(),
            i,
            j,
            f,
        }
    }
}

#[cfg(test)]
mod merge_by {
    use super::merge_by;

    #[test]
    fn merge_by_test() {
        let a = [1,3];
        let b = [2];
        assert_eq!(merge_by(a.into_iter(), b.into_iter(), |x, y| x.cmp(y)).collect::<Vec<_>>(), [1,2,3]);
        assert_eq!(merge_by(b.into_iter(), a.into_iter(), |x, y| x.cmp(y)).collect::<Vec<_>>(), [1,2,3]);
    }
}

impl<Basis: Clone + PartialEq + PartialOrd, Coeff: PID> std::iter::Sum for Vector<Coeff, Basis> {
    fn sum<I: std::iter::Iterator<Item=Self>>(iter: I) -> Self {
        let mut iters: Vec<_> = iter
            .map(|v| v.into_iter())
            .map(|mut i| (i.next(), i))
            .filter(|(val, _)| val!=&None )
            .map(|(val, i)| (val.unwrap(), i))
            .collect();
        iters.sort_by(|((_, basis1),_), ((_, basis2),_)| basis1.partial_cmp(basis2).unwrap() );
        
        let mut out = Vec::new();

        while !iters.is_empty() {
            let n = iters.iter().take_while(|((_,basis), _)| basis == &iters[0].0.1 ).count();
            let mut consumed_iter = iters.into_iter();
            let mut next_iters = Vec::new();

            let ((mut coeff, basis), mut iter) = consumed_iter.next().unwrap();
            if let Some(x) = iter.next() {
                next_iters.push((x, iter));
            }
            for _ in 0..n-1 {
                let ((c, _), mut iter) = consumed_iter.next().unwrap();
                if let Some(x) = iter.next() {
                    next_iters.push((x, iter));
                }
                coeff += c;
            }
            
            next_iters.sort_by(|((_,x),_),((_,y),_)| x.partial_cmp(y).unwrap());
            iters = merge_by(consumed_iter, next_iters.into_iter(), |((_,x),_),((_,y),_)| x.partial_cmp(y).unwrap() ).collect();

            if coeff != Coeff::zero() { 
                out.push((coeff, basis))
            }
        }
        Vector(out) 
    }
}


// index trait for retrieving coefficients
impl<Basis: Clone+PartialEq+Ord, Coeff: PID> std::ops::Index<&Basis> for Vector<Coeff, Basis> {
    type Output = Coeff;
    fn index(&self, index: &Basis) -> &Self::Output {
        &self.iter().find(|(_, x)| x==index ).map(|(i,_)| i).unwrap()
    }
}


impl<Basis: Clone+PartialEq+PartialOrd+Display, Coeff: PID+Display> std::fmt::Display for Vector<Coeff, Basis> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // if the members are empty, then print zero and return
        if self == &Vector::zero() {
            return write!(f, "0");
        }

        // else we print elements
        let mut iter = self.iter();

        // printing the first element
        let (i, x) = iter.next().unwrap();
        if i==&Coeff::one() {
            write!(f, "{}", x)?;
        } else if i==&(-Coeff::one()) {
            write!(f, "-{}", x)?;
        } else {
            write!(f, "{}·{}", i, x)?;
        };

        // printing rest of the elements
        for (i, x) in iter {
                write!(f, " + ({})·{}", i, x)?;
        };
        write!(f, "")
    }
}


#[macro_export]
macro_rules! formal_sum_impl {
    ($($t:ty)*) => ($(
        impl std::ops::Add for $t {
            type Output = Vector<i64, $t>;
            fn add(self, rhs: Self) -> Self::Output { 
                Vector::from(self) + Vector::from(rhs)
            }
        }

        impl std::ops::Add<Vector<i64, $t>> for $t {
            type Output = Vector<i64, $t>;
            fn add(self, rhs: Vector<i64, $t>) -> Self::Output { 
                Vector::from(self) + rhs
            }
        }

        impl std::ops::Sub for $t {
            type Output = Vector<i64, $t>;
            fn sub(self, rhs: Self) -> Self::Output { 
                Vector::from(self) - Vector::from(rhs)
            }
        }

        impl std::ops::Sub<Vector<i64, $t>> for $t {
            type Output = Vector<i64, $t>;
            fn sub(self, rhs: Vector<i64, $t>) -> Self::Output { 
                Vector::from(self) - rhs
            }
        }

        impl std::ops::Mul<$t> for i64 {
            type Output = Vector<i64, $t>;
            fn mul(self, rhs: $t) -> Vector<i64, $t> { 
                Vector::from(rhs) * self
            }
        }

        impl std::ops::Neg for $t {
            type Output = Vector<i64, $t>;
            fn neg(self) -> Self::Output { 
                -Vector::from(self)
            }
        }
    )*)
}


#[cfg(test)]
mod formal_sum_test {
    use super::Vector;

    #[test]
    fn init_test() {
        let _ = Vector::<i64, &str>::zero();
        let _ = Vector::<i64, &str>::from(["v", "w"]);
        let _ = Vector::<i64, &str>::from(["v"]);
    }

    #[test]
    fn arithmetic_test() {
        let v: Vector<i64, &str> = Vector::from("v");
        let w = Vector::from("w");
        let x = Vector::from("x");
        let s: Vector<i64, &str> = Vector::from(["v", "w"]);
        let t = Vector::from(["w", "x"]);
        let u = Vector::from_iter(["v", "w", "x"].into_iter());

        assert_eq!( v.clone()+w.clone(), s.clone() );
        assert_eq!( s.clone(), u.clone()-x.clone() );
        assert_eq!( s.clone()+t.clone(), u.clone()+w.clone() );
        assert_eq!( v.clone()+w.clone()*2_i64+x.clone()-s.clone()-t.clone(), Vector::zero() );
        assert_eq!( u*2 + w.clone(), [v,w,x,s,t].into_iter().sum() );
    }

    #[test]
    fn coefficient_test() {
        let t = Vector::<i64, &str>::from(vec!["w", "x", "y"]);
        let u = Vector::from(vec!["v", "w", "x"]);
        let s = t - u;

        assert_eq!(s[&"v"], -1);
        assert_eq!(s[&"y"], 1);
    }

    #[derive(Clone, PartialEq, Debug, Eq, PartialOrd, Ord)]
    struct Object {
        data: String
    }

    fn obj(data: &str) -> Object {
        Object { data: data.to_string() }
    }

    formal_sum_impl!{ Object }

    #[test]
    fn macro_test() {

        let s = 3 * obj("a") + obj("b");
        let t = -obj("a") + obj("b");
        assert_eq!( t*3+s, 4*obj("b") );
    }

}
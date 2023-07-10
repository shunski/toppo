use crate::commutative::PID;
use std::marker::PhantomData;
use std::ptr::NonNull;

pub struct Matrix<T: PID> {
    ptr: NonNull<T>,
    _marker: PhantomData<T>,
    
    size: (usize, usize),
    offsets: (usize, usize),
}

pub struct SubMatrix<T: PID>([T]);

pub mod matrix_impl;
mod matrix_test;

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
        let elementary_divisors: Vec<_> = elementary_divisors.into_iter().filter(|&(div, _)| div != Coeff::one()).collect();
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
impl<Coeff: PID> Display for Module<Coeff> {
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
        if let Some(&(div, mul)) = self.elementary_divisors.first() {
            if mul==1 {
                write!(f, "{}_{}", Coeff::symbol(), div)?;
            } else {
                write!(f, "({}_{})^{}", Coeff::symbol(), div, mul)?;
            }
        } else {
            return write!(f, "")
        }
        for &(divisor, multiplicity) in self.elementary_divisors.iter().skip(1) {
            if multiplicity==1 {
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
        let m = Module::<i128>::new(0,Vec::new());
        assert_eq!(m.to_string(), "0");

        let m = Module::<i128>::new(1,Vec::new());
        assert_eq!(m.to_string(), "Z");

        let m = Module::<i128>::new(0,vec![(2,1)]);
        assert_eq!(m.to_string(), "Z_2");

        let m = Module::<i128>::new(1,vec![(2,2)]);
        assert_eq!(m.to_string(), "Z + (Z_2)^2");

        let m = Module::<i128>::new(2,vec![(2,2), (8,1), (24, 2)]);
        assert_eq!(m.to_string(), "Z^2 + (Z_2)^2 + Z_8 + (Z_24)^2");

        let m = Module::<i128>::new(2,vec![(1,3), (3,2), (24, 2)]);
        assert_eq!(m.to_string(), "Z^2 + (Z_3)^2 + (Z_24)^2");
    }
}


#[derive(Debug, Clone)]
pub struct FormalSum<BasisT: Clone + PartialEq> {
    // the i-th element in 'coeffs' represents the coefficient of the i-the element of 'summands'. 
    // Hence 'summands' and ''coeffs' must be of the same length.
    summands: Vec<BasisT>,
    coeffs: Vec<i128>,
}

impl<BasisT: Clone + PartialEq> FormalSum<BasisT> {
    pub fn zero() -> Self {
        FormalSum { 
            summands: vec![], 
            coeffs: vec![]
        }
    }

    pub fn from(elem: BasisT) -> Self {
        FormalSum { 
            summands: vec![elem], 
            coeffs: vec![1]
        }
    }

    pub fn from_vec( summands: Vec<BasisT> ) -> Self {
        let len = summands.len();
        FormalSum {
            summands: summands,
            coeffs: vec![ 1; len ],
        }
    }
}

impl<'a, BasisT: Clone + PartialEq> FormalSum<BasisT> {
    pub fn iter(&'a self) -> std::iter::Zip<std::slice::Iter<'a, i128>, std::slice::Iter<'a, BasisT>> {
        self.coeffs.iter().zip( self.summands.iter() )
    }

    pub fn iter_mut(&'a mut self) -> std::iter::Zip<std::slice::IterMut<'a, i128>, std::slice::IterMut<'a, BasisT>> {
        self.coeffs.iter_mut().zip( self.summands.iter_mut() )
    }
}

impl<BasisT: Clone + PartialEq> FormalSum<BasisT> {
    pub fn into_iter(self) -> std::iter::Zip<std::vec::IntoIter<i128>, std::vec::IntoIter<BasisT>> {
        self.coeffs.into_iter().zip( self.summands.into_iter() )
    }
}

impl<BasisT: Clone + PartialEq> std::ops::Add for FormalSum<BasisT> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        for (j, s) in rhs.summands.into_iter().enumerate() {
            match self.summands.iter().position(|summand| summand == &s) {
                Some(i) => self.coeffs[i] += rhs.coeffs[j],
                None => {
                    self.summands.push( s );
                    self.coeffs.push( rhs.coeffs[j] );
                },
            };
        };

        // reduce zero elements
        let mut zero_pos = self.coeffs.iter().enumerate().filter(|&(_, &coeff)| coeff==0 ).map(|(i, _)| i).collect::<Vec<_>>();
        zero_pos.reverse();
        for i in zero_pos {
            self.summands.remove(i);
            self.coeffs.remove(i);
        }

        self
    }
}

impl<BasisT: Clone + PartialEq> std::ops::Sub for FormalSum<BasisT> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<BasisT: Clone + PartialEq> std::ops::Neg for FormalSum<BasisT> {
    type Output = Self;
    fn neg(mut self) -> Self::Output {
        self.coeffs.iter_mut().for_each(|coeff| *coeff *= -1);
        self
    }
}

impl<BasisT: Clone + PartialEq> std::ops::Mul<FormalSum<BasisT>> for i128 {
    type Output = FormalSum<BasisT>;
    fn mul(self, mut rhs: FormalSum<BasisT>) -> Self::Output {
        if self == 0 { return FormalSum::zero(); };
        rhs.coeffs.iter_mut().for_each(|coeff| *coeff *= self);
        rhs
    }
}

impl<BasisT: Clone + PartialEq> std::ops::Add<BasisT> for FormalSum<BasisT> {
    type Output = Self;
    fn add(self, rhs: BasisT) -> Self::Output {
        self + FormalSum::from( rhs )
    }
}

impl<BasisT: Clone+PartialEq> std::ops::Sub<BasisT> for FormalSum<BasisT> {
    type Output = Self;
    fn sub(self, rhs: BasisT) -> Self::Output {
        self - FormalSum::from( rhs )
    }
}

impl<BasisT: Clone + PartialEq> std::cmp::PartialEq for FormalSum<BasisT> {
    fn eq(&self, other: &Self) -> bool {
        for (i, x) in self.iter() {
            if let Some((j,_)) = other.iter().find(|(_, y)| x == *y ) {
                if i!=j { return false; } 
            } else {
                return false;
            }
        }

        true
    }
}


impl<BasisT: Clone+PartialEq> std::iter::Sum for FormalSum<BasisT> {
    fn sum<I: std::iter::Iterator<Item=Self>>(iter: I) -> Self {
        iter.fold(FormalSum::zero(), |accum, s| accum+s)
    }
}


// index trait for retrieving coefficients
impl<BasisT: Clone+PartialEq> std::ops::Index<&BasisT> for FormalSum<BasisT> {
    type Output = i128;
    fn index(&self, index: &BasisT) -> &Self::Output {
        &self.iter().find(|(_, x)| *x==index ).map(|(i,_)| i).unwrap()
    }
}


impl<BasisT: Clone+PartialEq+Display > std::fmt::Display for FormalSum<BasisT> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // if the members are empty, then print zero and return
        if self.coeffs.is_empty() {
            return write!(f, "0");
        }

        // else we print elements
        let mut iter = self.iter();

        // printing the initial element
        let (&i, x) = iter.next().unwrap();
        if i==1 {
            write!(f, "{}", x)?;
        } else if i==-1 {
            write!(f, "-{}", x)?;
        } else {
            write!(f, "{}·{}", i, x)?;
        };

        // printing rest of the elements
        for (&i, x) in iter {
            if i==1 {
                write!(f, " + {}", x)?;
            } else if i==-1 {
                write!(f, " - {}", x)?;
            } else if i<0 {
                write!(f, " - {}·{}", -i, x)?;
            } else {
                write!(f, " + {}·{}", i, x)?;
            };
        };
        write!(f, "")
    }
}


#[macro_export]
macro_rules! formal_sum_impl {
    ($($t:ty)*) => ($(
        impl std::ops::Add for $t {
            type Output = FormalSum<$t>;
            fn add(self, rhs: Self) -> Self::Output { 
                FormalSum::from(self) + FormalSum::from(rhs)
            }
        }

        impl std::ops::Add<FormalSum<$t>> for $t {
            type Output = FormalSum<$t>;
            fn add(self, rhs: FormalSum<$t>) -> Self::Output { 
                FormalSum::from(self) + rhs
            }
        }

        impl std::ops::Sub for $t {
            type Output = FormalSum<$t>;
            fn sub(self, rhs: Self) -> Self::Output { 
                FormalSum::from(self) - FormalSum::from(rhs)
            }
        }

        impl std::ops::Sub<FormalSum<$t>> for $t {
            type Output = FormalSum<$t>;
            fn sub(self, rhs: FormalSum<$t>) -> Self::Output { 
                FormalSum::from(self) - rhs
            }
        }

        impl std::ops::Mul<$t> for i128 {
            type Output = FormalSum<$t>;
            fn mul(self, rhs: $t) -> FormalSum<$t> { 
                self * FormalSum::from(rhs)
            }
        }

        impl std::ops::Neg for $t {
            type Output = FormalSum<$t>;
            fn neg(self) -> Self::Output { 
                -FormalSum::from(self)
            }
        }
    )*)
}


#[cfg(test)]
mod formal_sum_test {
    use super::FormalSum;

    #[test]
    fn init_test() {
        FormalSum::<&str>::zero();
        FormalSum::from_vec(vec!["v", "w"]);
        FormalSum::from(vec!["v"]);
    }

    #[test]
    fn arithmetic_test() {
        let v = FormalSum::from("v");
        let w = FormalSum::from("w");
        let x = FormalSum::from("x");
        let s = FormalSum::from_vec(vec!["v", "w"]);
        let t = FormalSum::from_vec(vec!["w", "x"]);
        let u = FormalSum::from_vec(vec!["v", "w", "x"]);

        assert_eq!( v.clone()+w.clone(), s.clone() );
        assert_eq!( s.clone(), u.clone()-x.clone() );
        assert_eq!( s.clone()+t.clone(), u.clone()+w.clone() );
        assert_eq!( v+2*w+x-s-t, FormalSum::zero() );
    }

    #[test]
    fn coefficient_test() {
        let t = FormalSum::from_vec(vec!["w", "x", "y"]);
        let u = FormalSum::from_vec(vec!["v", "w", "x"]);
        let s = t - u;

        assert_eq!(s[&"v"], -1);
        assert_eq!(s[&"y"], 1);
    }

    #[test]
    fn macro_test() {
        #[derive(Clone, PartialEq, Debug)]
        struct Object {
            data: String
        }

        fn obj(data: &str) -> Object {
            Object { data: data.to_string() }
        }

        formal_sum_impl!{ Object };

        let s = 3 * obj("a") + obj("b");
        let t = -obj("a") + obj("b");
        assert_eq!( 3*t+s, 4*obj("b") );
    }

}
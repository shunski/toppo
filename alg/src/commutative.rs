use std::ops::{Add, Sub, Mul, Div, Neg, AddAssign, SubAssign, MulAssign, DivAssign};
use std::iter::{Sum, Product};
use std::cmp::Eq;

pub trait Zero {
    fn zero() -> Self;
}

macro_rules! zero_impl {
    ($($t:ty)*) => ($(
        impl Zero for $t {
            #[inline]
            fn zero() -> Self { 0 as $t }
        }
    )*)
}

zero_impl! { isize i8 i16 i32 i64 i128 }

pub trait One {
    fn one() -> Self;
}

macro_rules! one_impl {
    ($($t:ty)*) => ($(
        impl One for $t {
            #[inline]
            fn one() -> Self { 1 as $t }
        }
    )*)
}

one_impl! { isize i8 i16 i32 i64 i128 }

pub trait Symbol {
    fn symbol()->String;
}

macro_rules! integers_symbol_impl {
    ($($t:ty)*) => ($(
        impl Symbol for $t {
            #[inline]
            fn symbol() -> String { "Z".to_string() }
        }
    )*)
}

macro_rules! reals_symbol_impl {
    ($($t:ty)*) => ($(
        impl Symbol for $t {
            #[inline]
            fn symbol() -> String { "R".to_string() }
        }
    )*)
}

integers_symbol_impl! { isize i8 i16 i32 i64 i128 }
reals_symbol_impl! { f64 }


// Rational numbers
#[derive(PartialEq, Eq, Copy, Clone)]
pub struct Rational {
    numerator: u64,
    denumerator: u64,
    sign: bool,
}

#[macro_export]
macro_rules! rational {
    ($a:expr; $b:expr) => { {
        Rational::new($a as i64, $b as i64)
    } };

    ($a:expr) => {{
        Rational::new($a,1)
    }};
}

mod rational_impl;

// Real numbers
macro_rules! reals_impl {
    ($($ty:ident)*) => ($(

        impl Zero for $ty {
            #[inline]
            fn zero() -> Self { 0.0 }
        }

        impl One for $ty {
            #[inline]
            fn one() -> Self { 1.0 }
        }

        // impl std::fmt::Display for $ty {
        //     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        //         std::fmt::Display::fmt(&self, f)
        //     }
        // }
    )*)
}
reals_impl!{ f64 }


// PID implementations
pub trait PID: 
    Zero + One
    + Add + Sub + Mul + Neg
    + AddAssign + SubAssign + MulAssign
    + Add<Output=Self> + Sub<Output=Self> + Mul<Output=Self> + Neg<Output=Self>
    + PartialEq + Sum + Product
    + Sized + Copy 
    + std::fmt::Debug
    + std::fmt::Display
    + Symbol
{
    fn is_prime(self) -> bool;
    fn n_prime_factors(self) -> usize;
    fn euclid_div(self, d: Self) -> (Self, Self);
    fn simplify(self) -> (Self, Self);
    fn divides(self, n: Self) -> bool {
        n.euclid_div(self).1 == Self::zero()
    }
}

macro_rules! pid_impl {
    ($($t:ty)*) => ($(
        impl PID for $t {
            fn is_prime(self) -> bool {
                if self == Self::zero() || self == Self::one() { return false; }
                (2..self).all(|x| self/x != 0)
            }

            fn n_prime_factors(self) -> usize {
                if self.is_prime() { return 1; }
                let mut count = 0;
                (0..self).filter(|&x| x.is_prime()).filter(|&x| self/x == 0).for_each(|x|{
                    let mut num = self;
                    while num/x==0 {
                        count += 1;
                        num /= x;
                    }
                });
                count
            }

            fn euclid_div(self, d: Self) -> (Self, Self) {
                let mut q = self/d;
                let mut r = self%d;
                if r < 0 {
                    if d > 0 {
                        r += d;
                        q -= 1;
                    } else {
                        r -= d;
                        q += 1;
                    }
                }
                (q, r)
            }

            fn simplify(self) -> (Self, Self) {
                let sign = if self < 0 {
                    -1 
                } else {
                    1
                };
                (sign, sign*self)
            }
        }
    )*)
}

pid_impl!{ isize i8 i16 i32 i64 i128 }

macro_rules! pid_impl_for_fields {
    ($($t:ty)*) => ($(
        impl PID for $t {
            fn is_prime(self) -> bool { true }
            fn n_prime_factors(self) -> usize { 1 }
            fn euclid_div(self, d: Self) -> (Self, Self) { (self/d, Self::zero()) }
            fn simplify(self) -> (Self, Self) { 
                (<$t>::one(), self)
            }
        }
    )*)
}

pid_impl_for_fields!{ Rational f64 }

#[allow(unused)]
fn gcd<T: PID>(a: T, b: T) -> T {
    if a == T::zero() { return b.simplify().1 };
    let (_, r) = b.euclid_div(a);
    gcd(r, a)
} 

pub fn bezout_identity<T: PID>(a: T, b: T) -> (T, T, T) { // output = (gcd(a, b), x, y)
    bezout_identity_rec(a, b, T::one(), T::zero(), T::zero(), T::one())
}

fn bezout_identity_rec<T: PID>(a: T, b: T, x1: T, x2: T, y1: T, y2: T) -> (T, T, T) {
    if a == T::zero() { 
        let (unit, b) = b.simplify();
        return (b, x2.euclid_div(unit).0, y2.euclid_div(unit).0)
    };
    let (q, r) = b.euclid_div(a);
    bezout_identity_rec(r, a, x2 - x1*q, x1, y2 - y1*q, y1)
}

#[cfg(test)]
mod pid_test {
    use crate::commutative::*;
    use rand::Rng;

    #[test]
    fn gcd_test() {
        assert_eq!(gcd(12, 15), 3);
        assert_eq!(gcd(2, 6), 2);
        assert_eq!(gcd(-12, 15), 3);
        assert_eq!(gcd(2, -6), 2);
        assert_eq!(gcd(-15, -12), 3);
        assert_eq!(gcd(-2, -6), 2);

        let mut rng = rand::thread_rng();
        let n_test_cases = 1000;
        for _ in 0..n_test_cases {
            let a: i64 = rng.gen(); let b: i64 = rng.gen();
            let g = gcd(a, b);
            assert!(g.divides(a));
            assert!(g.divides(b));
            assert!(g >= 0);
        }
    }

    #[test]
    fn bezout_identity_test() {
        assert_eq!(bezout_identity(14, 6), (2, 1, -2));

        assert_eq!(bezout_identity(6, 33), (3, -5, 1));

        assert_eq!(bezout_identity(1, 0), (1, 1, 0));

        assert_eq!(bezout_identity(1, 1), (1, 1, 0));

        assert_eq!(bezout_identity(0, 1), (1, 0, 1));

        assert_eq!(bezout_identity(1, -1), (1, 1, 0));

        let mut rng = rand::thread_rng();
        let n_test_cases = 1000;
        for _ in 0..n_test_cases {
            let range = -10000..10000;
            let a: i64 = rng.gen_range(range.clone()); 
            let b: i64 = rng.gen_range(range);
            let (g, x, y) = bezout_identity(a, b);
            assert_eq!(g, gcd(a, b));
            assert_eq!(a.euclid_div(g).0 * x + b.euclid_div(g).0 * y, 1);
        }
    }

    #[test]
    fn divides() {
        assert!(1.divides(0));
        assert!(1.divides(1));
        assert!(1.divides(-1));
        assert!((-1).divides(1));
        assert!((-1).divides(-1));
        assert!(1.divides(2));
        assert!(3.divides(27));
        assert!(!27.divides(3));
        assert!(!2.divides(5));
    }
}

pub  trait Field: PID + Div + DivAssign + Div<Output=Self> {}

macro_rules! field_impl {
    ($($t:ty)*) => ($(
        impl Field for $t {}
    )*)
}

field_impl!{ Rational f64 }
use std::ops::Add;
use std::ops::Sub;
use std::ops::Mul;
use std::ops::Div;
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

zero_impl! { isize i8 i16 i32 i64 i128 f32 f64 }

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

one_impl! { isize i8 i16 i32 i64 i128 f32 f64 }

pub trait NegativeOne {
    fn negative_one() -> Self;
}

macro_rules! negative_one_impl {
    ($($t:ty)*) => ($(
        impl NegativeOne for $t {
            #[inline]
            fn negative_one() -> Self { -1 as $t }
        }
    )*)
}

negative_one_impl! { isize i8 i16 i32 i64 i128 f32 f64 }

pub mod rational;
pub trait PID: 
    Zero + One + NegativeOne
    + Add + Sub + Mul
    + Add<Output=Self> + Sub<Output=Self> + Mul<Output=Self>
    + Eq + Sized + Copy 
    + std::fmt::Debug
    + std::fmt::Display
{
    fn is_prime(self) -> bool;
    fn n_prime_factors(self) -> usize;
    fn euclid_div(self, d: Self) -> (Self, Self);
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
                (self/d, self%d)
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
        }
    )*)
}

pid_impl_for_fields!{ rational::Rational }

#[allow(unused)]
fn gcd<T: PID>(a: T, b: T) -> T {
    if a == T::zero() { return b }
    let (_, r) = b.euclid_div(a);
    gcd(r, a)
} 

pub fn bezout_identity<T: PID>(a: T, b: T) -> (T, T, T) { // output = (gcd(a, b), x, y)
    bezout_identity_rec(a, b, T::one(), T::zero(), T::zero(), T::one())
}

fn bezout_identity_rec<T: PID>(a: T, b: T, x1: T, x2: T, y1: T, y2: T) -> (T, T, T) {
    if a == T::zero() { return (b, x2, y2) }
    let (q, r) = b.euclid_div(a);
    bezout_identity_rec(r, a, x2 - x1*q, x1, y2 - y1*q, y1)
}

#[cfg(test)]
mod pid_test {
    use crate::commutative::*;
    #[test]
    fn gcd_test() {
        assert_eq!(gcd(12, 15), 3);
        assert_eq!(gcd(2, 6), 2);
        assert_eq!(gcd(-12, 15), 3);
        assert_eq!(gcd(2, -6), 2);
        assert_eq!(gcd(-12, -15), -3);
        assert_eq!(gcd(-2, -6), -2);
    }

    #[test]
    fn bezout_identity_test() {
        assert_eq!(bezout_identity(14, 6), (2, 1, -2));

        assert_eq!(bezout_identity(6, 33), (3, -5, 1));

        let a: i64 = -1284719872; let b: i64 = -1293816;
        let (g, x, y) = bezout_identity(a, b);
        assert_eq!(g, gcd(a, b));
        assert_eq!(a*x+b*y, g);
    }

    #[test]
    fn divides() {
        assert!(1.divides(0));
        assert!(1.divides(1));
        assert!(1.divides(2));
        assert!(3.divides(27));
        assert!(!27.divides(3));
        assert!(!2.divides(5));
    }
}

pub struct Module<T: PID> {
    _val: T,
}

pub  trait Field: PID + Div + Div<Output=Self> {}

macro_rules! field_impl {
    ($($t:ty)*) => ($(
        impl Field for $t {}
    )*)
}

field_impl!{ rational::Rational }
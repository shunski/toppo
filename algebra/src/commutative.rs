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

pub mod rational;
pub trait PID: 
    Zero + One
    + Add + Sub + Mul 
    + Add<Output=Self> + Sub<Output=Self> + Mul <Output=Self>
    + Eq + Sized + Copy 
    + std::fmt::Debug
{}

macro_rules! pid_impl {
    ($($t:ty)*) => ($(
        impl PID for $t {}
    )*)
}

pid_impl!{ isize i8 i16 i32 i64 i128 rational::Rational }


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
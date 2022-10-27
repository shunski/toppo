use std::ops::Add;
use std::ops::Sub;
use std::ops::Mul;
use std::ops::Div;
use std::cmp::Ordering;
use num::integer::gcd;

use crate::commutative::Zero;
use crate::commutative::One;
use crate::commutative::NegativeOne;
use crate::commutative::Symbol;

#[derive(PartialEq, Eq, Copy, Clone)]
pub struct Rational {
    numerator: u64,
    denumerator: u64,
    sign: bool,
}

impl Rational {
    pub fn new(n: i64, d: i64)-> Rational {
        if d == 0 { panic!("Denumerators cannot be zero!") };

        let mut r = Rational {
            numerator: n.abs() as u64,
            denumerator: d.abs() as u64,
            sign: (n>=0) == (d>=0),
        };

        r.optimize();
        r
    }

    pub fn from(n: i64) -> Rational {
        Rational::new(n, 1)
    }

    fn optimize(&mut self) {

        if self.numerator == 0 {
            self.denumerator = 1;
            self.sign = true;
            return;
        }


        let gcd = gcd(self.numerator, self.denumerator);
        self.numerator = self.numerator / gcd;
        self.denumerator = self.denumerator / gcd;
    }

    pub fn abs(mut self) -> Rational {
        self.sign = true;
        self
    }

    pub fn get_numerator(&self) -> usize {
        self.numerator as usize
    }

    pub fn get_denumerator(&self) -> usize {
        self.denumerator as usize
    }
}

impl Ord for Rational {
    fn cmp(&self, rhs: &Self) -> Ordering {
        let gcd_of_n = gcd(self.numerator, rhs.numerator);
        let gcd_of_d = gcd(self.denumerator, rhs.denumerator);
        let a = (self.numerator/gcd_of_n) * (rhs.denumerator/gcd_of_d);
        let b = (rhs.numerator/gcd_of_n) * (self.denumerator/gcd_of_d);
        if self.sign & rhs.sign { 
            a.cmp(&b) 
        } else if !self.sign & !rhs.sign {
            b.cmp(&a)
        } else if self.sign & !rhs.sign {
            Ordering::Greater
        } else {
            Ordering::Less
        }
    }
}

impl PartialOrd for Rational {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}

pub trait RationalFromInteger {
    fn as_rational(&self) -> Rational;
    fn over(&self, d: Self) -> Rational;
}

macro_rules! rational_from_impl {
    ($($t:ty)*) => ($(
        impl RationalFromInteger for $t {
            #[inline]
            fn as_rational(&self) -> Rational { 
                let n = *self as i64;
                Rational::new(n, 1)
            }
            fn over(&self, d: Self) -> Rational { 
                let n = *self as i64;
                let d = d as i64;
                Rational::new(n, d)
             }
        }
    )*)
}

rational_from_impl! { isize i8 i16 i32 i64 usize }


impl Add for Rational {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let n: u64;
        let s: bool;
        let gcd_of_d = gcd(self.denumerator, rhs.denumerator);
        if self.sign == rhs.sign {
            n = (self.numerator * (rhs.denumerator/gcd_of_d)) + ((self.denumerator/gcd_of_d) * rhs.numerator);
            s = self.sign;
        } else if self.abs() >= rhs.abs() {
            n = (self.numerator * (rhs.denumerator/gcd_of_d)) - ((self.denumerator/gcd_of_d) * rhs.numerator);
            s = self.sign;
        } else {
            n = ((self.denumerator/gcd_of_d) * rhs.numerator) - (self.numerator * (rhs.denumerator/gcd_of_d));
            s = rhs.sign;
        }
        let mut r = Rational { 
            numerator: n, 
            denumerator: self.denumerator * (rhs.denumerator/gcd_of_d),
            sign: s,
        };
        r.optimize();
        r
    }
}

impl Sub for Rational {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + rhs * (-1).as_rational()
    }
}

impl Mul for Rational {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut r = Rational { 
            numerator: self.numerator * rhs.numerator, 
            denumerator: self.denumerator * rhs.denumerator,
            sign: self.sign == rhs.sign,
        };
        r.optimize();
        r
    }
}

impl Div for Rational {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        if rhs == 0.as_rational() { panic!("Cannot divide by zero.") };
        let mut r = Rational { 
            numerator: self.numerator * rhs.denumerator, 
            denumerator: self.denumerator * rhs.numerator,
            sign: self.sign == rhs.sign,
        };
        r.optimize();
        r
    }
}

impl Zero for Rational {
    fn zero() -> Self {
        Rational::new(0, 1)
    }
}

impl One for Rational {
    fn one() -> Self {
        Rational::new(1, 1)
    }
}

impl NegativeOne for Rational {
    fn negative_one() -> Self {
        Rational::new(-1, 1)
    }
}

impl Symbol for Rational {
    fn symbol() -> String {
        "Q".to_string()
    }
}

impl std::fmt::Debug for Rational {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let sign = if self.sign {
            ""
        } else {
            "-"
        };
        if self.denumerator == 1 {
            write!(f, "{}{}", sign, self.numerator)
        } else {
            write!(f, "{}{}/{}", sign, self.numerator, self.denumerator)
        }
    }
}

impl std::fmt::Display for Rational {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let sign = if self.sign {
            ""
        } else {
            "-"
        };
        if self.denumerator == 1 {
            write!(f, "{}{}", sign, self.numerator)
        } else {
            write!(f, "{}{}/{}", sign, self.numerator, self.denumerator)
        }
    }
}

#[cfg(test)]
mod rational_test {
    #[test] 
    fn init_tests() {
        use crate::commutative::rational::*;

        assert_eq!(Rational::new(-1, 1), Rational{ numerator: 1, denumerator: 1, sign: false });
        assert_eq!(Rational::new(0, 10), Rational{ numerator: 0, denumerator: 1, sign: true });
        assert_eq!(Rational::new(9, -6), Rational{ numerator: 3, denumerator: 2, sign: false });
        assert_eq!(Rational::new(-9, -3), Rational{ numerator: 3, denumerator: 1, sign: true });

        assert_eq!(Rational::new(-4, -6), Rational::new(2, 3));
        assert_eq!(Rational::new(15, -12), Rational::new(5, 4) * (-1).as_rational());
        assert_eq!(Rational::new(-4, 6).abs(), Rational::new(2, 3));

        assert_eq!(Rational::new(-4, 6), (-2).over(3));
        assert_eq!(Rational::zero(), (0).over(1));
    }

    #[test] 
    fn arithmetic_test() {
        use crate::commutative::rational::*;

        let r1 = Rational::new(-2, 2);
        assert_eq!(r1, (-1).as_rational());
        assert_eq!(Rational::zero()-Rational::one(), (-1).as_rational());

        let r2 = Rational::from(-2);
        assert_eq!(r2, Rational::new(-2, 1));

        assert_eq!(r2 < r1, true);

        let r3 = 3.as_rational() / r2;
        assert_eq!(r3, Rational::new(-3, 2));

        let r4 = r3 + Rational::new(2, 3);
        assert_eq!(r4, Rational::new(-5, 6));

        let r5 = r3 * r4;
        assert_eq!(r5, Rational::new(5, 4));

        let r6 = r5 - Rational::new(2, 3);
        assert_eq!(r6, Rational::new( 7, 12));

        let r7 = 0.over(3); 
        assert_eq!(r7 * r7, Rational::new( 0, 1)); 
        assert_eq!(r6+r7, r6);

    }

    #[test]
    #[should_panic]
    fn init_with_denumerator_zero() {
        use crate::commutative::rational::*;

        let _a = Rational::new(100, 0);
    }

    #[test]
    #[should_panic]
    fn division_by_zero_1() {
        use crate::commutative::rational::*;

        let a = 3.as_rational();
        let b = a - a;
        let b = b * Rational::from(-1000);
        let _ = a / b;
    }
}
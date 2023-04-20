use super::Rational;

use std::iter::Product;
use std::iter::Sum;
use std::ops::Add;
use std::ops::Neg;
use std::ops::Sub;
use std::ops::Mul;
use std::ops::Div;
use std::ops::AddAssign;
use std::ops::SubAssign;
use std::ops::MulAssign;
use std::ops::DivAssign;
use std::cmp::Ordering;
use num::integer::gcd;

use crate::commutative::Zero;
use crate::commutative::One;
use crate::commutative::Symbol;

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
        self + rhs * -Self::one()
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
        if rhs == Self::zero() { panic!("Cannot divide by zero.") };
        let mut r = Rational { 
            numerator: self.numerator * rhs.denumerator, 
            denumerator: self.denumerator * rhs.numerator,
            sign: self.sign == rhs.sign,
        };
        r.optimize();
        r
    }
}

macro_rules! op_assign_impl {
    ($t:ty, $op_assign_trait:ident, $op_assign:ident, $op_trait:ident, $op:ident) => {
        impl $op_assign_trait for $t {
            fn $op_assign(&mut self, rhs: $t) {
                *self =  $op_trait::$op( *self, rhs);
            }
        }
    };
}

op_assign_impl!( Rational, AddAssign, add_assign, Add, add );
op_assign_impl!( Rational, SubAssign, sub_assign, Sub, sub );
op_assign_impl!( Rational, MulAssign, mul_assign, Mul, mul );
op_assign_impl!( Rational, DivAssign, div_assign, Div, div );

impl Neg for Rational {
    type Output = Rational;
    fn neg(mut self) -> Self::Output {
        self.sign = !self.sign;
        self
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

impl Symbol for Rational {
    fn symbol() -> String {
        "Q".to_string()
    }
}


impl Sum for Rational {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Rational::zero(), |accum, a| accum + a)
    }
}

impl Product for Rational {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Rational::zero(), |accum, a| accum * a)
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
    use crate::commutative::{Rational, Zero, One};
    use crate::rational;
    #[test] 
    fn init_tests() {
        assert_eq!(Rational::new(-1, 1), Rational{ numerator: 1, denumerator: 1, sign: false });
        assert_eq!(Rational::new(0, 10), Rational{ numerator: 0, denumerator: 1, sign: true });
        assert_eq!(Rational::new(9, -6), Rational{ numerator: 3, denumerator: 2, sign: false });
        assert_eq!(Rational::new(-9, -3), Rational{ numerator: 3, denumerator: 1, sign: true });

        assert_eq!(Rational::new(-4, -6), Rational::new(2, 3));
        assert_eq!(Rational::new(15, -12), Rational::new(5, 4) * rational!(-1));
        assert_eq!(Rational::new(-4, 6).abs(), Rational::new(2, 3));

        assert_eq!(Rational::new(-4, 6), rational!(-2; 3));
        assert_eq!(Rational::zero(), rational!(0; 1));
    }

    #[test] 
    fn arithmetic_test() {
        let r1 = Rational::new(-2, 2);
        assert_eq!(r1, rational!(-1));
        assert_eq!(Rational::zero()-Rational::one(), rational!(-1));

        let r2 = Rational::from(-2);
        assert_eq!(r2, Rational::new(-2, 1));

        assert_eq!(r2 < r1, true);

        let r3 = rational!(3) / r2;
        assert_eq!(r3, Rational::new(-3, 2));

        let r4 = r3 + Rational::new(2, 3);
        assert_eq!(r4, Rational::new(-5, 6));

        let r5 = r3 * r4;
        assert_eq!(r5, Rational::new(5, 4));

        let r6 = r5 - Rational::new(2, 3);
        assert_eq!(r6, Rational::new( 7, 12));

        let r7 = rational!(0;3); 
        assert_eq!(r7 * r7, Rational::new( 0, 1)); 
        assert_eq!(r6+r7, r6);

    }

    #[test]
    #[should_panic]
    fn init_with_denumerator_zero() {
        let _a = Rational::new(100, 0);
    }

    #[test]
    #[should_panic]
    fn division_by_zero_1() {
        let a = rational!(3);
        let b = a - a;
        let b = b * Rational::from(-1000);
        let _ = a / b;
    }
}
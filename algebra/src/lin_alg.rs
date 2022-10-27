use crate::commutative::PID;

pub mod matrix;

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
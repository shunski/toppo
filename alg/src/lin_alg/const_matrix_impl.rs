use super::{ConstVector, ConstMatrix};
use crate::commutative::{PID, Field};

use std::ops::{
    Deref,
    DerefMut,
    Add,
    AddAssign,
    Sub,
    SubAssign,
    Neg,
    Mul,
    MulAssign,
    Div,
    DivAssign,
};

impl<T: PID, const N: usize> Deref for ConstVector<T, N> {
    type Target = [T; N];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: PID, const N: usize> DerefMut for ConstVector<T, N> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T: PID, const N_ROW: usize, const N_COL: usize> Deref for ConstMatrix<T, N_ROW, N_COL> {
    type Target = [[T; N_COL]; N_ROW];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: PID, const N_ROW: usize, const N_COL: usize> DerefMut for ConstMatrix<T, N_ROW, N_COL> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}


// ToDo: implemenet the into_iter for the ConstVector types, and implement this trat for non-copiable 'T'.
// ToDo: integrate the implimentations using macros
impl <T: PID + Copy, const N: usize> Add for ConstVector<T, N> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
} 

// ToDo: implemenet the into_iter for the ConstVector types, and implement this trat for non-copiable 'T'.
impl <T: PID + Copy, const N: usize> Sub for ConstVector<T, N> {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
} 

// ToDo: integrate the implimentations using macros
impl <T: PID, const N: usize> Mul<T> for ConstVector<T, N> {
    type Output = Self;
    fn mul(mut self, rhs: T) -> Self::Output {
        self *= rhs;
        self
    }
} 

impl <T: PID, const N_ROW: usize, const N_COL: usize> Mul<T> for ConstMatrix<T, N_ROW, N_COL> {
    type Output = Self;
    fn mul(mut self, rhs: T) -> Self::Output {
        self *= rhs;
        self
    }
} 

// ToDo: implemenet the into_iter for the ConstVector types, and implement this trat for non-copiable 'T'.
impl <T: Field, const N: usize> Div<T> for ConstVector<T, N> {
    type Output = Self;
    fn div(mut self, rhs: T) -> Self::Output {
        self /= rhs;
        self
    }
} 


impl <T: Field, const N_ROW: usize, const N_COL: usize> Div<T> for ConstMatrix<T, N_ROW, N_COL> {
    type Output = Self;
    fn div(mut self, rhs: T) -> Self::Output {
        self /= rhs;
        self
    }
} 


// ToDo: integrate the implimentations using macros
impl <T: PID + Copy, const N: usize> AddAssign for ConstVector<T, N> {
    fn add_assign(&mut self, rhs: Self) {
        self.iter_mut()
            .zip(rhs.into_iter())
            .for_each(|(x, y)| *x += y );
    }
} 

// ToDo: implemenet the into_iter for the ConstVector types, and implement this trat for non-copiable 'T'.
impl <T: PID + Copy, const N: usize> SubAssign for ConstVector<T, N> {
    fn sub_assign(&mut self, rhs: Self) {
        self.iter_mut()
            .zip(rhs.into_iter())
            .for_each(|(x, y)| *x -= y );
    }
} 

impl<T: PID, const N: usize> MulAssign<T> for ConstVector<T, N> {
    fn mul_assign(&mut self, rhs: T) {
        self.iter_mut()
            .for_each(|x| *x *= rhs.clone() );
    }
}

impl<T: Field, const N: usize> DivAssign<T> for ConstVector<T, N> {
    fn div_assign(&mut self, rhs: T) {
        self.iter_mut()
            .for_each(|x| *x /= rhs.clone() );
    }
}


// ToDo: Implement the into_iter() for ConstVector and remove the clone() in the implementation
impl<T: PID, const N: usize> Neg for ConstVector<T, N> {
    type Output = Self;
    fn neg(mut self) -> Self::Output {
        self.iter_mut()
            .for_each(|x| *x = x.clone().neg());
        self
    }
}

// ToDo: integrate the implimentations using macros
impl <T: PID + Copy, const N_ROW: usize, const N_COL: usize> AddAssign for ConstMatrix<T, N_ROW, N_COL> {
    fn add_assign(&mut self, rhs: Self) {
        self.iter_mut()
            .map(|r| r.iter_mut() )
            .flatten()
            .zip(rhs.into_iter().map(|r| r.into_iter() ).flatten())
            .for_each(|(x, y)| *x += y );
    }
} 

// ToDo: implemenet the into_iter for the ConstVector types, and implement this trait for non-copiable 'T'.
impl <T: PID + Copy, const N_ROW: usize, const N_COL: usize> SubAssign for ConstMatrix<T, N_ROW, N_COL> {
    fn sub_assign(&mut self, rhs: Self) {
        self.iter_mut()
            .map(|r| r.iter_mut() )
            .flatten()
            .zip(rhs.into_iter().map(|r| r.into_iter() ).flatten())
            .for_each(|(x, y)| *x -= y );
    }
} 

impl<T: PID, const N_ROW: usize, const N_COL: usize> MulAssign<T> for ConstMatrix<T, N_ROW, N_COL> {
    fn mul_assign(&mut self, rhs: T) {
        self.iter_mut()
            .map(|x| x.iter_mut() )
            .flatten()
            .for_each(|x| *x *= rhs.clone() );
    }
}

impl<T: Field, const N_ROW: usize, const N_COL: usize> DivAssign<T> for ConstMatrix<T, N_ROW, N_COL> {
    fn div_assign(&mut self, rhs: T) {
        self.iter_mut()
            .map(|x| x.iter_mut() )
            .flatten()
            .for_each(|x| *x /= rhs.clone() );
    }
}


// ToDo: Implement the into_iter() for ConstVector and remove the clone() in the implementation
impl<T: PID, const N_ROW: usize, const N_COL: usize> Neg for ConstMatrix<T, N_ROW, N_COL> {
    type Output = Self;
    fn neg(mut self) -> Self::Output {
        self.iter_mut()
            .map(|x| x.iter_mut() )
            .flatten()
            .for_each(|x| *x = x.clone().neg());
        self
    }
}

impl<T: PID + Copy, const N: usize> Copy for ConstVector<T, N> {}
impl<T: PID + Copy, const N_ROW: usize, const N_COL: usize> Copy for ConstMatrix<T, N_ROW, N_COL> {}

impl<T: PID, const N: usize> From<[T; N]> for ConstVector<T, N> {
    fn from(value: [T; N]) -> Self {
        Self(value)
    }
}

impl <T: PID + Copy, const N: usize> ConstVector<T, N> {
    pub fn zero() -> Self {
        Self([T::zero(); N])
    }
} 

impl <const N: usize> ConstVector<f64, N> {
    pub fn two_norm(self) -> f64 {
        self.into_iter().map(|x| x*x ).sum::<f64>().sqrt()
    }
}

impl <T: PID + Copy, const N_ROW: usize, const N_COL: usize> ConstMatrix<T, N_ROW, N_COL> {
    pub fn zero() -> Self {
        Self([[T::zero(); N_COL]; N_ROW])
    }
} 

impl <T: PID + Copy, const N: usize> ConstMatrix<T, N, N> {
    pub fn identity() -> Self {
        let mut m = Self::zero();
        for i in 0..N {
            m[i][i] = T::one()
        }
        m
    }
} 

impl <T: PID + Copy, const N_ROW: usize, const N_COL: usize> From<[ConstVector<T, N_ROW>; N_COL]> for ConstMatrix<T, N_ROW, N_COL> {
    fn from(value: [ConstVector<T, N_ROW>; N_COL]) -> Self {
        let mut m = Self::zero();
        for j in 0..N_COL {
            for i in 0..N_ROW {
                m[i][j] = value[j][i];
            }
        }
        m
    }
} 

impl <T: PID + Copy, const N_ROW: usize, const N_COL: usize> From<[[T; N_COL]; N_ROW]> for ConstMatrix<T, N_ROW, N_COL> {
    fn from(value: [[T; N_COL]; N_ROW]) -> Self {
        Self(value)
    }
} 

impl<T: PID + Eq, const N: usize> Eq for ConstVector<T, N> {}
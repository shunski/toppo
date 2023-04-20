use crate::commutative::{PID, Field};

use std::alloc;
use std::marker::PhantomData;
use std::mem;
use std::ptr;
use std::fmt;
use std::cmp;
use std::ops;

use super::Matrix;

#[derive(Clone, Copy)]
pub struct SubMatrix<'a, T: PID> {
    pub size: (usize, usize),
    offsets: (usize, usize),
    data: &'a [T],
}

pub struct SubMatrixMut<'a, T: PID> {
    pub size: (usize, usize),
    offsets: (usize, usize),
    data: &'a mut [T],
}

impl<'a, T: PID> SubMatrix<'a, T> {
    pub fn transpose(mut self) -> Self {
        let tmp = (self.offsets.0, self.size.0);
        (self.offsets.0, self.size.0) = (self.offsets.1, self.size.1);
        (self.offsets.1, self.size.1) = tmp;

        self
    }

    pub fn dot(self, rhs: Self) -> T {
        // This function assume that 'self' is a row vector and 'rhs' is a col vector
        debug_assert!(self.size.0 == 1 && rhs.size.1 == 1);

        // This function makes sure that the arguments are of the same size
        debug_assert!(self.size.1 == rhs.size.0);

        (0..self.size.1).map(|i| self.data[i*self.offsets.1] * rhs.data[i*rhs.offsets.0]).sum()
    }

    pub fn as_matrix(&self) -> Matrix<T> {
        let mut out = Matrix::new(self.size.0, self.size.1);
        for i in 0..self.size.0 {
            for j in 0..self.size.1 {
                out[(i,j)] = self[(i, j)];
            }
        }

        out
    }
}

impl<'a, T: PID> SubMatrixMut<'a, T> {
    pub fn swap_rows(&mut self, r1: usize, r2: usize) {
        let mut r1_offset = self.offsets.0 * r1;
        let mut r2_offset = self.offsets.0 * r2;

        for _ in 0..self.size.1 {
            let tmp = self.data[r1_offset];
            self.data[r1_offset] = self.data[r2_offset];
            self.data[r2_offset] = tmp;

            r1_offset+=self.offsets.1; 
            r2_offset+=self.offsets.1;
        }
    }

    pub fn swap_cols(&mut self, c1: usize, c2: usize) {
        let mut c1_offset = self.offsets.1 * c1;
        let mut c2_offset = self.offsets.1 * c2;

        for _ in 0..self.size.0 {
            let tmp = self.data[c1_offset];
            self.data[c1_offset] = self.data[c2_offset];
            self.data[c2_offset] = tmp;

            c1_offset+=self.offsets.0;
            c2_offset+=self.offsets.0;
        }
    }

    pub fn row_operation(mut self, r1: usize, r2: usize, scalar: T) {
        let r2 = SubMatrix{
            size: self.size,
            offsets: self.offsets,
            data: self.data,
        }.sub(r2, ..).as_matrix() * scalar;
        (0..self.size.1).for_each(|i| self[(r1, i)] += r2[(0, i)]);
    }

    pub fn col_operation(mut self, c1: usize, c2: usize, scalar: T) {
        let c2 = SubMatrix{
            size: self.size,
            offsets: self.offsets,
            data: self.data,
        }.sub(.., c2).as_matrix() * scalar;
        (0..self.size.0).for_each(|j| self[(j, c1)] += c2[(j, 0)]);
    }

    pub fn write(mut self, m: Matrix<T>) {
        assert_eq!(self.size, m.size);
        for i in 0..self.size.0 {
            for j in 0..self.size.1 {
                self[(i,j)] = m[(i,j)];
            }
        }
    }
}

use std::ops::RangeBounds;
pub trait ModifiedUsizeRangeBounds {
    fn modified_start_bound(&self)->ops::Bound<&usize>;
    fn modified_end_bound(&self)->ops::Bound<&usize>;
}

impl ModifiedUsizeRangeBounds for usize {
    fn modified_start_bound(&self)->ops::Bound<&usize> {
        ops::Bound::Included(&self)
    }
    fn modified_end_bound(&self)->ops::Bound<&usize> {
        ops::Bound::Included(&self)
    }
}

impl ModifiedUsizeRangeBounds for ops::RangeFull {
    fn modified_start_bound(&self)->ops::Bound<&usize> {
        self.start_bound()
    }
    fn modified_end_bound(&self)->ops::Bound<&usize> {
        self.end_bound()
    }
}

macro_rules! modified_usize_range_bounds_impl {
    ($($range_ty: ident)*) => ($(
        impl ModifiedUsizeRangeBounds for ops::$range_ty<usize> {
            fn modified_start_bound(&self)->ops::Bound<&usize> {
                self.start_bound()
            }
            fn modified_end_bound(&self)->ops::Bound<&usize> {
                self.end_bound()
            }
        }
    )*)
}

modified_usize_range_bounds_impl!( Range RangeFrom RangeInclusive RangeTo RangeToInclusive );

macro_rules! unpack_bounds_for_sub {
    ($range: expr, $size: expr) => {{
        use ops::Bound;
        let start = match $range.modified_start_bound() {
            Bound::Included(&x) => x,
            Bound::Excluded(_) => panic!("input range not valid"),
            Bound::Unbounded => 0,
        };

        let end = match $range.modified_end_bound() {
            Bound::Included(&x) => x+1,
            Bound::Excluded(&x) => x,
            Bound::Unbounded => $size,
        };

        (start, end)
    }}
}

impl<'a, T: PID> SubMatrix<'a, T> {
    pub fn sub(
        &self, 
        rows: impl ModifiedUsizeRangeBounds, 
        cols: impl ModifiedUsizeRangeBounds
    ) -> Self {
        let (row_start, row_end) = unpack_bounds_for_sub!(rows, self.size.0);
        let (col_start, col_end) = unpack_bounds_for_sub!(cols, self.size.1);

        assert!(row_end - row_start <= self.size.0, "row range not valid: specified value is {row_start}..{row_end} but this matrix has row length {}.", self.size.0);
        assert!(col_end - col_start <= self.size.1, "col range not valid: specified value is {col_start}..{col_end} but this matrix has col length {}.", self.size.1);

        let start_offset = self.offsets.0 * row_start + self.offsets.1 * col_start;
        let end_offset = self.offsets.0 * (row_end-1) + self.offsets.1 * (col_end-1) + 1;

        SubMatrix {
            size: (row_end - row_start, col_end - col_start),
            offsets: self.offsets,
            data: &self.data[start_offset..end_offset],
        }
    }
}

#[allow(unused)]
impl<T: PID> SubMatrixMut<'_, T> { 
    pub fn sub(
        &self, 
        rows: impl ModifiedUsizeRangeBounds, 
        cols: impl ModifiedUsizeRangeBounds
    ) -> SubMatrix<'_, T> {
        let sub = SubMatrix {
            size: self.size,
            offsets: self.offsets,
            data: &*self.data,
        };
        sub.sub(rows, cols)
    }
    pub fn sub_mut(
        &mut self, 
        rows: impl ModifiedUsizeRangeBounds, 
        cols: impl ModifiedUsizeRangeBounds
    ) -> SubMatrixMut<'_, T> {
        let (row_start, row_end) = unpack_bounds_for_sub!(rows, self.size.0);
        let (col_start, col_end) = unpack_bounds_for_sub!(cols, self.size.1);

        assert!(row_end - row_start <= self.size.0, "row range not valid: specified value is {row_start}..{row_end} but this matrix has row length {}.", self.size.0);
        assert!(col_end - col_start <= self.size.1, "col range not valid: specified value is {col_start}..{col_end} but this matrix has col length {}.", self.size.1);

        let start_offset = self.offsets.0 * row_start + self.offsets.1 * col_start;
        let end_offset = self.offsets.0 * (row_end-1) + self.offsets.1 * (col_end-1) + 1;

        SubMatrixMut {
            size: (row_end - row_start, col_end - col_start),
            offsets: self.offsets,
            data: &mut self.data[start_offset..end_offset],
        }
    }
}

impl<'a, T: PID> ops::Index<(usize, usize)> for SubMatrix<'a, T> {
    type Output = T;
    fn index(&self, idx: (usize, usize)) -> &'a T {
        &self.data[idx.0 * self.offsets.0 + idx.1 * self.offsets.1]
    }
}

impl<'a, T: PID> ops::Add for SubMatrix<'a, T> {
    type Output = Matrix<T>;
    fn add(self, rhs: Self) -> Self::Output {
        if self.size != rhs.size { panic!(); }
        
        let mut out = Matrix::new(self.size.0, self.size.1);
        for i in 0..self.size.0 {
            for j in 0..self.size.1 {
                out[(i,j)] = self[(i,j)] + rhs[(i,j)];
            }
        }
        out
    }
}

impl<'a, T: PID> ops::Sub for SubMatrix<'a, T> 
{
    type Output = Matrix<T>;
    fn sub(self, rhs: Self) -> Self::Output {
        if self.size != rhs.size { panic!(); }

        let mut out = Matrix::new(self.size.0, self.size.1);
        for i in 0..self.size.0 {
            for j in 0..self.size.1 {
                out[(i,j)] = self[(i,j)] - rhs[(i,j)];
            }
        }
        out
    }
}

impl<'a, T: PID> ops::Mul<Self> for SubMatrix<'a, T> {
    type Output = Matrix<T>;
    fn mul(self, rhs: SubMatrix<T>) -> Self::Output {
        let mut prod = Matrix::new(self.size.0, rhs.size.1);
        for i in 0..prod.size.0 {
            for j in 0..prod.size.1 {
                prod[(i, j)] = self.sub(i, ..).dot(rhs.sub(.., j));
            }
        };
        prod
    }
}

macro_rules! binary_op_between_original_and_sub_matrix_impl {
    ($(($bin_op_tr: ident, $bin_op_fn: ident)) *) => {$(
        impl<'a, T: PID> ops::$bin_op_tr<SubMatrix<'a, T>> for Matrix<T> {
            type Output = Self;
            fn $bin_op_fn(self, rhs: SubMatrix<T>) -> Self {
                let lhs = Matrix::<T>::sub(&self, .., ..);
                ops::$bin_op_tr::$bin_op_fn(lhs, rhs)
            }
        }

        impl<'a, T: PID> ops::$bin_op_tr<Matrix<T>> for SubMatrix<'a, T> {
            type Output = Matrix<T>;
            fn $bin_op_fn(self, rhs: Matrix<T>) -> Matrix<T> {
                let rhs = Matrix::<T>::sub(&rhs, .., ..);
                ops::$bin_op_tr::$bin_op_fn(self, rhs)
            }
        }
    )*}
}

binary_op_between_original_and_sub_matrix_impl!{ 
    (Add, add) 
    (Sub, sub) 
    (Mul, mul) 
}

impl<'a, T: PID> ops::Mul<T> for SubMatrix<'a, T> {
    type Output = Matrix<T>;
    fn mul(self, rhs: T) -> Self::Output {
        let mut out = Matrix::new(self.size.0, self.size.1);
        for i in 0..out.size.0 {
            for j in 0..out.size.1 {
                out[(i,j)] = self[(i,j)] * rhs;
            }
        };
        out
    }
}

macro_rules! scalar_mul_impl {
    ($($pid: ty) *) => ($(
        impl<'a> ops::Mul<SubMatrix<'a, $pid>> for $pid {
            type Output = Matrix<$pid>;
            fn mul(self, rhs: SubMatrix<'a, $pid>) -> Self::Output {
                ops::Mul::mul(rhs, self)
            }
        }
    )*)
}

scalar_mul_impl!{ isize i8 i16 i32 i64 i128 crate::commutative::Rational f64 }

impl<'a, T: Field> ops::Div<T> for SubMatrix<'a, T> {
    type Output = Matrix<T>;
    fn div(self, rhs: T) -> Self::Output {
        let mut out = Matrix::new(self.size.0, self.size.1);
        for i in 0..out.size.0 {
            for j in 0..out.size.1 {
                out[(i,j)] = self[(i,j)] / rhs;
            }
        };
        out
    }
}


impl<'a, T: PID> ops::MulAssign<T> for SubMatrixMut<'a, T> {
    fn mul_assign(&mut self, rhs: T) {
        for i in 0..self.size.0 {
            for j in 0..self.size.1{
                self[(i,j)] *= rhs;
            }
        }
    }
}

impl<'a, T: Field> ops::DivAssign<T> for SubMatrixMut<'a, T> {
    fn div_assign(&mut self, rhs: T) {
        for i in 0..self.size.0 {
            for j in 0..self.size.1{
                self[(i,j)] /= rhs;
            }
        };
    }
}

impl<'a, T: PID> ops::Index<(usize, usize)> for SubMatrixMut<'a, T> {
    type Output = T;
    fn index(&self, idx: (usize, usize)) -> &T {
        &self.data[idx.0 * self.offsets.0 + idx.1 * self.offsets.1]
    }
}

impl<'a, T: PID> ops::IndexMut<(usize, usize)> for SubMatrixMut<'a, T> {
    fn index_mut(&mut self, idx: (usize, usize)) -> &mut T {
        &mut self.data[idx.0 * self.offsets.0 + idx.1 * self.offsets.1]
    }
}

impl<'a, T: PID> cmp::PartialEq for SubMatrix<'a, T> {
    fn eq(&self, rhs: &Self) -> bool {
        if self.size != rhs.size {
            return false;
        }
        for i in 0..self.size.0 {
            for j in 0..self.size.1 {
                if self[(i, j)] != rhs[(i, j)] {
                    return false;
                }
            }
        }
        true
    }
}

#[allow(unused)]
impl<T: PID> Matrix<T> {
    fn new(n: usize, m: usize) -> Self {
        let mut mat = Matrix {
            ptr: ptr::NonNull::dangling(),
            _marker: PhantomData,
            size: (n, m),
            offsets: (m, 1),
        };
        mat.alloc();
        mat
    }

    fn alloc(&mut self) {
        let layout = alloc::Layout::array::<T>(self.size.0 * self.size.1).unwrap();
        
        // Ensure that the new allocation doesn't exceed `isize::MAX` bytes.
        assert!(
            layout.size() <= isize::MAX as usize,
            "Allocation too large"
        );
        
        let new_ptr = unsafe { alloc::alloc(layout) };
        
        // If allocation fails, `new_ptr` will be null, in which case we abort.
        self.ptr = match ptr::NonNull::new(new_ptr as *mut T) {
            Some(p) => p,
            None => alloc::handle_alloc_error(layout),
        };
    }

    fn is_diagonal(&self) -> bool {
        for i in 0..self.size.0 {
            for j in 0..self.size.1 {
                if self[(i,j)] != T::zero() && i!=j {
                    return false;
                };
            };
        };
        true
    }
    
    pub fn swap_rows(&mut self, r1: usize, r2: usize) {
        self.sub_mut(.., ..).swap_rows(r1, r2);
    }

    pub fn swap_cols(&mut self, c1: usize, c2: usize) {
        self.sub_mut(.., ..).swap_cols(c1, c2);
    }

    pub fn row_operation(&mut self, r1: usize, r2: usize, scalar: T) {
        self.sub_mut(.., ..).row_operation(r1, r2, scalar);
    }

    pub fn col_operation(&mut self, c1: usize, c2: usize, scalar: T) {
        self.sub_mut(.., ..).col_operation(c1, c2, scalar);
    }

    fn ptr(&self) -> *mut T {
        self.ptr.as_ptr()
    }

    // needed for the contructor macros
    pub fn from_array(v: Vec<Vec<T>>) -> Self {
        if v.len()==0 || v[0].len()==0 {
            panic!("Inputs for matrix constructor macro is invalid. Matrix must have positive dimension!")
        }
        let size = (v.len(), v[0].len());
        if v.iter().map(|x| x.len()).all(|len| len != size.1) {
            panic!("Inputs for matrix constructor macro is invalid");
        }

        let mut out = Matrix::new(size.0, size.1);
        for i in 0..size.0 {
            for j in 0..size.1 {
                out[(i, j)] = v[i][j];
            }
        }

        out
    }

    pub fn transpose(mut self) -> Self {
        let tmp = (self.size.0, self.offsets.0);
        (self.size.0, self.offsets.0) = (self.size.1, self.offsets.1);
        (self.size.1, self.offsets.1) = tmp;

        self
    }
}


use rand::prelude::*;
use rand::distributions::*;
impl<T> Matrix<T> 
    where T: PID, Standard: Distribution<T>
{
    pub fn random(n: usize, m: usize) -> Matrix<T> {
        let mut out = Matrix::new(n ,m);
        let mut rng = rand::thread_rng();
        for i in 0..out.size.0 {
            for j in 0..out.size.1 {
                out[(i,j)] = rng.gen();
            }
        }
        out
    }
}

impl<'a, T: PID> Matrix<T> {
    pub fn sub(&self, rows: impl ModifiedUsizeRangeBounds, cols: impl ModifiedUsizeRangeBounds) -> SubMatrix<'a, T> {
        let data = unsafe{ std::slice::from_raw_parts(self.ptr(), self.size.0 * self.size.1) };
        let orig_matrix_as_sub = SubMatrix {
            data: data,
            size: self.size,
            offsets: self.offsets,
        };
        orig_matrix_as_sub.sub(rows, cols)
    }

    pub fn sub_mut(&mut self, rows: impl ModifiedUsizeRangeBounds, cols: impl ModifiedUsizeRangeBounds) -> SubMatrixMut<'a, T> {
        let (row_start, row_end) = unpack_bounds_for_sub!(rows, self.size.0);
        let (col_start, col_end) = unpack_bounds_for_sub!(cols, self.size.1);
        let start_offset =  row_start * self.offsets.0 + col_start * self.offsets.1;
        let len = (row_end-row_start) * self.offsets.0 + (col_end-col_start) * self.offsets.1;

        let data = unsafe{ std::slice::from_raw_parts_mut(
            self.ptr().add(start_offset), 
            len
        )};
        
        SubMatrixMut {
            data: data,
            size: (row_end-row_start, col_end-col_start),
            offsets: self.offsets,
        }
    }
}


impl<T: PID> Drop for Matrix<T> {
    fn drop(&mut self) {
        let elem_size = mem::size_of::<T>();

        if elem_size != 0 {
            unsafe {
                alloc::dealloc(
                    self.ptr.as_ptr() as *mut u8,
                    alloc::Layout::array::<T>(self.size.0 * self.size.1).unwrap(),
                );
            }
        }
    }
}


// immutable element access
impl<T: PID> ops::Index<(usize, usize)> for Matrix<T> {
    type Output = T;
    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        if i>=self.size.0 || j>=self.size.1 { 
            panic!("Idx out of bounds: index: ({}, {}) while (row, col): ({}, {})", i, j, self.size.0, self.size.1); 
        }
        unsafe{ &*self.ptr.as_ptr().add(i*self.offsets.0 + j*self.offsets.1) }
    }
}


// mutable element access
impl<T: PID> ops::IndexMut<(usize, usize)> for Matrix<T> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        if i>=self.size.0 || j>=self.size.1 { 
            panic!("Idx out of bounds: index: ({}, {}) while (row, col): ({}, {})", i, j, self.size.0, self.size.1); 
        }
        unsafe{ &mut *self.ptr.as_ptr().add(i*self.offsets.0 + j*self.offsets.1) }
    }
}


impl<T: PID> Clone for Matrix<T> {
    fn clone(&self) -> Self {
        let mut clone = Matrix::new( self.size.0, self.size.1 );
        for i in 0..self.size.0 {
            for j in 0..self.size.1 {
                clone[(i,j)] = self[(i,j)];
            };
        };
        clone
    }
}


impl<T: PID> cmp::PartialEq for Matrix<T> {
    fn eq(&self, rhs: &Self) -> bool {
        if self.size != rhs.size {
            return false;
        }
        self.sub(.., ..) == rhs.sub(.., ..)
    }
}

impl<T: PID + std::fmt::Debug + std::fmt::Display> fmt::Debug for Matrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\n")?;
        for i in 0..self.size.0 {
            write!(f, "| ")?;
            for j in 0..self.size.1 {
                write!(f, "{}", self[(i,j)])?;
                if j != self.size.1-1 {
                    write!(f, "\t")?;
                }
            }
            write!(f, " |\n")?;
        }
        write!(f, "")
    }
}


#[macro_export]
macro_rules! vector {
    ( $scalar_type:ty, $( $scalar:expr ), * ) => {
        {
            let mut v = Vec::new();
            $(
                v.push($scalar);
            )*
            Matrix::from_vec(v)
        }
    };
}

macro_rules! matrix_op_trait_impl {
    ($(($op_trait: ident, $op_fn: ident)) *) => {$(
        impl<T: PID> ops::$op_trait for Matrix<T>
        {
            type Output = Self;

            fn $op_fn(self, rhs: Self) -> Self::Output {
                ops::$op_trait::$op_fn(Self::sub(&self, .., ..), Self::sub(&rhs, .., ..))
            }
        }
    )*}
}
matrix_op_trait_impl!{ (Add, add) (Sub, sub) (Mul, mul) }

impl<T: PID> ops::Mul<T> for Matrix<T>
{
    type Output = Self;
    fn mul(mut self, rhs: T) -> Self::Output {
        let mut sub = self.sub_mut(.., ..);
        sub *= rhs;
        self
    }
}

impl<T: Field> ops::Div<T> for Matrix<T>
{
    type Output = Self;
    fn div(mut self, rhs: T) -> Self::Output {
        let mut sub = self.sub_mut(.., ..);
        sub /= rhs;
        self
    }
}

impl<T: PID> ops::MulAssign<T> for Matrix<T>
{
    fn mul_assign(&mut self, rhs: T) {
        let mut sub = self.sub_mut(.., ..);
        sub *= rhs;
    }
}

impl<T: Field> ops::DivAssign<T> for Matrix<T>
{
    fn div_assign(&mut self, rhs: T) {
        let mut sub = self.sub_mut(.., ..);
        sub /= rhs;
    }
}

#[allow(unused)]
impl <T: PID> Matrix<T> {
    pub fn zero(n: usize, m: usize) -> Self {
        let mut out = Matrix::new(n, m);
        for i in 0..out.size.0 {
            for j in 0..out.size.1 {
                out[(i, j)] = T::zero();
            }
        }
        out
    }

    pub fn identity(size: usize) -> Self {
        let mut mat = Matrix::zero(size, size);
        (0..size).for_each(|i| mat[(i,i)] = T::one());
        mat
    }

    pub fn smith_normal_form(mut self) -> (Self, Self, Self, Self, Self) {
        use crate::commutative::bezout_identity;

        let mut r_op = Matrix::identity(self.size.0);
        let mut r_op_inv = Matrix::identity(self.size.0);
        let mut c_op = Matrix::identity(self.size.1);
        let mut c_op_inv = Matrix::identity(self.size.1);

        let (mut i, mut j) = (0, 0);

        // // row reduce some of the trivial cases
        // let mut non_zero_rows = self.size.0;
        // for k in 0..self.size.0 {
        //     for r in k+1..self.size.0 {
        //         if self.sub(k, ..) == self.sub(r, ..) {
        //             self.row_operation(r, k, -T::one());
        //             r_op.col_operation(r, k, -T::one());
        //             r_op_inv.col_operation(r, k, T::one());
        //         }
        //     }
        // }

        // // col reduce some of the trivial cases
        // for l in 0..self.size.1 {
        //     for s in l+1..self.size.1 {
        //         if self.sub(.., l) == self.sub(.., l) {
        //             self.col_operation(s, l, -T::one());
        //             c_op.col_operation(s, l, -T::one());
        //             c_op_inv.row_operation(l, s, T::one());
        //         }
        //     }
        // }

        // // before running the gereral algorithm we run bring the pivot element to the up left
        // let mut pivot_found = true;
        // while pivot_found {
        //     pivot_found = false;
        //     for k in i..self.size.0 {
        //         let mut pivot = None;
        //         for l in (j..self.size.1).filter( |&l| self[(k,l)]!=T::zero()) {
        //             if let Some(_) = pivot {
        //                 pivot = None;
        //                 break;
        //             } else {
        //                 pivot = Some(l);
        //             }
        //         };
        //         if let Some(l) = pivot {
        //             if (0..self.size.0).filter(|&x| self[(x, l)]!=T::zero()).count() == 1 {
        //                 // Then (i,j) is a pivot.
        //                 self.swap_rows(i, k);
        //                 r_op.swap_rows(i, k);
        //                 r_op_inv.swap_cols(i, k);
        //                 self.swap_cols(j, l);
        //                 c_op.swap_cols(j, l);
        //                 c_op_inv.swap_rows(j, l);
        //                 i+=1; j+=1;
        //                 pivot_found = true;
        //                 break;
        //             }
        //         }
        //     }
        // }

        // now run the general algorithm
        while( i < self.size.0 && j < self.size.1) {
            // choosing a pivot
            j = match (j..self.size.1).find(|&l| (i..self.size.0).any(|k| self[(k, l)] != T::zero()) ){
                Some(l) => l,
                None => break,
            };
            match (j..self.size.1).find(|&l| (i..self.size.0).any(|k| self[(k, l)] != T::zero()) ) {
                Some(l) => {
                    self.swap_cols(j, l);
                    c_op.swap_cols(j, l);
                    c_op_inv.swap_rows(j, l);
                },
                None => panic!(),
            };
            if self[(i,j)] == T::zero() {
                let k = (i..self.size.0).find(|&k| self[(k,j)] != T::zero() ).unwrap();
                self.swap_rows(i, k);
                r_op.swap_rows(i, k);
                r_op_inv.swap_cols(i, k);
            }

            debug_assert!(self[(i,j)]!=T::zero());

            // Improving the pivot
            // after this while loop, the entry in (i,j) must devides all the elements below and right
            while (i+1..self.size.0).any(|k| !self[(i,j)].divides(self[(k,j)])) || (j+1..self.size.1).any(|l| !self[(i,j)].divides(self[(i,l)])) {
                // modifying 'self' to improve self[i+1.., j]
                for k in i+1..self.size.0 {
                    if !self[(i,j)].divides(self[(k,j)]) {
                        let a = self[(i,j)];
                        let b = self[(k,j)];
                        let (gcd, x, y) = bezout_identity(a, b);
                        // Now, a*x+b*y==gcd
                        let s = a.euclid_div(gcd).0;
                        let t = b.euclid_div(gcd).0;

                        let ith_row = self.sub(i, j..)*x + self.sub(k, j..)*y;
                        let kth_row = self.sub(i, j..)*(-t) + self.sub(k, j..)*s;
                        self.sub_mut(i, j..).write( ith_row );
                        self.sub_mut(k, j..).write( kth_row );

                        let ith_row = r_op.sub(i, ..)*x + r_op.sub(k, ..)*y;
                        let kth_row = r_op.sub(i, ..)*(-t) + r_op.sub(k, ..)*s;
                        r_op.sub_mut(i, ..).write( ith_row );
                        r_op.sub_mut(k, ..).write( kth_row );

                        let ith_col = r_op_inv.sub(.., i)*s + r_op_inv.sub(.., k)*t;
                        let kth_col = r_op_inv.sub(.., i)*(-y) + r_op_inv.sub(.., k)*x;
                        r_op_inv.sub_mut(.., i).write( ith_col );
                        r_op_inv.sub_mut(.., k).write( kth_col );
                        
                        debug_assert!(self[(i,j)]==gcd);
                    }
                }
                // now the "self[i, j]" must divides all the elements in "self[i+1.., j]"
                debug_assert!( (i+1..self.size.0).all(|k| self[(i,j)].divides(self[(k,j)])),
                    "'self[{i}, {j}]' must divides all the elements in 'self[{i}+1.., {j}]', but it does not: \n{:?}", self
                );

                // modifying 'self' to improve self[i, j+1..]
                for l in j+1..self.size.1 {
                    if !self[(i,j)].divides(self[(i,l)]) {
                        let a = self[(i,j)];
                        let b = self[(i,l)];
                        let (gcd, x, y) = bezout_identity(a, b);
                        let s = a.euclid_div(gcd).0;
                        let t = b.euclid_div(gcd).0;

                        let jth_row = self.sub(i.., j) * x + self.sub(i..,l) * y;
                        let lth_row = self.sub(i.., j) * (-t) + self.sub(i.., l) * s;
                        self.sub_mut(i..,j).write( jth_row );
                        self.sub_mut(i..,l).write( lth_row );
                        
                        let jth_row = c_op.sub(.., j)*x + c_op.sub(.., l)*y;
                        let lth_row = c_op.sub(.., j)*(-t) + c_op.sub(.., l)*s;
                        c_op.sub_mut(..,j).write( jth_row );
                        c_op.sub_mut(..,l).write( lth_row );
                        
                        let jth_col = c_op_inv.sub(j, ..) * s + c_op_inv.sub(l, ..)*t;
                        let lth_col = c_op_inv.sub(j, ..) * (-y) + c_op_inv.sub(l, ..) * x;
                        c_op_inv.sub_mut(j,..).write( jth_col );
                        c_op_inv.sub_mut(l,..).write( lth_col );

                        debug_assert!(self[(i,j)]==gcd, "self[(i,j)]={}, gcd={}", self[(i,j)], gcd);
                    }
                }
                // now the "self[i, j]" must divides all the elements in "self[i, j+1..]"
                debug_assert!( (j+1..self.size.1).all(|l| self[(i,j)].divides(self[(i,l)])),
                    "'self[{i}, {j}]' must divides all the elements in 'self[{i}, {j}+1]', but it does not: \n{:?}", self
                );
            }

            // Eliminating entries
            for k in i+1..self.size.0 {
                if self[(k,j)] == T::zero() { continue; }

                let ratio = self[(k,j)].euclid_div(self[(i,j)]).0;
                self.row_operation(k, i, -ratio);
                r_op.row_operation(k, i, -ratio);
                r_op_inv.col_operation(i, k, ratio);
            }

            for l in j+1..self.size.1 {
                if self[(i,l)] == T::zero() { continue; }

                let ratio = self[(i,l)].euclid_div(self[(i,j)]).0;
                self.col_operation(l, j, -ratio);
                c_op.col_operation(l, j, -ratio);
                c_op_inv.row_operation(j, l, ratio);
            }

            i += 1; j+= 1;
        };


        // Bringing all the pivots to the diagonal
        for r in 0..std::cmp::min(self.size.0, self.size.1) {
            // if the entry at (r, r) is zero then we do nothing 
            if self[(r,r)] != T::zero() { continue; }

            // otherwise, we fetch the nonzero element and swap it with the entry at (r, r)
            let s = (r+1..self.size.1).find(|&s| self[(r,s)] != T::zero() );
            let s = match s {
                Some(s) => s,
                None => break,
            };

            self[(r, r)] = self[(r,s)];
            self[(r, s)] = T::zero();
            c_op.swap_cols(r, s);
            c_op_inv.swap_rows(s, r);
        }


        // Finally, force this diagonal to satisfy a_{i, i} | a_{i+1, i+1} for each i
        let rank = match (0..std::cmp::min(self.size.0, self.size.1)).find(|&r| self[(r,r)] == T::zero()) {
            Some(r) => r,
            None => std::cmp::min(self.size.0, self.size.1),
        };

        for r in 0..rank {
            for s in r+1..rank {
                if self[(r,r)].divides(self[(s,s)]) {
                    continue;
                }

                // This col operation makes the entry at (s, r) nonzero 
                self.col_operation(r, s, T::one());
                c_op.col_operation(r, s, T::one());
                c_op_inv.row_operation(s, r, -T::one());

                // row operations
                let a = self[(r,r)];
                let b = self[(s,r)];
                let (gcd, x, y) = bezout_identity(a, b);
                let z = a.euclid_div(gcd).0;
                let w = b.euclid_div(gcd).0;

                let rth_row = self.sub(r, ..)*x + self.sub(s, ..)*y;
                let sth_row = self.sub(r, ..)*(-w) + self.sub(s, ..)*z;
                self.sub_mut(r, ..).write( rth_row );
                self.sub_mut(s, ..).write( sth_row );
                
                let rth_row = r_op.sub(r, ..)*x + r_op.sub(s, ..)*y;
                let sth_row = r_op.sub(r, ..)*(-w) + r_op.sub(s, ..)*z;
                r_op.sub_mut(r, ..).write( rth_row );
                r_op.sub_mut(s, ..).write( sth_row );

                let rth_col = r_op_inv.sub(.., r)*z + r_op_inv.sub(.., s)*w;
                let sth_col = r_op_inv.sub(.., r)*(-y) + r_op_inv.sub(.., s)*x;
                r_op_inv.sub_mut(.., r).write( rth_col );
                r_op_inv.sub_mut(.., s).write( sth_col );

                // A column operation to eliminate the entry at (r, s)
                self.col_operation(s, r, -w * y);
                c_op.col_operation(s, r, -w * y);
                c_op_inv.row_operation(s, r, w * y);

                debug_assert!(self[(r,r)]==gcd);
                debug_assert!(self[(r,s)]==T::zero());
                debug_assert!(self[(s,r)]==T::zero());
            }
        }
        debug_assert!((0..rank).all( |i| (i+1..rank).all(|k| self[(i,i)].divides(self[(k,k)]))), " Diagonal is not sorted: {:?}", self);
        // simplify r_op


        (r_op_inv, r_op, self, c_op, c_op_inv)
    }
}

#[allow(unused)]
impl <T:Field> Matrix<T> {
    pub fn rank_as_linear_map(&self) -> usize {
        let mut clone = self.clone();
        let shorter_side = std::cmp::min(self.size.0, self.size.1);
        let mut skipped = 0;
        for k in 0..shorter_side {
            if clone[(k,k)] == T::zero() {
                let result = (k+2..clone.size.1).find(|&i| clone[(k, i)] != T::zero() );
                match result {
                    Some(j) => clone.col_operation(k, j, T::one()),
                    None => (),
                }
            }
            if clone[(k,k)] == T::zero() {
                let result = (k+2..clone.size.0).find(|&j| clone[(j, k)] != T::zero() );
                match result {
                    Some(i) => clone.row_operation(k, i, T::one()),
                    None => {
                        skipped += 1;
                        continue;
                    },
                }
            }

            // Now we have the nonzero element in (k, k)
            let coeff = T::one() / clone[(k,k)];
            clone.row_operation(k, k, coeff - T::one()); // Fix

            for l in k+1..self.size.0 {
                let minus_one = T::zero() - T::one();
                let a_lk = clone[(l,k)];
                if a_lk == T::zero() { continue; } 
                clone.row_operation(l, k, minus_one * a_lk);
            }
        };

        shorter_side - skipped
    }
}


use super::FormalSum;
// Change of basis
impl<Basis> ops::Mul<Vec<Basis>> for Matrix<i128> 
    where Basis: Clone + PartialEq, i128: ops::Mul<Basis, Output = FormalSum<Basis>>
{
    type Output = Vec<FormalSum<Basis>>;

    fn mul(self, rhs: Vec<Basis>) -> Self::Output {
        assert!(self.size.1==rhs.len(), "Multiplication failed: dimensions do not match: the left has size {:?} but the right has size {}", self.size, rhs.len());
        let mut out = Vec::new();
        for i in 0..self.size.0 {
            out.push( (0..self.size.1).map(|j| self[(i, j)]*rhs[j].clone() ).sum::<FormalSum<_>>() );
        }
        out
    }
}

// numerical linear algebra
impl Matrix<f64> {
    pub fn det(mut self) -> f64 {
        if self.size.0 != self.size.1 {
            return 0_f64;
        };

        // this implementation of determinant is based of LU decomposition
        self.lu();
        (0..self.size.0)
            .map(|i| self[(i,i)])
            .product()
    }

    pub fn random_orthogonal(n: usize) -> Matrix<f64> {
        let mut out = Matrix::zero(n ,n);
        for i in 0..n {
            let mut random = Matrix::random(n-i,1);
            random /= random.two_norm();
            out.sub_mut(i..,i).write( random );
        }
        for i in (0..n-1).rev() {
            let (v, b) = out.sub(i..,i).householder_vec();
            let v = v.sub(..,..);
            let update = out.sub(i.., i+1..) - v*b * (v.transpose()*out.sub(i.., i+1..));
            out.sub_mut(i.., i+1..).write( update );
        }
        out
    }

    pub fn householder_vec(mut self) -> (Matrix<f64>, f64) {
        if self.size.1 > 1 {
            panic!("this function takes a column vector only");
        }
        
        let a: f64 = self.sub(1.., 0).transpose().dot( self.sub(1.., 0) );
        let b = self[(0,0)];
        self[(0,0)] = 1.0;

        if a==0.0 && b>=0.0 {
            (self, 0.0)
        } else if a==0.0 && b<0.0 {
            (self, -2.0)
        } else {
            let u = (a+b*b).sqrt();
            self[(0,0)] = if b<=0.0 {
                b-u
            } else {
                -a/(b+u)
            };
            let d = self[(0,0)]*self[(0,0)]; 
            self /= self[(0,0)];
            (self, 2.0*d/(a+d))
        }
    }

    pub fn solve(mut self, mut b: Matrix<f64>) -> Self {
        assert_eq!(
            self.size.0, 
            self.size.1, 
            "The first input to this function must be a square matrix,but it is not: size={:?}.",
            self.size
        );
        assert_eq!(
            self.size.1, 
            b.size.0, 
            "The number of rows of the first matrix and the number of cols of the second matrix have to match, but self.size = {:?}, b.size = {:?}.",
            self.size,
            b.size
        );
        let n = self.size.0;
        self.sub_mut(..,..).lu();
        
        // Solve Ly = b by forward substitution
        for i in 1..n {
            b[(i, 0)] = b[(i,0)] - self.sub(i, 0..i).dot( b.sub(0..i, 0) );
        }

        // Solve Ux = y by backward substitution
        b[(n-1,0)] = b[(n-1,0)] / self[(n-1,n-1)];
        for i in (0..n-1).rev() {
            b[(i, 0)] = b[(i,0)] - self.sub(i, i+1..n).dot( b.sub(i+1..n, 0));
            b[(i, 0)] /= self[(i,i)];
        }
        
        b
    }

    pub fn spectrum(mut self) -> Matrix<f64> {
        assert_eq!( self.size.0, self.size.1,
            "Input to this function must be a square matrix."
        );
        let n = self.size.0;
        self.qr();

        let mut out = vec!{0_f64; n};
        let mut i=0;
        while i < n {
            if i == n-1 || self[(i+1,i)].abs() < 0.00000000001 {
                out[i] = self[(i,i)];
                i+=1;
            } else {
                (out[i], out[i+1]) = self.sub(i..i+2, i..i+2).real_spectrum_of_two_by_two().unwrap();
                i+=2;
            }
        }

        out.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let out = {
            let mut a = Matrix::new(n,1);
            (0..n).for_each( |i| a[(i,0)]=out[i] );
            a
        };
        out
    }

    pub fn spectrum_with_invariant_space(self) -> (Matrix<f64>, Matrix<f64>) {
        assert_eq!( self.size.0, self.size.1,
            "Input to this function must be a square matrix."
        );
        let n = self.size.0;

        // initialize the eigenvectors
        let mut eigen_vectors = Matrix::<f64>::new(n, n);
        let elem = 1.0 / (n as f64).sqrt();
        for i in 0..n {
            for j in 0..n {
                eigen_vectors[(i,j)] = elem;
            }
        }

        // get the eigenvalues
        let spectrum = self.clone().spectrum();

        // run the inverse iteration to get the eigenvector for each of the eigenvalues
        for i in 0..n {
            let lambda = spectrum[(i,0)];

            let mut m = self.clone();
            for j in 0..n {
                m[(j,j)] -= lambda;
            }

            let mut update = m.clone().solve( eigen_vectors.sub(.., i).as_matrix() );
            update /= update.two_norm();
            eigen_vectors.sub_mut(.., i).write( update );
            let mut update = m.solve( eigen_vectors.sub(.., i).as_matrix() );
            update /= update.two_norm();
            eigen_vectors.sub_mut(.., i).write( update );
        }

        (spectrum, eigen_vectors)
    }
}

// impl<'a> Matrix<f64> {
//     pub fn one_norm(&self) -> f64 {
//         Matrix::<f64>::sub(self, ..,..).one_norm()
//     }
// }

impl<'a> SubMatrix<'a, f64> {
    pub fn det(self) -> f64 {
        // this implementation of determinant is based on the LU decomposition
        self.as_matrix().det()
    }

    pub fn one_norm(self) -> f64 {
        (0..self.size.1)
            .map(|j| (0..self.size.0).map(|i| self[(i ,j)].abs() ).sum() )
            .fold(0.0, |max, val| if max<=val { val } else { max }  )
    }

    pub fn two_norm(self) -> f64 {
        if self.size.1 ==  1 {
            ((0..self.size.0).map(|i| self[(i,0)]*self[(i,0)] ).sum::<f64>()).sqrt()
        } else {
            panic!("2-norm for matrices are not supported yet.");
        }
    }

    pub fn inf_norm(self) -> f64 {
        (0..self.size.0)
            .map(|i| (0..self.size.1).map(|j| self[(i ,j)].abs() ).sum() )
            .fold(0.0, |max, val| if max<=val { val } else { max } )
    }

    pub fn frobenius_norm(self) -> f64 {
        (0..self.size.0)
            .map( |i|
                (0..self.size.1).map(|j| self[(i,j)]*self[(i,j)] ).sum::<f64>()
            )
            .sum::<f64>()
            .sqrt()
    }

    pub fn householder_vec(self) -> (Matrix<f64>, f64) {
        self.as_matrix().householder_vec()
    }

    fn real_spectrum_of_two_by_two(self) -> Option<(f64, f64)> {
        debug_assert_eq!( self.size, (2,2) );
        let det = self[(0,0)] * self[(1,1)] - self[(0,1)] * self[(1,0)];
        let m = (self[(0,0)] + self[(1,1)]) / 2.0;
        if m*m-det > 0.0 {
            let rad = (m*m-det).sqrt();
            Some((m-rad, m+rad))
        } else {
            None
        }
    }
}

impl<'a> SubMatrixMut<'a, f64> {
    pub fn lu(&mut self) {
        assert_eq!(self.size.0, self.size.1, "Input has to this function has to be a square matrix.");

        let n = self.size.0;
        for i in 0..n-1 {
            let update = self.sub(i+1..n, i) / self[(i,i)];
            self.sub_mut(i+1..n,i).write( update );
            let update = self.sub(i+1..n,i+1..n) - self.sub(i+1..n,i) * self.sub(i, i+1..n);
            self.sub_mut(i+1..n, i+1..n).write( update );
        }
    }


    pub fn hessenberg_form(&mut self) {
        if self.size.0 != self.size.1 {
            panic!("Input to this function must be a square matrix.");
        }

        let n = self.size.0; // = self.size.1
        for k in 0..n-2 {
            let (v, b) = self.sub(k+1..n, k).householder_vec();
            let v = v.sub(.., ..);
            let v_t = v.transpose();

            let update = self.sub(k+1..n, k..n) - v*(v_t * self.sub(k+1..n, k..n))*b;
            self.sub_mut(k+1..n, k..n).write( update );

            let update = self.sub(0..n, k+1..n) - (self.sub(0..n, k+1..n)*v)*v_t*b;
            self.sub_mut(0..n, k+1..n).write( update );
        };
    }

    pub fn raw_francis_step(&mut self, return_op: bool) -> Matrix<f64> {
        assert_eq!( self.size.0, self.size.1,
            "Input to this function must be a square matrix."
        );
        let n = self.size.0;

        let op = if return_op {
            Matrix::identity(n)
        } else {
            Matrix::identity(1)
        };

        let s = self[(n-2,n-2)] + self[(n-1,n-1)];
        let t = self[(n-2,n-2)]*self[(n-1,n-1)] - self[(n-2,n-1)]*self[(n-1,n-2)];
        let mut x = self[(0,0)]*self[(0,0)] + self[(0,1)]*self[(1,0)] -s*self[(0,0)] + t;
        let mut y = self[(1,0)] * (self[(0,0)] + self[(1,1)] - s);
        let mut z = self[(1,0)]*self[(2,1)];

        for k in 0..n-2 {
            let (v, b) = crate::matrix!{f64; [[x,y,z]]}.transpose().householder_vec();
            let v = v.sub(.., ..);
            let q = if k>0 { k-1 } else { 0 };
            let update = self.sub(k..k+3, q..n) -
                v * (v.transpose() * self.sub(k..k+3, q..n) * b);
            self.sub_mut(k..k+3, q..n).write( 
                update
            );

            let r = if k<n-3 { k+4 } else { n };
            let update = self.sub(0..r, k..k+3) -
                (self.sub(0..r, k..k+3) * v) * (v.transpose()*b);
            self.sub_mut(0..r, k..k+3).write(
                update
            );
            x = self[(k+1, k)];
            y = self[(k+2, k)];
            if k < n-3 {
                z = self[(k+3, k)];
            }
        }

        let (v, b) = crate::matrix!{f64; [[x,y]]}.transpose().householder_vec();
        let v = v.sub(.., ..);
        let update = self.sub(n-2..n, n-3..n) -
            v * (v.transpose() * self.sub(n-2..n, n-3..n) * b);
        self.sub_mut(n-2..n, n-3..n).write( update );

        let update = self.sub(0..n, n-2..n) -
            (self.sub(0..n, n-2..n) * v) * (v.transpose()*b);
        self.sub_mut(0..n, n-2..n).write( update );

        op
    }

    pub fn francis_step(&mut self) {
        self.raw_francis_step(false);
    }

    pub fn francis_step_with_op(&mut self) -> Matrix<f64> {
        self.raw_francis_step(true)
    }

    fn raw_qr(&mut self, returns_op: bool) -> Matrix<f64> {
        assert_eq!( self.size.0, self.size.1,
            "Input to this function must be a square matrix."
        );
        let n = self.size.0;

        self.hessenberg_form();
        let mut op = if returns_op {
            Matrix::identity(n)
        } else {
            Matrix::identity(1)
        };

        let tolerance = 0.000000000000001;
        let mut q=0;
        let mut p;

        while q!=n {
            for i in 1..n {
                if self[(i,i-1)].abs() < tolerance * (self[(i-1,i-1)].abs() + self[(i,i)].abs()) {
                    self[(i,i-1)] = 0.0;
                }
            }

            q = (0..n-2)
                .find( |&i| self[(n-i-1,n-i-2)].abs()>tolerance && self[(n-i-2,n-i-3)].abs()>tolerance )
                .unwrap_or(n);

            p = (1..n-q)
                .rev()
                .find(|&i| self[(i,i-1)].abs()<tolerance )
                .unwrap_or(0);

            if q<n {
                debug_assert!( (n-q)-p > 2);
                if returns_op {
                    op = op * self.sub_mut(p..n-q, p..n-q).francis_step_with_op();
                } else {
                    self.sub_mut(p..n-q, p..n-q).francis_step();
                };
            }
        };

        op
    }

    pub fn qr(&mut self) {
        self.raw_qr(false);
    }

    pub fn qr_with_op(&mut self) -> Matrix<f64> {
        self.raw_qr(true)
    }
}

impl std::fmt::Display for Matrix<f64> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision().unwrap_or(3);
        write!(f, "\n")?;
        for i in 0..self.size.0 {
            write!(f, "| ")?;
            for j in 0..self.size.1 {
                write!(f, "{:.decimals$}", self[(i,j)])?;
                if j != self.size.1-1 {
                    write!(f, "\t")?;
                }
            }
            write!(f, " |\n")?;
        }
        write!(f, "")
    }
}

macro_rules! immut_ref_fn_cast_impl_for_matrix {
    ( $(($fn: ident, $out: ty)) * ) => {
        impl<'a> Matrix<f64> {
            $(
                pub fn $fn(&self) -> f64 {
                    self.sub(.., ..).$fn()
                }
            ) *
        }
    }
}
immut_ref_fn_cast_impl_for_matrix!{ 
    (one_norm, ())
    (two_norm, ())
    (inf_norm, ())
    (frobenius_norm, ())
}

macro_rules! mut_ref_fn_cast_impl_for_matrix {
    ( $(($fn: ident, $out: ty)) * ) => {
        impl Matrix<f64> {
            $(
                pub fn $fn(&mut self) -> $out {
                    self.sub_mut(.., ..).$fn();
                }
            ) *
        }
    }
}
mut_ref_fn_cast_impl_for_matrix!{ 
    (lu, ())
    (hessenberg_form, ())
    (francis_step, ())
    (qr, ()) 
    (francis_step_with_op, ())
    (qr_with_op, ()) 
}

// macro_rules! fn_cast_impl_for_matrix {
//     ( $(($fn: ident, $out: ty)) * ) => {
//         impl Matrix<f64> {
//             $(
//                 pub fn $fn(mut self) -> $out {
//                     self.sub_mut(.., ..).$fn()
//                 }
//             ) *
//         }
//     }
// }

// fn_cast_impl_for_matrix!{ 
// }
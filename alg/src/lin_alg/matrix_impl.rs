use crate::commutative::{PID, Field};
use crate::matrix;

use std::fmt::{Display, Debug};
use std::{
    alloc,
    marker::PhantomData,
    mem,
    ptr,
    fmt,
    cmp,
    ops
};

use super::Matrix;
use super::SubMatrix;

#[inline]
pub const fn metadata<T: ?Sized>(ptr: *const T) -> usize {
    unsafe { PtrRepr { const_ptr: ptr }.components.metadata }
}

#[inline]
pub const fn compress(size: (usize, usize), offsets: (usize, usize)) -> usize {
    let data_offset = mem::size_of::<usize>() * 8 / 4;
    size.0 << data_offset*3
        | size.1 << data_offset*2
        | offsets.0 << data_offset*1
        | offsets.1
}

#[inline]
pub const fn extract(data: usize) -> ((usize, usize), (usize, usize)) {
    let data_offset = mem::size_of::<usize>() * 8 / 4;
    let window: usize = !(!0 << data_offset);
    let size = ((data>>data_offset*3) & window, (data>>data_offset*2) & window);
    let offsets = ( (data>>data_offset*1) & window, data & window);
    (size, offsets)
}

#[inline]
pub const fn submatrix_from_raw_parts<T: ?Sized>(
    data_address: *const (),
    size: (usize, usize),
    offsets: (usize, usize),
) -> *const T {
    let metadata = compress(size, offsets);
    unsafe { PtrRepr { components: PtrComponents { data_address, metadata } }.const_ptr }
}

#[inline]
pub const fn submatrix_from_raw_parts_mut<T: ?Sized>(
    data_address: *mut (),
    size: (usize, usize),
    offsets: (usize, usize),
) -> *mut T {
    let metadata = compress(size, offsets);
    unsafe { PtrRepr { components: PtrComponents { data_address, metadata } }.mut_ptr }
}

#[inline]
pub fn slice_from_sub_matrix<T: PID>(m: &SubMatrix<T>) -> *const [T] {
    let (size, offsets) = extract( metadata(m));
    let metadata = (size.0 * offsets.0 - 1) + (size.1 * offsets.1 - 1) + 1;
    let data_address = unsafe{ PtrRepr{const_ptr: m as *const SubMatrix<T>}.components.data_address};
    unsafe { PtrRepr { components: PtrComponents { data_address, metadata } }.const_ptr }
}

#[inline]
pub fn slice_from_sub_matrix_mut<T: PID>(m: &mut SubMatrix<T>) -> *mut [T] {
    let (size, offsets) = extract( metadata(m));
    let metadata = (size.0 * offsets.0 - 1) + (size.1 * offsets.1 - 1) + 1;
    let data_address = unsafe{ PtrRepr{ mut_ptr: m as *mut SubMatrix<T> }.components.data_address };
    unsafe { PtrRepr { components: PtrComponents { data_address, metadata } }.mut_ptr }
}

#[repr(C)]
union PtrRepr<T: ?Sized> {
    const_ptr: *const T,
    mut_ptr: *mut T,
    components: PtrComponents,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct PtrComponents {
    data_address: *const (),
    metadata: usize,
}

impl<T: PID> SubMatrix<T> {
    #[inline]
    pub fn size(&self) -> (usize, usize) {
        extract( metadata(self) ).0
    }

    #[inline]
    pub fn size_as_range(&self) -> (ops::Range<usize>, ops::Range<usize>) {
        (0..self.size().0, 0..self.size().1)
    }

    #[inline]
    fn offsets(&self) -> (usize, usize) {
        extract( metadata(self) ).1
    }

    #[inline]
    fn as_ptr(&self) -> *const T {
        unsafe{ (&*slice_from_sub_matrix(self)).as_ptr() }
    }

    #[inline]
    fn as_mut_ptr(&mut self) -> *mut T {
        unsafe{ (&mut *slice_from_sub_matrix_mut(self)).as_mut_ptr() }
    }

    #[inline]
    pub fn transpose(&self) -> &Self {
        let (size, offsets) = extract(metadata(self));
        let new_size = (size.1, size.0);
        let new_offsets = (offsets.1, offsets.0);
        
        unsafe{ &*submatrix_from_raw_parts((self as *const Self).cast(), new_size, new_offsets) }
    }

    pub fn dot(&self, rhs: &Self) -> T {
        // This function makes sure that 'self' is a row vector and 'rhs' is a col vector
        assert!(self.size().0 == 1 && rhs.size().1 == 1);

        // This function makes sure that the arguments are of the same size
        assert!(self.size().1 == rhs.size().0, 
            "The dot product cannot be performed between the matrices of different size. left has {:?} cols, but the right has {:?} rows",
            self.size().1,
            rhs.size().0
        );

        (0..self.size().1).map(|i| self[(0, i)].clone() * rhs[(i,0)].clone()).sum()
    }

    pub fn as_matrix(&self) -> Matrix<T> {
        let mut out = Matrix::new(self.size().0, self.size().1);
        for i in 0..self.size().0 {
            for j in 0..self.size().1 {
                out.write( (i,j), self[(i, j)].clone() );
            }
        }

        out
    }

    pub fn swap_rows(&mut self, r1: usize, r2: usize) {
        for j in 0..self.size().1 {
            let tmp = self[(r1, j)].clone();
            self[(r1, j)] = self[(r2, j)].clone();
            self[(r2, j)] = tmp;
        }
    }

    pub fn swap_cols(&mut self, c1: usize, c2: usize) {
        for i in 0..self.size().0 {
            let tmp = self[(i, c1)].clone();
            self[(i, c1)] = self[(i, c2)].clone();
            self[(i, c2)] = tmp;
        }
    }

    pub fn row_operation(&mut self, r1: usize, r2: usize, scalar: T) {
        let r2 = &self[(r2, ..)] * scalar;
        (0..self.size().1).for_each(|i| self[(r1, i)] += r2[(0, i)].clone());
    }

    pub fn col_operation(&mut self, c1: usize, c2: usize, scalar: T) {
        let c2 = &self[(.., c2)] * scalar;
        (0..self.size().0).for_each(|j| self[(j, c1)] += c2[(j, 0)].clone());
    }

    pub fn write(&mut self, m: &SubMatrix<T>) {
        assert_eq!(self.size(), m.size(), "cannot write a matrix to a (sub)matrix of different size");
        for i in 0..self.size().0 {
            for j in 0..self.size().1 {
                self[(i,j)] = m[(i,j)].clone();
            }
        }
    }

    pub fn write_and_move(&mut self, m: Matrix<T>) {
        self.write(&m);
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

use ops::Bound;
fn unpack_bounds_for_sub(range: &impl ModifiedUsizeRangeBounds, size: usize) -> (usize, usize) {
    let start = match range.modified_start_bound() {
        Bound::Included(&x) => x,
        Bound::Excluded(_) => panic!("input range not valid"),
        Bound::Unbounded => 0,
    };

    let end = match range.modified_end_bound() {
        Bound::Included(&x) => x+1,
        Bound::Excluded(&x) => x,
        Bound::Unbounded => size,
    };

    (start, end)
}

fn check_bounds_for_sub<const ROW: bool>(start: usize, end: usize, size: usize) {
    assert!(end - start <= size, "range not valid: specified range is {start}..{end} but this matrix has {size} {}s.", 
        if ROW {"row"} else {"column"}
    );
    assert!(end <= size, "range not valid: specified range is {start}..{end} but this matrix has {size} {}s.",
        if ROW {"row"} else {"column"}
    );
}

macro_rules! submatrix_index_by_ranges_impl {
    ($(($range1: ty, $range2: ty)) *) => {$(
        impl<T: PID> ops::Index<($range1, $range2)> for SubMatrix<T> {
            type Output = SubMatrix<T>;
            fn index(&self, (rows, cols): ($range1, $range2)) -> &Self::Output {
                let (row_start, row_end) = unpack_bounds_for_sub(&rows, self.size().0);
                let (col_start, col_end) = unpack_bounds_for_sub(&cols, self.size().1);

                check_bounds_for_sub::<true>(row_start, row_end, self.size().0);
                check_bounds_for_sub::<false>(col_start, col_end, self.size().1);

                let start_offset = self.offsets().0 * row_start + self.offsets().1 * col_start;

                let size = (row_end - row_start, col_end - col_start);

                unsafe {
                    &*submatrix_from_raw_parts( self.as_ptr().add(start_offset).cast(), size, self.offsets() ) 
                }
            }
        }

        impl<T: PID> ops::IndexMut<($range1, $range2)> for SubMatrix<T> {
            fn index_mut(&mut self, (rows, cols): ($range1, $range2)) -> &mut Self::Output {
                let (row_start, row_end) = unpack_bounds_for_sub(&rows, self.size().0);
                let (col_start, col_end) = unpack_bounds_for_sub(&cols, self.size().1);

                check_bounds_for_sub::<true>(row_start, row_end, self.size().0);
                check_bounds_for_sub::<false>(col_start, col_end, self.size().1);

                let start_offset = self.offsets().0 * row_start + self.offsets().1 * col_start;

                let size = (row_end - row_start, col_end - col_start);

                unsafe {
                    &mut *submatrix_from_raw_parts_mut( self.as_mut_ptr().add(start_offset).cast(), size, self.offsets() ) 
                }
            }
        }
    ) *}
}

use std::ops::*;
submatrix_index_by_ranges_impl!{
    (Range<usize>, Range<usize>) (Range<usize>, RangeFrom<usize>) (Range<usize>, RangeFull) (Range<usize>, RangeInclusive<usize>) (Range<usize>, RangeTo<usize>) (Range<usize>, RangeToInclusive<usize>) (Range<usize>, usize)
    (RangeFrom<usize>, Range<usize>) (RangeFrom<usize>, RangeFrom<usize>) (RangeFrom<usize>, RangeFull) (RangeFrom<usize>, RangeInclusive<usize>) (RangeFrom<usize>, RangeTo<usize>) (RangeFrom<usize>, RangeToInclusive<usize>) (RangeFrom<usize>, usize)
    (RangeFull, Range<usize>) (RangeFull, RangeFrom<usize>) (RangeFull, RangeFull) (RangeFull, RangeInclusive<usize>) (RangeFull, RangeTo<usize>) (RangeFull, RangeToInclusive<usize>) (RangeFull, usize)
    (RangeInclusive<usize>, Range<usize>) (RangeInclusive<usize>, RangeFrom<usize>) (RangeInclusive<usize>, RangeFull) (RangeInclusive<usize>, RangeInclusive<usize>) (RangeInclusive<usize>, RangeTo<usize>) (RangeInclusive<usize>, RangeToInclusive<usize>) (RangeInclusive<usize>, usize)
    (RangeTo<usize>, Range<usize>) (RangeTo<usize>, RangeFrom<usize>) (RangeTo<usize>, RangeFull) (RangeTo<usize>, RangeInclusive<usize>) (RangeTo<usize>, RangeTo<usize>) (RangeTo<usize>, RangeToInclusive<usize>) (RangeTo<usize>, usize)
    (RangeToInclusive<usize>, Range<usize>) (RangeToInclusive<usize>, RangeFrom<usize>) (RangeToInclusive<usize>, RangeFull) (RangeToInclusive<usize>, RangeInclusive<usize>) (RangeToInclusive<usize>, RangeTo<usize>) (RangeToInclusive<usize>, RangeToInclusive<usize>) (RangeToInclusive<usize>, usize)
    (usize, Range<usize>) (usize, RangeFrom<usize>) (usize, RangeFull) (usize, RangeInclusive<usize>) (usize, RangeTo<usize>) (usize, RangeToInclusive<usize>)
}

#[inline]
fn check_bounds_for_index<T: PID>(r: usize, c: usize, m: &SubMatrix<T> ) {
    let (row_range, col_range) = m.size_as_range();
    assert!( 
        row_range.contains( &r ) && col_range.contains( &c ),
        "Index out of bounds: index=({r}, {c}), but the matrix is of size {:?}", m.size()
    );
}

impl<T: PID> ops::Index<(usize, usize)> for SubMatrix<T> {
    type Output = T;
    fn index(&self, (r, c): (usize, usize)) -> &T {
        check_bounds_for_index(r, c, self);
        unsafe{ &*self.as_ptr().add( r*self.offsets().0 + c*self.offsets().1 ) }
    }
}

impl<T: PID> ops::IndexMut<(usize, usize)> for SubMatrix<T> {
    fn index_mut(&mut self, (r, c): (usize, usize)) -> &mut T {
        check_bounds_for_index(r, c, self);
        unsafe{ &mut *self.as_mut_ptr().add( r*self.offsets().0 + c*self.offsets().1 )}
    }
}

macro_rules! submatrix_add_sub_impl {
    ($(($bin_op_tr: ident, $bin_op_fn: ident)) *) => {$(
        impl<T: PID> ops::$bin_op_tr for &SubMatrix<T> {
            type Output = Matrix<T>;
            fn $bin_op_fn(self, rhs: Self) -> Self::Output {
                if self.size() != rhs.size() { panic!(); }

                let mut out = Matrix::new(self.size().0, self.size().1);
                for i in 0..self.size().0 {
                    for j in 0..self.size().1 {
                        out[(i,j)] = $bin_op_tr::$bin_op_fn(self[(i,j)].clone(), rhs[(i,j)].clone());
                    }
                }
                out
            }
        }
    )*}
}
submatrix_add_sub_impl!((Add, add) (Sub, sub));


macro_rules! submatrix_add_sub_assign_impl {
    ($(($bin_op_tr: ident, $bin_op_fn: ident)) *) => {$(
        impl<T: PID> ops::$bin_op_tr<&SubMatrix<T>> for Matrix<T> {
            fn $bin_op_fn(&mut self, rhs: &SubMatrix<T>) {
                assert_eq!( self.size(), rhs.size(), 
                    "Addition or subtraction cannnot be performed between matrices of different size." 
                );

                for i in 0..self.size().0 {
                    for j in 0..self.size().1 {
                        $bin_op_tr::$bin_op_fn(&mut self[(i,j)], rhs[(i,j)].clone());
                    }
                };
            }
        }

        impl<T: PID> ops::$bin_op_tr<Matrix<T>> for Matrix<T> {
            fn $bin_op_fn(&mut self, rhs: Matrix<T>) {
                assert_eq!( self.size(), rhs.size(), 
                    "Addition or subtraction cannnot be performed between matrices of different size." 
                );

                for i in 0..self.size().0 {
                    for j in 0..self.size().1 {
                        $bin_op_tr::$bin_op_fn(&mut self[(i,j)], rhs[(i,j)].clone());
                    }
                };
            }
        }
    )*}
}
submatrix_add_sub_assign_impl!((AddAssign, add_assign) (SubAssign, sub_assign));



impl<T: PID> ops::Mul for &SubMatrix<T> {
    type Output = Matrix<T>;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut prod = Matrix::new(self.size().0, rhs.size().1);
        for i in 0..prod.size().0 {
            for j in 0..prod.size().1 {
                prod[(i, j)] = self[(i, ..)].dot(&rhs[(.., j)]);
            }
        };
        prod
    }
}



/// multiplication between
///     (1) &SubMatrix<T> and &SubMatrix<T>
///     (2) Matrix<T> and &SubMatrix<T>
///     (3) &SubMatrix<T> and Matrix<T>
///     (4) Matrix<T> and Matrix<T>
///  are suppoerted

impl<T: PID> ops::Mul<&SubMatrix<T>> for Matrix<T> {
    type Output = Matrix<T>;
    fn mul(self, rhs: &SubMatrix<T>) -> Matrix<T> {
        let lhs = &*self;
        lhs*rhs
    }
}

impl<T: PID> ops::Mul<Matrix<T>> for &SubMatrix<T> {
    type Output = Matrix<T>;
    fn mul(self, rhs: Matrix<T>) -> Matrix<T> {
        let rhs = &*rhs;
        self*rhs
    }
}

impl<T: PID> ops::Mul<Matrix<T>> for Matrix<T> {
    type Output = Matrix<T>;
    fn mul(self, rhs: Matrix<T>) -> Matrix<T> {
        let lhs = &*self;
        let rhs = &*rhs;
        lhs*rhs
    }
}

/// Scalar multiplication, addition/subtraction of matrices consumes the ownership of Matrix<T>
macro_rules! binary_op_impl {
    ($(($bin_op_tr: ident, $bin_op_fn: ident, $op_assign: tt)) *) => {$(
        impl<T: PID> ops::$bin_op_tr<&SubMatrix<T>> for Matrix<T> {
            type Output = Matrix<T>;
            fn $bin_op_fn(mut self, rhs: &SubMatrix<T>) -> Matrix<T> {
                self $op_assign rhs;
                self
            }
        }

        impl<T: PID> ops::$bin_op_tr<Matrix<T>> for &SubMatrix<T> {
            type Output = Matrix<T>;
            fn $bin_op_fn(self, rhs: Matrix<T>) -> Matrix<T> {
                let mut lhs =self.as_matrix();
                lhs $op_assign rhs;
                lhs
            }
        }

        impl<T: PID> ops::$bin_op_tr<Matrix<T>> for Matrix<T> {
            type Output = Matrix<T>;
            fn $bin_op_fn(mut self, rhs: Matrix<T>) -> Matrix<T> {
                self $op_assign &rhs[(..,..)];
                self
            }
        }
    )*}
}

binary_op_impl!{ 
    (Add, add, +=)
    (Sub, sub, -=)
}

impl<T: PID> ops::Mul<T> for Matrix<T> {
    type Output = Matrix<T>;
    fn mul(mut self, rhs: T) -> Self::Output {
        *self *= rhs;
        self
    }
}

impl<T: Field> ops::Div<T> for Matrix<T> {
    type Output = Matrix<T>;
    fn div(mut self, rhs: T) -> Self::Output {
        *self /= rhs;
        self
    }
}

// scalar multiplication
impl<T: PID> ops::Mul<T> for &SubMatrix<T> {
    type Output = Matrix<T>;
    fn mul(self, rhs: T) -> Self::Output {
        let mut out = Matrix::new(self.size().0, self.size().1);
        for i in 0..out.size().0 {
            for j in 0..out.size().1 {
                out[(i,j)] = self[(i,j)].clone() * rhs.clone();
            }
        };
        out
    }
}

macro_rules! scalar_mul_impl {
    ($($pid: ty) *) => ($(
        impl ops::Mul<&SubMatrix<$pid>> for $pid {
            type Output = Matrix<$pid>;
            fn mul(self, rhs: &SubMatrix<$pid>) -> Self::Output {
                ops::Mul::mul(rhs, self)
            }
        }
    )*)
}

scalar_mul_impl!{ isize i8 i16 i32 i64 i128 crate::commutative::Rational f64 }

impl<T: Field> ops::Div<T> for &SubMatrix<T> {
    type Output = Matrix<T>;
    fn div(self, rhs: T) -> Self::Output {
        let mut out = Matrix::new(self.size().0, self.size().1);
        for i in 0..out.size().0 {
            for j in 0..out.size().1 {
                out[(i,j)] = self[(i,j)].clone() / rhs.clone();
            }
        };
        out
    }
}


impl<T: PID> ops::MulAssign<T> for SubMatrix<T> {
    fn mul_assign(&mut self, rhs: T) {
        for i in 0..self.size().0 {
            for j in 0..self.size().1{
                self[(i,j)] *= rhs.clone();
            }
        }
    }
}

impl<T: PID> ops::MulAssign<T> for Matrix<T> {
    fn mul_assign(&mut self, rhs: T) {
        **self *= rhs;
    }
}

impl<T: Field> ops::DivAssign<T> for SubMatrix<T> {
    fn div_assign(&mut self, rhs: T) {
        for i in 0..self.size().0 {
            for j in 0..self.size().1{
                self[(i,j)] /= rhs.clone();
            }
        };
    }
}

impl<T: Field> ops::DivAssign<T> for Matrix<T> {
    fn div_assign(&mut self, rhs: T) {
        **self /= rhs;
    }
}

impl<T: PID> cmp::PartialEq for SubMatrix<T> {
    fn eq(&self, rhs: &Self) -> bool {
        if self.size() != rhs.size() {
            return false;
        }
        for i in 0..self.size().0 {
            for j in 0..self.size().1 {
                if self[(i, j)] != rhs[(i, j)] {
                    return false;
                }
            }
        }
        true
    }
}

impl<T: PID> cmp::PartialEq<Matrix<T>> for SubMatrix<T> {
    fn eq(&self, rhs: &Matrix<T>) -> bool {
        self == rhs
    }
}

impl<T: PID> cmp::PartialEq<Matrix<T>> for Matrix<T> {
    fn eq(&self, rhs: &Matrix<T>) -> bool {
        **self == **rhs
    }
}

#[allow(unused)]
impl<T: PID> Matrix<T> {
    #[inline]
    fn size(&self) -> (usize, usize) {
        self.size
    }

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
        let layout = alloc::Layout::array::<T>(self.size().0 * self.size().1).unwrap();
        
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
        for i in 0..self.size().0 {
            for j in 0..self.size().1 {
                if self[(i,j)] != T::zero() && i!=j {
                    return false;
                };
            };
        };
        true
    }
    
    pub fn swap_rows(&mut self, r1: usize, r2: usize) {
        self[(.., ..)].swap_rows(r1, r2);
    }

    pub fn swap_cols(&mut self, c1: usize, c2: usize) {
        self[(.., ..)].swap_cols(c1, c2);
    }

    pub fn row_operation(&mut self, r1: usize, r2: usize, scalar: T) {
        self[(.., ..)].row_operation(r1, r2, scalar);
    }

    pub fn col_operation(&mut self, c1: usize, c2: usize, scalar: T) {
        self[(.., ..)].col_operation(c1, c2, scalar);
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
                out[(i, j)] = v[i][j].clone();
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

    fn write(&mut self, (i,j): (usize, usize), val: T) {
        unsafe{ std::ptr::write( self.ptr().add(i*self.offsets().0 + j*self.offsets().1 ), val); }
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
        for i in 0..out.size().0 {
            for j in 0..out.size().1 {
                out[(i,j)] = rng.gen();
            }
        }
        out
    }

    pub fn random_symmetric(n: usize, m: usize) -> Matrix<T> {
        let mut out = Matrix::new(n ,m);
        let mut rng = rand::thread_rng();
        for i in 0..out.size().0 {
            for j in i..out.size().1 {
                out[(i,j)] = rng.gen();
                out[(j,i)] = out[(i,j)].clone();
            }
        }
        out
    }
}

impl<T: PID> ops::Deref for Matrix<T> {
    type Target = SubMatrix<T>;
    fn deref(&self) -> &Self::Target {
        unsafe{ 
            &*submatrix_from_raw_parts(self.ptr.as_ptr().cast(), self.size, self.offsets)
        }
    }
}

impl<T: PID> ops::DerefMut for Matrix<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        unsafe{ 
            &mut *submatrix_from_raw_parts_mut(self.ptr.as_ptr().cast(), self.size, self.offsets)
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
                    alloc::Layout::array::<T>(self.size().0 * self.size().1).unwrap(),
                );
            }
        }
    }
}

impl<T: PID> Clone for Matrix<T> {
    fn clone(&self) -> Self {
        (*self).as_matrix()
    }
}

impl<T: PID + std::fmt::Debug> fmt::Debug for SubMatrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\n")?;
        for i in 0..self.size().0 {
            write!(f, "| ")?;
            for j in 0..self.size().1 {
                write!(f, "{:?}", self[(i,j)])?;
                if j != self.size().1-1 {
                    write!(f, "\t")?;
                }
            }
            write!(f, " |\n")?;
        }
        write!(f, "")
    }
}

impl<T: PID + std::fmt::Display> fmt::Display for Matrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&(**self), f)
    }
}

impl<T: PID + std::fmt::Debug> fmt::Debug for Matrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(&(**self), f)
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

impl <T: PID> Matrix<T> {
    pub fn zero(n: usize, m: usize) -> Self {
        assert!(n>0 && m>0, "the size of the matrix must be a pair positive integers but size=({n}, {m}).");
        let mut out = Matrix::new(n, m);
        for i in 0..n {
            for j in 0..m {
                out.write((i,j), T::zero());
            }
        }
        out
    }

    pub fn identity(size: usize) -> Self {
        assert!(size>0 , "the size of the identity matrix must be a positive integer but size={size}).");
        let mut mat = Matrix::zero(size, size);
        (0..size).for_each(|i| mat[(i,i)] = T::one());
        mat
    }
}

impl <T: PID+Debug+Display+Copy> Matrix<T> {
    pub fn smith_normal_form(mut self) -> (Self, Self, Self, Self, Self) {
        use crate::commutative::bezout_identity;

        let mut r_op = Matrix::identity(self.size().0);
        let mut r_op_inv = Matrix::identity(self.size().0);
        let mut c_op = Matrix::identity(self.size().1);
        let mut c_op_inv = Matrix::identity(self.size().1);

        let (mut i, mut j) = (0, 0);

        // // row reduce some of the trivial cases
        // let mut non_zero_rows = self.size().0;
        // for k in 0..self.size().0 {
        //     for r in k+1..self.size().0 {
        //         if self[(k, ..) == self[(r, ..) {
        //             self.row_operation(r, k, -T::one());
        //             r_op.col_operation(r, k, -T::one());
        //             r_op_inv.col_operation(r, k, T::one());
        //         }
        //     }
        // }

        // // col reduce some of the trivial cases
        // for l in 0..self.size().1 {
        //     for s in l+1..self.size().1 {
        //         if self[(.., l) == self[(.., l) {
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
        //     for k in i..self.size().0 {
        //         let mut pivot = None;
        //         for l in (j..self.size().1).filter( |&l| self[(k,l)]!=T::zero()) {
        //             if let Some(_) = pivot {
        //                 pivot = None;
        //                 break;
        //             } else {
        //                 pivot = Some(l);
        //             }
        //         };
        //         if let Some(l) = pivot {
        //             if (0..self.size().0).filter(|&x| self[(x, l)]!=T::zero()).count() == 1 {
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
        while i < self.size().0 && j < self.size().1 {
            // choosing a pivot
            j = match (j..self.size().1).find(|&l| (i..self.size().0).any(|k| self[(k, l)] != T::zero()) ){
                Some(l) => l,
                None => break,
            };
            match (j..self.size().1).find(|&l| (i..self.size().0).any(|k| self[(k, l)] != T::zero()) ) {
                Some(l) => {
                    self.swap_cols(j, l);
                    c_op.swap_cols(j, l);
                    c_op_inv.swap_rows(j, l);
                },
                None => panic!(),
            };
            if self[(i,j)] == T::zero() {
                let k = (i..self.size().0).find(|&k| self[(k,j)] != T::zero() ).unwrap();
                self.swap_rows(i, k);
                r_op.swap_rows(i, k);
                r_op_inv.swap_cols(i, k);
            }

            debug_assert!(self[(i,j)]!=T::zero());

            // Improving the pivot
            // after this while loop, the entry in (i,j) must devides all the elements below and right
            while (i+1..self.size().0).any(|k| !self[(i,j)].divides(self[(k,j)])) || (j+1..self.size().1).any(|l| !self[(i,j)].divides(self[(i,l)])) {
                // modifying 'self' to improve self[i+1.., j]
                for k in i+1..self.size().0 {
                    if !self[(i,j)].divides(self[(k,j)]) {
                        let a = self[(i,j)];
                        let b = self[(k,j)];
                        let (gcd, x, y) = bezout_identity(a, b);
                        // Now, a*x+b*y==gcd
                        let s = a.euclid_div(gcd).0;
                        let t = b.euclid_div(gcd).0;

                        let ith_row = &self[(i, j..)]*x + &self[(k, j..)]*y;
                        let kth_row = &self[(i, j..)]*(-t) + &self[(k, j..)]*s;
                        self[(i, j..)].write_and_move( ith_row );
                        self[(k, j..)].write_and_move( kth_row );

                        let ith_row = &r_op[(i, ..)]*x + &r_op[(k, ..)]*y;
                        let kth_row = &r_op[(i, ..)]*(-t) + &r_op[(k, ..)]*s;
                        r_op[(i, ..)].write_and_move( ith_row );
                        r_op[(k, ..)].write_and_move( kth_row );

                        let ith_col = &r_op_inv[(.., i)]*s + &r_op_inv[(.., k)]*t;
                        let kth_col = &r_op_inv[(.., i)]*(-y) + &r_op_inv[(.., k)]*x;
                        r_op_inv[(.., i)].write_and_move( ith_col );
                        r_op_inv[(.., k)].write_and_move( kth_col );
                        
                        debug_assert!(self[(i,j)]==gcd);
                    }
                }
                // now the "self[i, j]" must divides all the elements in "self[i+1.., j]"
                debug_assert!( (i+1..self.size().0).all(|k| self[(i,j)].divides(self[(k,j)])),
                    "'self[{i}, {j}]' must divides all the elements in 'self[{i}+1.., {j}]', but it does not: \n{:?}", self
                );

                // modifying 'self' to improve self[i, j+1..]
                for l in j+1..self.size().1 {
                    if !self[(i,j)].divides(self[(i,l)]) {
                        let a = self[(i,j)];
                        let b = self[(i,l)];
                        let (gcd, x, y) = bezout_identity(a, b);
                        let s = a.euclid_div(gcd).0;
                        let t = b.euclid_div(gcd).0;

                        let jth_row = &self[(i.., j)] * x + &self[(i..,l)] * y;
                        let lth_row = &self[(i.., j)] * (-t) + &self[(i.., l)] * s;
                        self[(i..,j)].write_and_move( jth_row );
                        self[(i..,l)].write_and_move( lth_row );
                        
                        let jth_row = &c_op[(.., j)]*x + &c_op[(.., l)]*y;
                        let lth_row = &c_op[(.., j)]*(-t) + &c_op[(.., l)]*s;
                        c_op[(..,j)].write_and_move( jth_row );
                        c_op[(..,l)].write_and_move( lth_row );
                        
                        let jth_col = &c_op_inv[(j, ..)] * s + &c_op_inv[(l, ..)]*t;
                        let lth_col = &c_op_inv[(j, ..)] * (-y) + &c_op_inv[(l, ..)] * x;
                        c_op_inv[(j,..)].write_and_move( jth_col );
                        c_op_inv[(l,..)].write_and_move( lth_col );

                        debug_assert!(self[(i,j)]==gcd, "self[(i,j)]={}, gcd={}", self[(i,j)], gcd);
                    }
                }
                // now the "self[i, j]" must divides all the elements in "self[i, j+1..]"
                debug_assert!( (j+1..self.size().1).all(|l| self[(i,j)].divides(self[(i,l)])),
                    "'self[{i}, {j}]' must divides all the elements in 'self[{i}, {j}+1]', but it does not: \n{:?}", self
                );
            }

            // Eliminating entries
            for k in i+1..self.size().0 {
                if self[(k,j)] == T::zero() { continue; }

                let ratio = self[(k,j)].euclid_div(self[(i,j)]).0;
                self.row_operation(k, i, -ratio);
                r_op.row_operation(k, i, -ratio);
                r_op_inv.col_operation(i, k, ratio);
            }

            for l in j+1..self.size().1 {
                if self[(i,l)] == T::zero() { continue; }

                let ratio = self[(i,l)].euclid_div(self[(i,j)]).0;
                self.col_operation(l, j, -ratio);
                c_op.col_operation(l, j, -ratio);
                c_op_inv.row_operation(j, l, ratio);
            }

            i += 1; j+= 1;
        };


        // Bringing all the pivots to the diagonal
        for r in 0..std::cmp::min(self.size().0, self.size().1) {
            // if the entry at (r, r) is zero then we do nothing 
            if self[(r,r)] != T::zero() { continue; }

            // otherwise, we fetch the nonzero element and swap it with the entry at (r, r)
            let s = (r+1..self.size().1).find(|&s| self[(r,s)] != T::zero() );
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
        let rank = match (0..std::cmp::min(self.size().0, self.size().1)).find(|&r| self[(r,r)] == T::zero()) {
            Some(r) => r,
            None => std::cmp::min(self.size().0, self.size().1),
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

                let rth_row = &self[(r, ..)]*x + &self[(s, ..)]*y;
                let sth_row = &self[(r, ..)]*(-w) + &self[(s, ..)]*z;
                self[(r, ..)].write_and_move( rth_row );
                self[(s, ..)].write_and_move( sth_row );
                
                let rth_row = &r_op[(r, ..)]*x + &r_op[(s, ..)]*y;
                let sth_row = &r_op[(r, ..)]*(-w) + &r_op[(s, ..)]*z;
                r_op[(r, ..)].write_and_move( rth_row );
                r_op[(s, ..)].write_and_move( sth_row );

                let rth_col = &r_op_inv[(.., r)]*z + &r_op_inv[(.., s)]*w;
                let sth_col = &r_op_inv[(.., r)]*(-y) + &r_op_inv[(.., s)]*x;
                r_op_inv[(.., r)].write_and_move( rth_col );
                r_op_inv[(.., s)].write_and_move( sth_col );

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
impl <T:Field+Copy> Matrix<T> {
    pub fn rank_as_linear_map(&self) -> usize {
        let mut clone = self.clone();
        let shorter_side = std::cmp::min(self.size().0, self.size().1);
        let mut skipped = 0;
        for k in 0..shorter_side {
            if clone[(k,k)] == T::zero() {
                let result = (k+2..clone.size().1).find(|&i| clone[(k, i)] != T::zero() );
                match result {
                    Some(j) => clone.col_operation(k, j, T::one()),
                    None => (),
                }
            }
            if clone[(k,k)] == T::zero() {
                let result = (k+2..clone.size().0).find(|&j| clone[(j, k)] != T::zero() );
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

            for l in k+1..self.size().0 {
                let minus_one = T::zero() - T::one();
                let a_lk = clone[(l,k)];
                if a_lk == T::zero() { continue; } 
                clone.row_operation(l, k, minus_one * a_lk);
            }
        };

        shorter_side - skipped
    }
}


use super::Vector;
// Change of basis
impl<Basis> ops::Mul<Vec<Basis>> for Matrix<i64> 
    where Basis: Clone + PartialEq + PartialOrd, i64: ops::Mul<Basis, Output = Vector<i64, Basis>>
{
    type Output = Vec<Vector<i64, Basis>>;

    fn mul(self, rhs: Vec<Basis>) -> Self::Output {
        assert!(self.size().1==rhs.len(), "Multiplication failed: dimensions do not match: the left has size {:?} but the right has size {}", self.size, rhs.len());
        let mut out = Vec::new();
        for i in 0..self.size().0 {
            out.push( (0..self.size().1).map(|j| self[(i, j)]*rhs[j].clone() ).sum::<Vector<_, _>>() );
        }
        out
    }
}

// numerical linear algebra
impl Matrix<f64> {
    pub fn det(mut self) -> f64 {
        if self.size().0 != self.size().1 {
            return 0_f64;
        };

        // this implementation of determinant is based on LU decomposition
        self.lu();
        (0..self.size().0)
            .map(|i| self[(i,i)])
            .product()
    }

    pub fn random_orthogonal(n: usize) -> Matrix<f64> {
        let mut out = Matrix::zero(n ,n);
        for i in 0..n {
            let mut random = Matrix::random(n-i,1);
            random /= random.two_norm();
            out[(i..,i)].write_and_move( random );
        }
        for i in (0..n-1).rev() {
            let (v, b) = out[(i..,i)].householder_vec();
            let v = &*v;
            let update = &out[(i.., i+1..)] - v*b * (v.transpose() * &out[(i.., i+1..)]);
            out[(i.., i+1..)].write_and_move( update );
        }
        out
    }

    pub fn random_unit_vec(n: usize) -> Matrix<f64> {
        let mut out = Matrix::new(n ,1);
        let mut r = thread_rng();
        for i in 0..n {
            out[(i,0)] = r.gen_range(-1.0..1.0);
        }
        &*out / (&*out).two_norm()
    }

    pub fn random_sparse_symmetric(n: usize, m: usize) -> Matrix<f64> {
        let mut out = Matrix::random_symmetric(n ,m);
        let mut rng = rand::thread_rng();
        for i in 0..n {
            for j in i..m {
                let a: usize = rng.gen::<usize>() % 4;
                if a!=0 {
                    out[(i,j)] = 0.0;
                    out[(j,i)] = 0.0;
                }
            }
        }
        out
    }

    pub fn gram_schmidt(mut self) -> Matrix<f64> {
        // Normalize the first colmnn
        let factor = self[(..,0)].two_norm();
        self[(..,0)] /= factor;

        for i in 1..self.size().1 {
            for j in 0..i {
                let update = &self[(..,i)] - &self[(..,j)] * (self[(..,i)].transpose().dot(&self[(..,j)]));
                self[(..,i)].write( &update );
            }
            let factor = self[(..,i)].two_norm();
            self[(..,i)] /= factor;
        }

        self
    }

    pub fn gram_schmidt_without_normalize(mut self) -> Matrix<f64> {
        for i in 1..self.size().1 {
            for j in 0..i {
                let update = &self[(..,i)] - &self[(..,j)] * (self[(..,i)].transpose().dot(&self[(..,j)])) / (self[(..,j)].transpose().dot(&self[(..,j)]));
                self[(..,i)].write( &update );
            }
        }

        self
    }

    fn gram_schmidt_coeff_update(&mut self, basis: &Matrix<f64>, orthogonal_basis: &Matrix<f64>) {
        for i in 0..basis.size().1 {
            for j in 0..i {
                self[(i,j)] = basis[(..,i)].transpose().dot(&orthogonal_basis[(..,j)])
                    / orthogonal_basis[(..,j)].transpose().dot(&orthogonal_basis[(..,j)]);
            }
        }
    }

    pub fn lll_lattice_reduce(mut self, delta: f64) -> Matrix<f64> {
        assert!(self.size().0 >= self.size().1, "input basis invalid.");
        assert!(0.25 < delta && delta < 1.0, "'delta' must satisfy 0.25 < 'delta' < 1.0, but it is 'delta'");

        let mut orthogonal_basis = self.clone().gram_schmidt_without_normalize();
        let mut gram_schmidt_coeff = Self::zero(self.size().1, self.size().1);
        gram_schmidt_coeff.gram_schmidt_coeff_update(&self, &orthogonal_basis);

        let mut i = 1;
        while i < self.size().1 {
            for j in (0..i).rev() {
                if gram_schmidt_coeff[(i,j)].abs() <= 0.5 { continue; }

                let update = &self[(..,i)] - &self[(..,j)] * gram_schmidt_coeff[(i,j)].round();
                self[(..,i)].write( &update );

                // ToDo: change the following code so that it does the necessary updates only.
                orthogonal_basis = self.clone().gram_schmidt_without_normalize();
                gram_schmidt_coeff.gram_schmidt_coeff_update(&self, &orthogonal_basis)
            }

            let i_norm_squared = orthogonal_basis[(..,i)].transpose().dot( &orthogonal_basis[(..,i)] );
            let i_minus_1_norm_squared = orthogonal_basis[(..,i-1)].transpose().dot( &orthogonal_basis[(..,i-1)] );
            if i_norm_squared > i_minus_1_norm_squared * (delta - gram_schmidt_coeff[(i,i-1)].powi(2)) {
                i = i+1;
            } else {
                self.swap_cols(i-1, i);
                orthogonal_basis = self.clone().gram_schmidt_without_normalize();
                gram_schmidt_coeff.gram_schmidt_coeff_update(&self, &orthogonal_basis);
                i = std::cmp::max(i-1, 1);
            }
        }

        self
    }

    pub fn householder_vec(mut self) -> (Matrix<f64>, f64) {
        if self.size().1 > 1 {
            panic!("this function takes a column vector only, but the input has size {:?}.", self.size());
        }
        
        let a: f64 = self[(1.., 0)].transpose().dot( &self[(1.., 0)] );
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

    pub fn givens_rotation(a: f64, b: f64) -> (f64, f64) {
        if b.abs()<0.000000000001 {
            return (1.0, 0.0)
        }

        if a.abs() < b.abs() {
            let c = - a/b; 
            let y = 1.0/((1.0+c*c).sqrt());
            let x = c*y;
            (x, y)
        } else {
            let c = - b/a; 
            let x = 1.0/((1.0+c*c).sqrt());
            let y = c*x;
            (x,y)
        }
    }

    pub fn solve(mut self, mut b: Matrix<f64>) -> Self {
        assert_eq!(
            self.size().0, 
            self.size().1, 
            "The first input to this function must be a square matrix,but it is not: size={:?}.",
            self.size
        );
        assert_eq!(
            self.size().1, 
            b.size().0, 
            "The number of rows of the first matrix and the number of cols of the second matrix have to match, but self.size = {:?}, b.size = {:?}.",
            self.size,
            b.size
        );
        let n = self.size().0;
        self[(..,..)].lu();
        
        // Solve Ly = b by forward substitution
        for i in 1..n {
            b[(i, 0)] = b[(i,0)] - self[(i, 0..i)].dot( &b[(0..i, 0)] );
        }

        // Solve Ux = y by backward substitution
        b[(n-1,0)] = b[(n-1,0)] / self[(n-1,n-1)];
        for i in (0..n-1).rev() {
            b[(i, 0)] = b[(i,0)] - self[(i, i+1..n)].dot( &b[(i+1..n, 0)]);
            b[(i, 0)] /= self[(i,i)];
        }
        
        b
    }

    pub fn spectrum(mut self) -> Matrix<f64> {
        assert_eq!( self.size().0, self.size().1,
            "Input to this function must be a square matrix."
        );
        let n = self.size().0;
        self.qr();

        let mut out = vec!{0_f64; n};
        let mut i=0;
        while i < n {
            if i == n-1 || self[(i+1,i)].abs() < 0.00000000001 {
                out[i] = self[(i,i)];
                i+=1;
            } else {
                (out[i], out[i+1]) = self[(i..i+2, i..i+2)].real_spectrum_of_two_by_two().unwrap();
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
        assert_eq!( self.size().0, self.size().1,
            "Input to this function must be a square matrix."
        );
        let n = self.size().0;

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

            let mut update = m.clone().solve( eigen_vectors[(.., i)].as_matrix() );
            update /= update.two_norm();
            eigen_vectors[(.., i)].write_and_move( update );
            let mut update = m.solve( eigen_vectors[(.., i)].as_matrix() );
            update /= update.two_norm();
            eigen_vectors[(.., i)].write_and_move( update );
        }

        (spectrum, eigen_vectors)
    }


    pub fn spectrum_symmetric(mut self) -> Self {
        assert_eq!( self.size().0, self.size().1,
            "Input to this function must be a square matrix."
        );
        let n = self.size().0;
        self.symmetric_qr(false);

        let mut out = vec!{0_f64; n};
        let mut i=0;
        while i < n {
            if i == n-1 || self[(i+1,i)].abs() < 0.00000000000001 {
                out[i] = self[(i,i)];
                i+=1;
            } else {
                (out[i], out[i+1]) = self[(i..i+2, i..i+2)].real_spectrum_of_two_by_two().unwrap();
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

    pub fn spectrum_with_invariant_space_symmetric(self) -> (Matrix<f64>, Matrix<f64>) {
        assert_eq!( self.size().0, self.size().1,
            "Input to this function must be a square matrix."
        );
        let n = self.size().0;

        // initialize the eigenvectors
        let mut eigen_vectors = Matrix::<f64>::new(n, n);
        let elem = 1.0 / (n as f64).sqrt();
        for i in 0..n {
            for j in 0..n {
                eigen_vectors[(i,j)] = elem;
            }
        }

        // get the eigenvalues
        let spectrum = self.clone().spectrum_symmetric();

        // run the inverse iteration to get the eigenvector for each of the eigenvalues
        for i in 0..n {
            let lambda = spectrum[(i,0)];

            let mut m = self.clone();
            for j in 0..n {
                m[(j,j)] -= lambda;
            }

            let mut update = m.clone().solve( eigen_vectors[(.., i)].as_matrix() );
            update /= update.two_norm();
            eigen_vectors[(.., i)].write_and_move( update );
            let mut update = m.solve( eigen_vectors[(.., i)].as_matrix() );
            update /= update.two_norm();
            eigen_vectors[(.., i)].write_and_move( update );
        }

        (spectrum, eigen_vectors)
    }

    pub fn spectrum_with_n_smallest_eigenvecs_symmetric(self, n: usize) -> (Matrix<f64>, Matrix<f64>) {
        assert_eq!( self.size().0, self.size().1,
            "Input to this function must be a square matrix."
        );
        let size = self.size().0;

        // initialize the eigenvectors
        let mut eigen_vectors = Matrix::<f64>::new(size, n);
        let elem = 1.0 / (size as f64).sqrt();
        for i in 0..size {
            for j in 0..n {
                eigen_vectors[(i,j)] = elem;
            }
        }

        // get the eigenvalues
        let spectrum = self.clone().spectrum_symmetric();

        // run the inverse iteration to get the eigenvector for each of the eigenvalues
        for i in 0..n {
            let lambda = spectrum[(i,0)];

            let mut m = self.clone();
            for j in 0..size {
                m[(j,j)] -= lambda;
            }

            let mut update = m.clone().solve( eigen_vectors[(.., i)].as_matrix() );
            update /= update.two_norm();
            eigen_vectors[(.., i)].write_and_move( update );
            let mut update = m.solve( eigen_vectors[(.., i)].as_matrix() );
            update /= update.two_norm();
            eigen_vectors[(.., i)].write_and_move( update );
        }

        (spectrum, eigen_vectors)
    }

    // pub fn extremal_spectrum_symmetric(self, n_eigen: usize, small: bool) -> (Matrix<f64>, Matrix<f64>) {
    //     let tolerance = 0.000000000001;
    //     let (mut a, mut b, op) = self.lanzos_tridiagonalization();
        
    //     let n = a.size().0;

    //     let mut q=0;
    //     while q<n {
    //         for i in 0..n-1 {
    //             if b[(i,0)].abs()<tolerance {
    //                 b[(i,0)] = 0.0;
    //             }
    //         }

    //         q = (0..n-2)
    //             .find( |&i| b[(n-i-2, 0)].abs()>tolerance && b[(n-i-3, 0)].abs()>tolerance )
    //             .unwrap_or(n);

    //         let p = (1..n-q)
    //             .rev()
    //             .find(|&i| b[(i-1, 0)].abs()<tolerance )
    //             .unwrap_or(0);

    //         if q<n {
    //             debug_assert!( (n-q)-p > 2);

    //             let mut tmp = Matrix::zero(n-q-p, n-q-p);
    //             for i in 0..tmp.size().0 {
    //                 tmp[(i,i)] = a[(p+i, 0)];
    //                 if i<tmp.size().0-1 {
    //                     tmp[(i,i+1)] = b[(p+i, 0)];
    //                     tmp[(i+1,i)] = b[(p+i, 0)];
    //                 }
    //             }

    //             let _givens = tmp.implicit_symmetric_qr_with_wilkinson_shift(false);

    //             for i in 0..tmp.size().0 {
    //                 a[(p+i, 0)] = tmp[(i,i)];
    //                 if i<tmp.size().0-1 {
    //                     b[(p+i, 0)] = tmp[(i,i+1)];
    //                 }
    //             }
    //         }
    //     }

    //     let mut spectrum = vec!{0_f64; n};
    //     let mut i=0;
    //     while i < n {
    //         if i == n-1 || self[(i+1,i)].abs() < 0.00000000001 {
    //             spectrum[i] = self[(i,i)];
    //             i+=1;
    //         } else {
    //             (spectrum[i], spectrum[i+1]) = self[(i..i+2, i..i+2)].real_spectrum_of_two_by_two().unwrap();
    //             i+=2;
    //         }
    //     }

    //     spectrum.sort_by(|x, y| x.partial_cmp(y).unwrap());
    //     let spectrum = {
    //         let mut a = Matrix::new(n_eigen,1);
    //         let offset = if small{0} else {n-n_eigen};
    //         (0..n_eigen).for_each( |i| a[(i,0)]=spectrum[i+offset] );
    //         a
    //     };

    //     // initialize the eigenvectors
    //     let mut eigenvectors = Matrix::<f64>::new(n, n_eigen);
    //     let elem = 1.0 / (n as f64).sqrt();
    //     for i in 0..n {
    //         for j in 0..n_eigen {
    //             eigenvectors[(i,j)] = elem;
    //         }
    //     }

    //     // run the inverse iteration to get the eigenvector for each of the eigenvalues
    //     for i in 0..n_eigen {
    //         let offset = if small{0} else {n-n_eigen};
    //         let lambda = spectrum[(i+offset,0)];

    //         let mut m = self.clone();
    //         for j in 0..n {
    //             m[(j,j)] -= lambda;
    //         }

    //         let mut update = m.clone().solve( eigenvectors[(.., i)].as_matrix() );
    //         update /= update.two_norm();
    //         eigenvectors[(.., i)].write_and_move( update );
    //         let mut update = m.solve( eigenvectors[(.., i)].as_matrix() );
    //         update /= update.two_norm();
    //         eigenvectors[(.., i)].write_and_move( update );
    //     }

    //     // recover the original eigenvectors
    //     let eigenvectors = op * eigenvectors;

    //     (spectrum, eigenvectors)
    // }
}

impl SubMatrix<f64> {
    pub fn det(&self) -> f64 {
        // this implementation of determinant is based on the LU decomposition
        self.as_matrix().det()
    }

    pub fn one_norm(&self) -> f64 {
        (0..self.size().1)
            .map(|j| (0..self.size().0).map(|i| self[(i ,j)].abs() ).sum() )
            .fold(0.0, |max, val| if max<=val { val } else { max }  )
    }

    pub fn two_norm(&self) -> f64 {
        if self.size().1 ==  1 {
            ((0..self.size().0).map(|i| self[(i,0)]*self[(i,0)] ).sum::<f64>()).sqrt()
        } else {
            panic!("2-norm for matrices are not supported yet.");
        }
    }

    pub fn inf_norm(&self) -> f64 {
        (0..self.size().0)
            .map(|i| (0..self.size().1).map(|j| self[(i ,j)].abs() ).sum() )
            .fold(0.0, |max, val| if max<=val { val } else { max } )
    }

    pub fn frobenius_norm(&self) -> f64 {
        (0..self.size().0)
            .map( |i|
                (0..self.size().1).map(|j| self[(i,j)]*self[(i,j)] ).sum::<f64>()
            )
            .sum::<f64>()
            .sqrt()
    }

    pub fn householder_vec(&self) -> (Matrix<f64>, f64) {
        self.as_matrix().householder_vec()
    }

    fn real_spectrum_of_two_by_two(&self) -> Option<(f64, f64)> {
        debug_assert_eq!( self.size(), (2,2) );
        let det = self[(0,0)] * self[(1,1)] - self[(0,1)] * self[(1,0)];
        let m = (self[(0,0)] + self[(1,1)]) / 2.0;
        if m*m-det > 0.0 {
            let rad = (m*m-det).sqrt();
            Some((m-rad, m+rad))
        } else {
            None
        }
    }

    pub fn lu(&mut self) {
        assert_eq!(self.size().0, self.size().1, "Input has to this function has to be a square matrix.");

        let n = self.size().0;
        for i in 0..n-1 {
            let update = &self[(i+1..n, i)] / self[(i,i)];
            self[(i+1..n,i)].write_and_move( update );
            let update = &self[(i+1..n,i+1..n)] - &self[(i+1..n,i)] * &self[(i, i+1..n)];
            self[(i+1..n, i+1..n)].write_and_move( update );
        }
    }


    pub fn hessenberg_form(&mut self) {
        if self.size().0 != self.size().1 {
            panic!("Input to this function must be a square matrix.");
        }

        if self.size().0 < 2 {
            return;
        }

        let n = self.size().0; // = self.size().1
        for k in 0..n-2 {
            let (v, b) = self[(k+1..n, k)].householder_vec();
            let v = &*v;
            let v_t = v.transpose();

            let update = &self[(k+1..n, k..n)] - v*(v_t * &self[(k+1..n, k..n)])*b;
            self[(k+1..n, k..n)].write_and_move( update );

            let update = &self[(0..n, k+1..n)] - (&self[(0..n, k+1..n)]*v)*v_t*b;
            self[(0..n, k+1..n)].write_and_move( update );
        };
    }

    pub fn raw_francis_step(&mut self, return_op: bool) -> Matrix<f64> {
        assert_eq!( self.size().0, self.size().1,
            "Input to this function must be a square matrix."
        );
        let n = self.size().0;

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
            let v = &*v;
            let q = if k>0 { k-1 } else { 0 };
            let update = &self[(k..k+3, q..n)] -
                v * (v.transpose() * &self[(k..k+3, q..n)] * b);
            self[(k..k+3, q..n)].write_and_move( 
                update
            );

            let r = if k<n-3 { k+4 } else { n };
            let update = &self[(0..r, k..k+3)] -
                (&self[(0..r, k..k+3)] * v) * (v.transpose()*b);
            self[(0..r, k..k+3)].write_and_move(
                update
            );
            x = self[(k+1, k)];
            y = self[(k+2, k)];
            if k < n-3 {
                z = self[(k+3, k)];
            }
        }

        let (v, b) = crate::matrix!{f64; [[x,y]]}.transpose().householder_vec();
        let v = &*v;
        let update = &self[(n-2..n, n-3..n)] -
            v * (v.transpose() * &self[(n-2..n, n-3..n)] * b);
        self[(n-2..n, n-3..n)].write_and_move( update );

        let update = &self[(0..n, n-2..n)] -
            (&self[(0..n, n-2..n)] * v) * (v.transpose()*b);
        self[(0..n, n-2..n)].write_and_move( update );

        op
    }

    pub fn francis_step(&mut self) {
        self.raw_francis_step(false);
    }

    pub fn francis_step_with_op(&mut self) -> Matrix<f64> {
        self.raw_francis_step(true)
    }

    fn raw_qr(&mut self, returns_op: bool) -> Matrix<f64> {
        assert_eq!( self.size().0, self.size().1,
            "Input to this function must be a square matrix."
        );
        let n = self.size().0;

        if self.size().0<2 {
            return Matrix::identity(1);
        }

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
                    op = op * self[(p..n-q, p..n-q)].francis_step_with_op();
                } else {
                    self[(p..n-q, p..n-q)].francis_step();
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

    pub fn householder_tridiagonalization(&mut self, with_op: bool) {
        let n = self.size().0;
        for i in 0..n-2 {
            let (v, b) = self[(i+1.., i)].householder_vec();
            let v = &*v;
            let p = (&self[(i+1.., i+1..)] * v) * b;
            let p = &*p;
            let w = p - v * b/2.0*(v.transpose().dot(p));
            let w = &*w;
            self[(i+1, i)] = self[(i+1.., i)].two_norm();
            self[(i, i+1)] = self[(i+1, i)];


            // computing 
            //  - let update = &self[(i+1.., i+1..)] - v*w.transpose() - w*v.transpose();
            //  - self[(i+1.., i+1..)].write( &*update );
            // but with exploiting the symmetricity
            for j in 0..n-(i+1) {
                self[(j+i+1, j+i+1)] -= 2.0*v[(j,0)]*w[(j,0)];
                for k in j+1..n-(i+1) {
                    self[(j+i+1,k+i+1)] -= v[(j,0)]*w[(k,0)] + v[(k,0)]*w[(j,0)];
                    self[(k+i+1,j+i+1)] = self[(j+i+1,k+i+1)];
                }
            }
            if with_op {
                self[(i+2.., i)].write( &v[(1..,0)] );
            }
        }
    }

    pub fn extract_householder_vecs(&self) -> Matrix<f64> {
        assert_eq!( self.size().0, self.size().1,
            "Input to this function must be a square matrix."
        );
        let n = self.size().0;

        if n == 1 {
            let x = self[(0,0)];
            let s = 2.0/(x*x+1.0);
            matrix!(f64;
                [[1.0-s, -s*x],
                 [-s*x, 1.0-s*x*x]]
            )
        } else {
            let mut m = Matrix::identity(n+1);
            let v = &self[(..,0)];
            let s = 2.0 / (v.transpose().dot(v)+1.0);
            let q = self[(1.., 1..)].extract_householder_vecs();
            let q = &*q;

            m[(0,0)] = 1.0-s;
            m[(1..,0)].write(&(q * v *(-s)));
            m[(0,1..)].write(&(v.transpose()*(-s)));
            m[(1..,1..)].write(&(q-(q*v)*(v.transpose()*s)));
            m
        }
    }


    pub fn implicit_symmetric_qr_with_wilkinson_shift(&mut self, with_op: bool) -> Matrix<f64> {
        let n = self.size().0;
        let d = (self[(n-2, n-2)] - self[(n-1, n-1)])/2.0;
        let sign_d = if d>=0.0 { 1.0 } else { -1.0 };
        let m = self[(n-1, n-1)]-self[(n-1, n-2)].powi(2)/(d+sign_d*(d*d+self[(n-1, n-2)]*self[(n-1, n-2)]).sqrt());
        let mut x = self[(0,0)] - m;
        let mut z = self[(1,0)];
        
        let mut op = Matrix::zero(n-1, 2);

        for i in 0..n-1 {
            let (c,s) = Matrix::<f64>::givens_rotation(x,z);
            let givens =  matrix!(f64; [[c, s],[-s, c]]);
            let givens = &*givens;

            let start = if i==0 {0} else {i-1};
            let end = n.min(start+4);
            let update = givens.transpose() * &self[(i..i+2, start..end)];
            self[(i..i+2, start..end)].write(&update);
            let update = &self[(start..end, i..i+2)] * givens;
            self[(start..end, i..i+2)].write(&update);


            if i < n-2 {
                x = self[(i+1, i)];
                z = self[(i+2, i)];
            }

            if with_op {
                op[(i,0)] = c;
                op[(i,1)] = s;
            }
        }
        op
    }

    pub fn symmetric_qr(&mut self, with_op: bool) -> Matrix<f64> {
        let n = self.size().0;
        let tolerance = 0.0000000001;
        let mut op = Matrix::identity(n);
        self.householder_tridiagonalization(with_op);
        if with_op {
            op = Matrix::identity(n);
            op[(1..,1..)].write( &(self[(2..,..n-2)].extract_householder_vecs()) );
        };

        for i in 0..n-2 {
            for j in i+2..n{
                self[(i, j)] = 0.0;
                self[(j, i)] = 0.0;
            }
        }

        let mut q=0;
        while q<n {
            for i in 0..n-1 {
                if self[(i+1,i)].abs()<tolerance{
                    self[(i+1, i)] = 0.0;
                    self[(i, i+1)] = 0.0;
                }
            }

            q = (0..n-2)
                .find( |&i| self[(n-i-1,n-i-2)].abs()>tolerance && self[(n-i-2,n-i-3)].abs()>tolerance )
                .unwrap_or(n);

            let p = (1..n-q)
                .rev()
                .find(|&i| self[(i,i-1)].abs()<tolerance )
                .unwrap_or(0);

            if q<n {
                debug_assert!( (n-q)-p > 2);
                let givens = self[(p..n-q, p..n-q)].implicit_symmetric_qr_with_wilkinson_shift(with_op);
                if with_op {
                    for (i, c,s) in (0..givens.size.0).map(|i| (i, givens[(i,0)], givens[(i,1)])) {
                        let givens =  matrix!(f64; [[c, -s],[s, c]]);
                        let givens = &*givens;
                        let update = givens * &op[(i..i+2, ..)];
                        op[(i..i+2, ..)].write(&update);
                    }
                }
            }
        }
        op
    }

    // pub fn lanzos_tridiagonalization(&self) -> (Matrix<f64>, Matrix<f64>, Matrix<f64>) {
    //     let mut q = vec![Matrix::zero(self.size().0, 1)];
    //     q.reserve(self.size().0);
    //     q.push({
    //         let mut m = Matrix::new(self.size().0, 1);
    //         let val = 1.0/((self.size().0) as f64).sqrt();
    //         for i in 0..self.size().0 {
    //             m[(i,0)] = val;
    //         }
    //         m
    //     });

    //     let mut alpha: Vec<f64> = vec![0.0];
    //     alpha.reserve(self.size().0);
    //     let mut beta: Vec<f64> = vec![1.0];
    //     beta.reserve(self.size().0);
        
    //     let mut r = q[1].clone();
        
    //     let mut i = 0;
    //     while (beta[i]).abs() > 0.00001 && i<20 {
    //         if i!=0 {
    //             q.push(r / beta[i]);
    //         }
    //         i = i + 1;
    //         alpha.push(((&*q[i]).transpose() * self).dot( &*q[i] ));
    //         debug_assert_eq!(alpha.len(), i+1);
    //         r = {
    //             let mut tmp = self * &*q[i];
    //             for j in 0..tmp.size().0 {
    //                 tmp[(j,0)] -= alpha[i] * q[i][(j,0)] + beta[i-1]*q[i-1][(j,0)];
    //             }
    //             tmp
    //         };
    //         beta.push( r.two_norm() );
    //         debug_assert_eq!(beta.len(), i+1);

    //         println!("beta={}", beta[i]);
    //     }

    //     let mut q_out = Matrix::new(self.size().0, q.len()-1);
    //     for (j, v) in q.into_iter().skip(1).enumerate() {
    //         q_out[(.., j)].write(&v);
    //     }

    //     let mut alpha_out = Matrix::new(i, 1);
    //     for j in 1..=i {
    //         alpha_out[(j-1, 0)] = alpha[j];
    //     }

    //     let mut beta_out = Matrix::new(i-1, 1);
    //     for j in 1..i {
    //         beta_out[(j-1, 0)] = beta[j];
    //     }

    //     (alpha_out, beta_out, q_out)
    // }

    // pub fn lanzos_tridiagonalization(&self) -> (Matrix<f64>, Matrix<f64>, Matrix<f64>) {
    //     let n = self.size().0;
    //     let mut alpha = Vec::new();
    //     let mut beta = Vec::new();
    //     let mut q = Vec::new();

    //     // let mut w = {
    //     //     let mut m = Matrix::new(n, 1);
    //     //     let val = 1.0/ ((n as f64).sqrt());
    //     //     for i in 0..n {
    //     //         m[(i,0)] = val;
    //     //     }
    //     //     m
    //     // };

    //     let mut w = {
    //         let mut m = Matrix::zero(n, 1);
    //         m[(0,0)] = 1.0;
    //         m
    //     };

    //     q.push(w.clone());
    //     let mut v = self * &*w;

    //     alpha.push((&*w).transpose().dot(&v));

    //     v = v - (&*w) * *alpha.last().unwrap();

    //     beta.push(v.two_norm());

    //     let mut i=0;

    //     while *beta.last().unwrap() > 0.00000001 && i<self.size().0-1 {
    //         println!( "beta[{i}]={}", *beta.last().unwrap() );
    //         for j in 0..n {
    //             let t = w[(j,0)];
    //             w[(j,0)] = v[(j,0)] / *beta.last().unwrap();
    //             v[(j,0)] = - beta.last().unwrap() * t;
    //         }
    //         v = v + self * &*w;
    //         i+=1;
    //         alpha.push((&*w).transpose().dot(&v));
    //         v = v - &*w * *alpha.last().unwrap();
    //         beta.push(v.two_norm());

    //         q.push(w.clone());
    //     }

    //     let a_out = Matrix::from_array(vec![alpha]).transpose();
    //     beta.pop();
    //     let b_out = Matrix::from_array(vec![beta]).transpose();
    //     let q_out = {
    //         let mut q_out = Matrix::new(self.size().0, q.len());
    //         for i in 0..q.len() {
    //             q_out[(.., i)].write( &q[i] );
    //         }
    //         q_out
    //     };

    //     // for i in 0..q_out.size().1 {
    //     //     for j in i+1..q_out.size().1 {
    //     //         let p = q_out[(.., i)].transpose().dot(&q_out[(..,j)]);
    //     //         assert!(p<0.001, "p={p}");
    //     //     }
    //     // }
        
    //     (a_out, b_out, q_out)
    // }
}

impl<T: PID+Display> std::fmt::Display for SubMatrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision().unwrap_or(3);
        write!(f, "\n")?;
        for i in 0..self.size().0 {
            write!(f, "| ")?;
            for j in 0..self.size().1 {
                write!(f, "{:.decimals$}", self[(i,j)])?;
                if j != self.size().1-1 {
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
                    self[(.., ..)].$fn()
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
                    self[(.., ..)].$fn();
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
//                     self[(.., ..).$fn()
//                 }
//             ) *
//         }
//     }
// }

// fn_cast_impl_for_matrix!{ 
// }
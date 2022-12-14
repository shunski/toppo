use crate::commutative::{PID, Field, bezout_identity};

use std::alloc::{self, Layout};
use std::marker::PhantomData;
use std::mem;
use std::ops::{Deref, DerefMut};
use std::ptr::{self, NonNull};
use std::fmt;
use std::cmp::PartialEq;
use std::ops::Mul;
use std::ops::Add;


pub struct Matrix<T: PID> {
    ptr: NonNull<T>,
    _marker: PhantomData<T>,
    
    n_row: usize,
    n_col: usize,
}

unsafe impl<T: Send + PID> Send for Matrix<T> {}
unsafe impl<T: Sync + PID> Sync for Matrix<T> {}

#[allow(unused)]
impl<T: PID> Matrix<T> {
    fn new(n: usize, m: usize) -> Self {
        let mut mat = Matrix {
            ptr: NonNull::dangling(),
            _marker: PhantomData,
            n_row: n,
            n_col: m,
        };
        mat.alloc();
        mat
    }

    pub fn zero(n: usize, m: usize) -> Self {
        if n==0 || m==0 { panic!("Any matrix must have positive dimension."); }
        let mut mat = Matrix {
            ptr: NonNull::dangling(),
            _marker: PhantomData,
            n_row: n,
            n_col: m,
        };
        
        mat.alloc();
        mat.set_all_to(T::zero());
        mat
    }

    pub fn identity(size: usize) -> Self {
        let mut mat = Matrix::zero(size, size);
        (0..size).for_each(|i| mat.write(i,i,T::one()));
        mat
    }

    fn alloc(&mut self) {
        let layout = Layout::array::<T>(self.n_row * self.n_col).unwrap();
        
        // Ensure that the new allocation doesn't exceed `isize::MAX` bytes.
        assert!(
            layout.size() <= isize::MAX as usize,
            "Allocation too large"
        );

        let new_ptr = unsafe { alloc::alloc(layout) };

        // If allocation fails, `new_ptr` will be null, in which case we abort.
        self.ptr = match NonNull::new(new_ptr as *mut T) {
            Some(p) => p,
            None => alloc::handle_alloc_error(layout),
        };
    }

    fn set_all_to (&mut self, val: T) {
        let n_elem = self.n_row * self.n_col;
        (0..n_elem).for_each( |i|
            unsafe { ptr::write(self.ptr().add(i), val); }
        );
    }

    pub fn write (&mut self, i: usize, j: usize, val: T) {
        if i>=self.n_row || j>=self.n_col { 
            panic!("Idx out of bounds: index: ({}, {}) while (ros, col): ({}, {})", i, j, self.n_row, self.n_col); 
        }
        unsafe { ptr::write(self.ptr().add(i*self.n_col + j), val); }
    }

    pub fn get (&self, i: usize, j: usize) -> T {
        if i>=self.n_row || j>=self.n_col { 
            panic!("Idx out of bounds: index: ({}, {}) while (ros, col): ({}, {})", i, j, self.n_row, self.n_col); 
        }
        unsafe { ptr::read(self.ptr().add(i*self.n_col + j)) }
    }

    fn write_row(&mut self, i: usize, new_row: Matrix<T>) {
        if i>=self.n_row { panic!(); }

        new_row.iter().enumerate().for_each(|(j, &x)| self.write(i, j, x));
    }

    fn write_col(&mut self, j: usize, new_col: Matrix<T>) {
        if j>=self.n_col { panic!(); }

        new_col.iter().enumerate().for_each(|(i, &x)| self.write(i, j, x));
    }

    fn clone_row(&self, i: usize) -> Self {
        let mut m = Matrix::new( 1, self.n_col );
        self.row_iter(i).enumerate().for_each(|(j, &x)| m.write(0, j, x));
        m
    }

    fn clone_col(&self, j: usize) -> Self {
        let mut m = Matrix::new( self.n_row, 1 );
        self.col_iter(j).enumerate().for_each(|(i, &x)| m.write(i, 0, x));
        m
    }

    fn scalar_mul(mut self, scalar: T) -> Self {
        self.iter_mut().for_each(|x| *x = scalar * *x );
        self
    }

    fn is_row_zero(&self, i: usize) -> bool {
        self.row_iter(i).all(|&x| x == T::zero() )
    }

    fn is_col_zero(&self, j: usize) -> bool {
        self.col_iter(j).all(|&x| x == T::zero() )
    }

    fn row_operation(&mut self, row1:usize, s1: T, row2: usize, s2: T) {
        if row1 >= self.n_row || row2 >= self.n_row { panic!(); }
        let new_row = self.clone_row(row1).scalar_mul(s1) + self.clone_row(row2).scalar_mul(s2);
        self.write_row(row1, new_row);
    }

    fn col_operation(&mut self, col1:usize, s1: T, col2: usize, s2: T) {
        if col1 >= self.n_col || col2 >= self.n_col { panic!(); }
        let new_col = self.clone_col(col1).scalar_mul(s1) + self.clone_col(col2).scalar_mul(s2);
        self.write_col(col1, new_col);
    }

    fn swap_rows(&mut self, row1: usize, row2: usize) {
        let tmp = self.clone_row(row1);
        self.write_row(row1, self.clone_row(row2));
        self.write_row(row2, tmp);
    }

    fn swap_cols(&mut self, col1: usize, col2: usize) {
        let tmp = self.clone_col(col1);
        self.write_col(col1, self.clone_row(col2));
        self.write_col(col2, tmp);
    }

    fn ptr(&self) -> *mut T {
        self.ptr.as_ptr()
    }

    pub fn size(&self) -> (usize, usize) {
        (self.n_row, self.n_col)
    }

    fn from_vec(v: Vec<T>) -> Matrix<T> {
        if v.len() == 0 { panic!(); }

        let mut mat = Matrix {
            ptr: NonNull::dangling(),
            _marker: PhantomData,
            n_row: v.len(),
            n_col: 1,
        };
        
        mat.alloc();
        v.iter().enumerate().for_each(|(i, &x)| mat.write(i, 1, x));
        mat
    }

    fn from_array(mut v: Vec<Vec<T>>) -> Matrix<T> {
        if v.len() == 0 || v[0].len() == 0 { panic!(); }

        if v.iter().any(|x| x.len() != v[0].len()) { 
            panic!("Input array is not a rectagle"); 
        }

        let mut mat = Matrix {
            ptr: NonNull::dangling(),
            _marker: PhantomData,
            n_row: v.len(),
            n_col: v[0].len(),
        };

        mat.alloc();
        v.iter().flatten().enumerate().for_each(|(i, &x)| mat.write(i / mat.n_col, i % mat.n_col, x));
        mat
    }
}

impl<'a, T> Matrix<T> 
    where T: 'a + PID
{
    fn row_iter(&'a self, i: usize) -> MatIter<'a, T>  {
        if i > self.n_row { panic!(); }
        
        let start = i * self.n_col;
        let end = (i+1) * self.n_col - 1;

        MatIter::new(&self[start..=end], 1, self.n_col)
    }

    fn col_iter(&'a self, j: usize) -> MatIter<'a, T> {
        if j > self.n_col { panic!(); }
        
        let start = j;
        let end = (self.n_row-1) * self.n_col + j;

        MatIter::new(&self[start..=end], self.n_col, self.n_row)
    }
}

impl<T: PID> Drop for Matrix<T> {
    fn drop(&mut self) {
        let elem_size = mem::size_of::<T>();

        if elem_size != 0 {
            unsafe {
                alloc::dealloc(
                    self.ptr.as_ptr() as *mut u8,
                    Layout::array::<T>(self.n_row * self.n_col).unwrap(),
                );
            }
        }
    }
}

impl<T: PID> Deref for Matrix<T> {
    type Target = [T];
    fn deref(&self) -> &[T] {
        unsafe { std::slice::from_raw_parts(self.ptr(), self.n_row * self.n_col) }
    }
}

impl<T: PID> DerefMut for Matrix<T> {
    fn deref_mut(&mut self) -> &mut [T] {
        unsafe { std::slice::from_raw_parts_mut(self.ptr(), self.n_row * self.n_col) }
    }
}

impl<T: PID> Clone for Matrix<T> {
    fn clone(&self) -> Self {
        let mut clone = Matrix::new( self.n_row, self.n_col );
        self.iter().enumerate().for_each( |(k, &x)| {
            let i = k / self.n_col; 
            let j = k % self.n_col; 
            clone.write(i,j,x);
        });
        clone
    }
}


impl<T: PID> PartialEq for Matrix<T> {
    fn eq(&self, rhs: &Self) -> bool {
        if self.n_row != rhs.n_row || self.n_row != rhs.n_row {
            return false;
        }
        self.iter().zip(rhs.iter()).all(|(&x, &y)| x == y)
    }
}

impl<T: PID + std::fmt::Debug + std::fmt::Display> fmt::Debug for Matrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\n")?;
        for i in 0..self.n_row {
            write!(f, "| ")?;
            for (j, &x) in self.row_iter(i).enumerate() {
                write!(f, "{}", x)?;
                if j != self.n_col-1 {
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


#[macro_export]
macro_rules! matrix {
    ( integer; [$([$( $elem:expr ), *]), *]) => {
        {
            let mut v: Vec<Vec<i64>> = Vec::new();
            $(
                v.push(vec!($( $elem ), *));
            )*
            Matrix::<i64>::from_array(v)
        }
    };

    ( rational; [$([$( $elem:expr ), *]), *]) => {
        {
            use crate::commutative::rational::Rational;
            let mut v: Vec<Vec<Rational>> = Vec::new();
            $(
                v.push(vec!($( Rational::from($elem) ), *));
            )*
            Matrix::<Rational>::from_array(v)
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


#[cfg(test)]
mod init_test {
    use crate::lin_alg::matrix::Matrix;
    #[test]
    fn init_test() {
        let m1 = Matrix::<i64>::zero(5,4);
        let m2 = matrix!(integer; 
            [[0,0,0,0],
             [0,0,0,0],
             [0,0,0,0],
             [0,0,0,0],
             [0,0,0,0]]
        );
        assert_eq!(m1, m2);

        let m1 = Matrix::<i64>::identity(5);
        let m2 = matrix!(integer; 
            [[1,0,0,0,0],
             [0,1,0,0,0],
             [0,0,1,0,0],
             [0,0,0,1,0],
             [0,0,0,0,1]]
        );
        assert_eq!(m1, m2);
    }
}




impl<T: PID> IntoIterator for Matrix<T> {
    type Item = T;
    type IntoIter = MatIntoIter<T>;
    fn into_iter(self) -> MatIntoIter<T> {
        unsafe {
            let iter = MatIntoIter::new(&self, 1);
            mem::forget(self);
            iter
        }
    }
}

struct MatIter<'a, T> 
    where T: 'a + PID
{
    data: &'a[T],
    offset: usize,
    curr: usize,
    size: usize,
}

impl<'a, T> MatIter<'a, T> 
    where T: 'a + PID
{
    fn new(slice: &'a[T], offset_in: usize, n: usize) -> MatIter<T> {
        MatIter{
            data: slice,
            curr: 0,
            offset: offset_in,
            size: n,
        }
    }
}

impl<'a, T> Iterator for MatIter<'a, T> 
    where T: 'a + PID
{
    type Item = &'a T;
    fn next(&mut self) -> Option<Self::Item> {
        let out = if self.curr < self.size {
            Some(&self.data[self.curr * self.offset])
        } else {
            None
        };
        self.curr += 1;
        out
    }
}


pub struct MatIntoIter<T: PID> {
    start: *const T,
    end: *const T,
    offset: isize,
}

impl<T: PID> MatIntoIter<T> {
    unsafe fn new(slice: &[T], offset_in: isize) -> Self {
        MatIntoIter {
            start: slice.as_ptr(),
            end: if mem::size_of::<T>() == 0 {
                ((slice.as_ptr() as usize) + slice.len()) as *const _
            } else if slice.len() == 0 {
                slice.as_ptr()
            } else {
                slice.as_ptr().add(slice.len())
            },
            offset: offset_in,
        }
    }
}

impl<T: PID> Iterator for MatIntoIter<T> {
    type Item = T;
    fn next(&mut self) -> Option<T> {
        if self.start == self.end {
            None
        } else {
            unsafe {
                if mem::size_of::<T>() == 0 {
                    self.start = (self.start as usize + 1) as *const _;
                    Some(ptr::read(NonNull::<T>::dangling().as_ptr()))
                } else {
                    let old_ptr = self.start;
                    self.start = self.start.offset(self.offset);
                    Some(ptr::read(old_ptr))
                }
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let elem_size = mem::size_of::<T>();
        let len = (self.end as usize - self.start as usize)
                  / if elem_size == 0 { 1 } else { elem_size };
        (len, Some(len))
    }
}


#[cfg(test)]
mod basic_functionality_test {
    use crate::lin_alg::matrix::Matrix;
    #[test] 
    fn iter_test() {
        let m1 = matrix!(integer;
            [[1, 0, 1],
             [2, 1, -1],
             [2, -3, -1]]
        );
        let v: Vec<_> = m1.row_iter(1).map(|&x| x).collect();
        assert_eq!(v, vec![2,1,-1]);

        let v: Vec<_> = m1.col_iter(1).map(|&x| x).collect();
        assert_eq!(v, vec![0,1,-3]);
    }

    #[test] 
    fn ref_test() {
        let mut m1 = matrix!(integer;
            [[1, 0, 1],
             [2, 1, -1],
             [2, -3, -1]]
        );
        m1[1] = 3;
    }

    #[test]
    fn zero_test() {
        let m = matrix!(integer;
            [[1, 0, 0],
             [0, 0, 0]]
        );

        assert!(!m.is_col_zero(0));
        assert!(m.is_col_zero(1));
        assert!(m.is_col_zero(2));
        assert!(!m.is_row_zero(0));
        assert!(m.is_row_zero(1));
    }
}


impl<T: PID> Add for Matrix<T> 
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        if (self.n_row, self.n_col) != (rhs.n_row, rhs.n_col) { panic!(); }

        let mut sum = Matrix::new(self.n_row, self.n_col);

        self.iter()
            .zip(rhs.iter())
            .enumerate()
            .for_each( |(k, (&x, &y))|
                sum.write( k / sum.n_col, k % sum.n_col, x + y )
        );

        sum
    }
}

impl<T: PID> Mul for Matrix<T>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.n_col != rhs.n_row { panic!(); }

        let mut product = Matrix::<T>::new(self.n_row, rhs.n_col);

        (0..product.n_row * product.n_col).for_each( |n| {
            let i = n / product.n_col;
            let j = n % product.n_col;
            let r_iter = self.row_iter(i);
            let c_iter = rhs.col_iter(j);
            let val = r_iter
                .zip(c_iter)
                .map(|(&a, &b)| a*b)
                .fold(T::zero(), |acc, x| acc + x);
            product.write(i, j, val);
        } );

        product
    }
}


#[cfg(test)]
mod arithmetic_test {
    use crate::lin_alg::matrix::Matrix;
    #[test]
    fn addition_test() {
        let m1 = matrix!(integer;
            [[1, 1],
             [2, -1]]
        );
        let m2 = matrix!(integer;
            [[1, 0],
             [2, -1]]
        );
        let m = m1 + m2;
        let ans = matrix!(integer;
            [[2, 1],
             [4, -2]]
        );
        assert_eq!(m ,ans);


        use crate::commutative::rational::*;

        let m1 = matrix!(rational;
            [[-3, 1],
             [5, 0]]
        );
        let m2 = matrix!(
            [[0.over(1), 1.over(5)],
             [1.over(1), 3.over(5)]]
        );

        assert_eq!(m1*m2, Matrix::identity(2));
    }

    #[test]
    fn multiplication_test() {
        let m1 = matrix!(integer;
            [[1, 0, 1],
             [2, 1, -1]]
        );
        let m2 = matrix!(integer;
            [[1, 0],
             [1, 0],
             [2, -1]]
        );
        let m = m1 * m2;
        let ans = matrix!(integer;
            [[3, -1],
             [1, 1]]
        );
        assert_eq!(m ,ans);


        use crate::commutative::rational::*;

        let m1 = matrix!(rational;
            [[-3, 1],
             [5, 0]]
        );
        let m2 = matrix!(
            [[0.over(1), 1.over(5)],
             [1.over(1), 3.over(5)]]
        );

        assert_eq!(m1*m2, Matrix::identity(2));
    }
}


#[allow(unused)]
impl <T:Field> Matrix<T> {
    pub fn rank_as_linear_map(&self) -> usize {
        let mut clone = self.clone();
        let shorter_side = std::cmp::min(self.n_row, self.n_col);
        let mut skipped = 0;
        for k in 0..shorter_side {
            if clone.get(k,k) == T::zero() {
                let result = clone.row_iter(k).enumerate().skip(k+1).find(|(_, &x)| x != T::zero() );
                match result {
                    Some((j, _)) => clone.col_operation(k, T::one(), j, T::one()),
                    None => (),
                }
            }
            if clone.get(k,k) == T::zero() {
                let result = clone.col_iter(k).enumerate().skip(k+1).find(|(_, &x)| x != T::zero() );
                match result {
                    Some((i, _)) => clone.row_operation(k, T::one(), i, T::one()),
                    None => {
                        skipped += 1;
                        continue;
                    },
                }
            }

            // Now we have the nonzero element in (k, k)
            let coeff = T::one() / clone.get(k,k);
            clone.row_operation(k, coeff, k, T::zero());

            for l in k+1..self.n_row {
                let minus_one = T::zero() - T::one();
                let a_lk = clone.get(l,k);
                if a_lk == T::zero() { continue; } 
                clone.row_operation(l, T::one(), k, minus_one * a_lk);
            }
        };

        shorter_side - skipped
    }
}


#[allow(unused)]
impl <T: PID> Matrix<T> {
    pub fn smith_normal_form(mut self) -> (Self, Self, Self) {
        let mut r_op = Matrix::zero(self.n_row, self.n_row);
        let mut c_op = Matrix::zero(self.n_col, self.n_col);

        let (mut i, mut j) = (0, 0);

        while( i < self.n_row && j < self.n_col) {
            // choosing a pivot
            j = match (j..self.n_col).find(|&l| !self.is_col_zero(l)) {
                Some(l) => l,
                None => break,
            };
            if self.get(i, j) == T::zero() {
                let k = (i..self.n_row).find(|&k| self.get(k, j) != T::zero() ).unwrap();
                self.swap_rows(i, k);
            }

            debug_assert!(self.get(i, j)!=T::zero());

            // Improving the pivot
            // after this while loop, the entry in (i,j) must devides all the elements below and right
            while (i+1..self.n_row).any(|k| !self.get(i, j).divides(self.get(k, j))) || (j+1..self.n_col).any(|l| !self.get(i, j).divides(self.get(i, l))) {
                for k in i+1..self.n_row {
                    if !self.get(i, j).divides(self.get(k, j)) {
                        let a = self.get(i, j);
                        let b = self.get(k, j);
                        let (gcd, x, y) = bezout_identity(a, b);

                        // row operations
                        let i_th_row = self.clone_row(i).scalar_mul(x) + self.clone_row(k).scalar_mul(y);
                        let k_th_row = self.clone_row(i).scalar_mul(T::negative_one()* b.euclid_div(gcd).0) + self.clone_row(k).scalar_mul(a.euclid_div(gcd).0);

                        self.write_row(i, i_th_row);
                        self.write_row(k, k_th_row);

                        debug_assert!(self.get(i, j)==gcd);
                    }
                }

                for l in j+1..self.n_col {
                    if !self.get(i, j).divides(self.get(i, l)) {
                        let a = self.get(i, j);
                        let b = self.get(i, l);
                        let (gcd, x, y) = bezout_identity(a, b);

                        // col operations
                        let piv_th_col = self.clone_col(j).scalar_mul(x) + self.clone_col(l).scalar_mul(y);
                        let l_th_col = self.clone_col(j).scalar_mul(T::negative_one()* b.euclid_div(gcd).0) + self.clone_col(l).scalar_mul(a.euclid_div(gcd).0);

                        self.write_col(j, piv_th_col);
                        self.write_col(l, l_th_col);

                        debug_assert!(self.get(i, j)==gcd);
                    }
                }
            }

            // Eliminating entries
            for k in i+1..self.n_row {
                if self.get(k, j) == T::zero() { continue; }

                let ratio = self.get(k, j).euclid_div(self.get(i, j)).0;
                self.row_operation(k, T::negative_one(), i, ratio);
            }

            for l in j+1..self.n_col {
                if self.get(i, l) == T::zero() { continue; }

                let ratio = self.get(i, l).euclid_div(self.get(i, j)).0;
                self.col_operation(l, T::negative_one(), j, ratio);
            }

            i += 1; j+= 1;
        };


        // Bringing all the pivots to the diagonal
        for r in 0..std::cmp::min(self.n_row, self.n_col) {
            // if the entry at (r, r) is zero then we do nothing 
            if self.get(r, r) != T::zero() { continue; }

            // otherwise, we fetch the nonzero element and swap it with the entry at (r, r)
            let s = (r+1..self.n_col).find(|&s| self.get(r, s) != T::zero() );
            let s = match s {
                Some(s) => s,
                None => break,
            };

            self.write(r, r, self.get(r,s));
            self.write(r, s, T::zero());
        }


        // Finally, force this diagonal to satisfy a_{i, i} | a_{i+1, i+1} for each i
        let rank = match (0..std::cmp::min(self.n_row, self.n_col)).find(|&r| self.get(r,r) == T::zero()) {
            Some(r) => r,
            None => std::cmp::min(self.n_row, self.n_col),
        };

        for r in 0..rank {
            for s in r+1..rank {
                if self.get(r, r).divides(self.get(s, s)) {
                    continue;
                }

                // This col operation makes the entry at (k, i) nonzero 
                self.col_operation(r, T::one(), s, T::one());

                let a = self.get(r, r);
                let b = self.get(s, r);
                let (gcd, x, y) = bezout_identity(a, b);

                // row operations
                let r_th_row = self.clone_row(r).scalar_mul(x) + self.clone_row(s).scalar_mul(y);
                let s_th_row = self.clone_row(r).scalar_mul(T::negative_one()* b.euclid_div(gcd).0) + self.clone_row(s).scalar_mul(a.euclid_div(gcd).0);

                self.write_row(r, r_th_row);
                self.write_row(s, s_th_row);

                // A column operation to eliminate the entry at (i, k)
                self.col_operation(s, T::one(), r, T::negative_one() * b.euclid_div(gcd).0 * y);

                debug_assert!(self.get(r, r)==gcd);
                debug_assert!(self.get(r, s)==T::zero());
                debug_assert!(self.get(s, r)==T::zero());
            }
        }
        debug_assert!((0..rank).all( |i| (i+1..rank).all(|k| self.get(i, i).divides(self.get(k, k)))), " Diagonal is not sorted: {:?}", self);

        // Fianlly, we make all the elements "non_negative"
        for r in 0..rank {
            if self.get(r,r).is_nonnegative() { continue; }
            let val = T::negative_one() * self.get(r,r);
            self.write(r,r, val);
        }

        (r_op, self, c_op)
    }
}



#[cfg(test)]
mod operation_test {
    use crate::lin_alg::matrix::Matrix;
    use crate::commutative::rational::*;
    #[test]
    fn row_col_operation_test() {
        let mut m = matrix!(integer;
            [[1, 0, 1],
             [2, 3, 1],
             [-1, 2, -1]]
        );
        m.row_operation(0, 3, 2, 1);
        let m1 = matrix!(integer;
            [[2, 2, 2],
             [2, 3, 1],
             [-1, 2, -1]]
        );
        assert_eq!(m, m1);
        m.col_operation(1, 1, 0, -2);
        let m2 = matrix!(integer;
            [[2, -2, 2],
             [2, -1, 1],
             [-1, 4, -1]]
        );
        assert_eq!(m, m2);
    }

    #[test]
    fn rank_as_linear_map_test() {
        let m = Matrix::<Rational>::identity(100);
        let r = m.rank_as_linear_map();
        assert_eq!(r, 100);

        let m = matrix!(rational;
            [[2, -2, 2],
             [2, -1, 1],
             [-1, 4, -1]]
        );
        let r = m.rank_as_linear_map();
        assert_eq!(r, 3);

        let m = matrix!(rational;
            [[-1, 0, -1],
             [1, -1, 0],
             [0, 1, 1]]
        );
        let r = m.rank_as_linear_map();
        assert_eq!(r, 2);

        let m = matrix!(rational;
            [[0, 0, -1],
             [0, -2, 0],
             [0, 0, 0]]
        );
        let r = m.rank_as_linear_map();
        assert_eq!(r, 2);

        let m = matrix!(rational;
            [[0, 0, -1],
             [0, 1, 0],
             [1, 0, 0]]
        );
        let r = m.rank_as_linear_map();
        assert_eq!(r, 3);

        let m = matrix!(rational;
            [[0, 0, 0],
             [0, 0, 0],
             [0, 0, -1]]
        );
        let r = m.rank_as_linear_map();
        assert_eq!(r, 1);

        let m = Matrix::<Rational>::identity(100);
        let r = m.rank_as_linear_map();
        assert_eq!(r, 100);

        let m = Matrix::<Rational>::zero(50, 60);
        let r = m.rank_as_linear_map();
        assert_eq!(r, 0);

        let mut m = Matrix::<Rational>::zero(60, 50);
        let r = m.rank_as_linear_map();
        assert_eq!(r, 0);

        m.write(25, 25, 1.over(25));
        let r = m.rank_as_linear_map();
        assert_eq!(r, 1);
    }

    #[test]
    fn smith_normal_form_test() {
        // test 1
        let m = matrix!(integer;
            [[2, 4, 4],
             [-6, 6, 12],
             [10, 4, 16]]
        );
        let (_, m, _) = m.smith_normal_form();
        let ans = matrix!(integer;
            [[2, 0, 0],
             [0, 2, 0],
             [0, 0, 156]]
        );
        assert_eq!(m, ans);

        // test 2
        let (n_row, n_col) = (2, 3);
        let mut  m = Matrix::<i64>::zero(n_row,n_col);
        m.write(0, 0, 1);
        for i in 0..n_row {
            for j in 0..n_col {
                println!("testing i={}, j={}", i, j);
                let mut n = Matrix::<i64>::zero(n_row,n_col);
                n.write(i, j, 1);
                let (_, n, _) = n.smith_normal_form();
                assert_eq!(m, n);
            }
        }

        // test 3
        let (n_row, n_col) = (3, 2);
        let mut  m = Matrix::<i64>::zero(n_row,n_col);
        m.write(0, 0, 1);
        for i in 0..n_row {
            for j in 0..n_col {
                println!("testing i={}, j={}", i, j);
                let mut n = Matrix::<i64>::zero(n_row,n_col);
                n.write(i, j, 1);
                let (_, n, _) = n.smith_normal_form();
                assert_eq!(m, n);
            }
        }

        // test 4
        let m = Matrix::<i64>::zero(1,1);
        let (_, n, _) = m.clone().smith_normal_form();
        assert_eq!(m, n);
        let m = Matrix::<i64>::identity(1);
        let (_, n, _) = m.clone().smith_normal_form();
        assert_eq!(m, n);
        let m = Matrix::<i64>::identity(15);
        let (_, n, _) = m.clone().smith_normal_form();
        assert_eq!(m, n);
        let m = Matrix::<i64>::zero(5, 8);
        let (_, n, _) = m.clone().smith_normal_form();
        assert_eq!(m, n);
        let m = Matrix::<i64>::zero(8, 5);
        let (_, n, _) = m.clone().smith_normal_form();
        assert_eq!(m, n);

        // test 5
        let m = matrix!(integer;
            [[2, 4, 4],
            [-6, 6, 12]]
        );
        let (_, m, _) = m.smith_normal_form();
        let ans = matrix!(integer;
            [[2, 0, 0],
            [0, 6, 0]]
        );
        assert_eq!(m, ans);

        let m = matrix!(integer;
            [[6, -6],
             [4, 6],
             [4, 12]]
        );
        let (_, m, _) = m.smith_normal_form();
        let ans = matrix!(integer;
            [[2, 0],
             [0, 6],
             [0, 0]]
        );
        assert_eq!(m, ans);
    }
}
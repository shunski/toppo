use std::{rc::Rc, fmt::Display};
use crate::lin_alg::Matrix;
use super::{SparseNode, SparseMatrix};
use crate::commutative::PID;

impl<T: PID> SparseNode<T> {
    #[inline]
    fn row_idx(&self) -> usize {
        match self {
            Self::Entry((i,_),_, _) => *i,
            _ => panic!("Failed to unwwrap: The node is not an entry.")
        }
    }

    #[inline]
    fn col_idx(&self) -> usize {
        match self {
            Self::Entry((_,j),_, _) => *j,
            _ => panic!("Failed to unwwrap: The node is not an entry.")
        }
    }

    #[inline]
    fn start_ptr(self) -> Rc<SparseNode<T>> {
        match self {
            Self::Start(_, r) => r,
            _ => panic!("Failed to unwwrap: The node is not a start node.")
        }
    }

    #[inline]
    fn start_node_idx(&self) -> usize {
        match self {
            Self::Start(i, _) => *i,
            _ => panic!("Failed to unwwrap: The node is not a start node.")
        }
    }
    
    #[inline]
    fn start_ref(&self) -> &SparseNode<T> {
        match self {
            Self::Start(_, r) => &*r,
            _ => panic!("Failed to unwwrap: The node is not a start node."),
        }
    }
}


impl<T: PID, F: Fn((usize, usize)) -> Option<T>> From<((usize, usize), F)> for SparseMatrix<T> {
    fn from((size, f): ((usize, usize), F)) -> Self {
        let mut m = SparseMatrix{rows: Vec::new(), cols: Vec::new(), size, transposed: false};

        for i in (0..size.0).rev() {
            let mut curr_node =  Rc::new(SparseNode::End);
            for j in (0..size.1).rev() {
                let val = if let Some(val) = f((i,j)) {
                    val
                } else {
                    continue;
                };

                // check if there is already an element in the 'j'th column.
                let new_node_ptr = if let Some(start_node) = m.cols.iter_mut().find(|node| node.start_node_idx()==j) {
                    // create a dummy node and swap it with the 'start_node', so that 'start_node' won't have to move out.
                    let mut dummy = SparseNode::End;
                    std::mem::swap(&mut dummy, start_node);

                    let new_node = SparseNode::Entry((i,j), val, (curr_node, dummy.start_ptr()));
                    let new_node_ptr = Rc::new(new_node);
                    
                    // Note that '&mut start_node' is a mut reference to the element of 'm.cols', with the dummy value in it.
                    // We replace this dummy value with the pointer to the 'new_node'
                    *start_node = SparseNode::Start(j, Rc::clone(&new_node_ptr));
                    new_node_ptr
                } else { // this means that it is the first entry in 'j'th col
                    let new_node = SparseNode::Entry((i,j), val, (curr_node, Rc::new(SparseNode::End)));
                    let new_node_ref = Rc::new(new_node);
                    m.cols.push( SparseNode::Start(j, Rc::clone(&new_node_ref)));
                    new_node_ref
                };

                // update the 'curr_node' and go to the next iteration
                curr_node = new_node_ptr;
            }

            // if 'curr_node' is not 'End', then it means that there was some none-zero entry in the row,
            // so in this case add the row. Otherwise, do nothing.
            if *curr_node != SparseNode::End {
                m.rows.push(SparseNode::Start(i, curr_node));
            }
        }

        // sort the 'rows' and 'cols' so that we can run the binary search on them later
        m.rows.sort_by(|x, y| x.start_node_idx().cmp( &y.start_node_idx() ) );
        m.cols.sort_by(|x, y| x.start_node_idx().cmp( &y.start_node_idx() ) );


        // if 'rows' and 'cols' are empty, then add some elements to 'rows' so that the iter() does not panic.
        if m.rows.is_empty() {
            m.rows.push(SparseNode::End);
            m.cols.push(SparseNode::End);
        }

        m
    }
}



pub struct Iter<'a, T: PID> {
    data: &'a SparseNode<T>,
    is_done: bool,
    dir: bool,
}

impl<T: PID> SparseMatrix<T> {
    pub fn row_iter(&self, i: usize) -> Iter<'_, T> {
        // let row = self.rows.iter().find(|f| i == f.start_node_idx() );
        let row = self.rows.binary_search_by(|f| f.start_node_idx().cmp(&i) );
        match row {
            Ok(x) => Iter{
                data: &self.rows[x].start_ref(),
                is_done: false,
                dir: !self.transposed,
            },
            Err(_) => Iter{
                data: &self.rows[0], // Recall that 'rows' is not empty.
                is_done: true,
                dir: !self.transposed
            }
        }
    }

    pub fn col_iter(&self, i: usize) -> Iter<'_, T> {
        let col = self.cols.binary_search_by(|f| f.start_node_idx().cmp(&i) );
        match col {
            Ok(x) => Iter {
                data: &self.cols[x].start_ref(),
                is_done: false,
                dir: self.transposed,
            },
            Err(_) => Iter {
                data: &self.rows[0], // Recall that 'rows' is not empty.
                is_done: true,
                dir: self.transposed,
            }
        }
    }

    #[inline]
    pub fn size(&self) -> (usize, usize) {
        self.size
    }

    #[inline]
    pub fn transpose(mut self) -> Self {
        self.transposed = !self.transposed;
        std::mem::swap(&mut self.rows, &mut self.cols);
        std::mem::swap(&mut self.size.0, &mut self.size.1);
        self
    }
}

impl<'a, T: PID> Iterator for Iter<'a, T> {
    type Item = (usize, &'a T);
    fn next(&mut self) -> Option<Self::Item> {
        if self.is_done {
            return None
        }

        match self.data {
            SparseNode::Start(..) => panic!(),
            SparseNode::Entry(idx, val, r) => {
                let r = if self.dir {&r.0} else {&r.1};
                self.data = &*r;
                let idx = if self.dir {idx.1} else {idx.0};
                Some((idx, val))
            },
            SparseNode::End => {
                self.is_done = true;
                return None;
            }
        }
    }
}

impl<T: PID> std::cmp::PartialEq<Matrix<T>> for SparseMatrix<T> {
    fn eq(&self, other: &Matrix<T>) -> bool {
        if self.size() != other.size() {
            return false;
        }

        for i in 0..self.size().0 {
            let mut j = 0;
            for (k, val) in self.row_iter(i) {
                if (j..k).any(|l| other[(i,l)]!=T::zero()) || &other[(i,k)] != val {
                    return false;
                }
                j = k+1;
            }
            if (j..self.size().1).any(|k| other[(i,k)]!=T::zero() ) {
                return false;
            }
        }

        return true
    }
}

impl<T: PID> std::ops::Mul for SparseMatrix<T> {
    type Output = Matrix<T>;
    fn mul(self, rhs: Self) -> Self::Output {
        if self.size().1 != rhs.size().0 {
            panic!("size not compatible");
        }

        let mut m = Matrix::zero(self.size().0, rhs.size().1);

        for i in 0..self.size().0 {
            for j in 0..rhs.size().1 {
                let mut sum = T::zero();
                let mut it1 = self.row_iter(i);
                let mut it2 = rhs.col_iter(j);
                let (mut k2, mut v2) = if let Some(next) = it2.next() {
                    next
                } else {
                    continue;
                };

                while let Some((mut k1, mut v1)) = it1.next() {
                    if k1 == k2 {
                        sum += v1.clone() * v2.clone();
                    } else if k1 > k2 {
                        std::mem::swap(&mut it1, &mut it2);
                        std::mem::swap(&mut k1, &mut k2);
                        std::mem::swap(&mut v1, &mut v2);
                    }
                }

                m[(i,j)] = sum;
            }
        }
        m
    }
}


impl<T: PID, S: PID> std::ops::Mul<Matrix<S>> for SparseMatrix<T> 
    where T: std::ops::Mul<S>, <T as std::ops::Mul<S>>::Output: PID,
{
    type Output = Matrix<<T as std::ops::Mul<S>>::Output>;
    fn mul(self, rhs: Matrix<S>) -> Self::Output {
        if self.size().1 != rhs.size().0 {
            panic!("size not compatible");
        }

        let mut m = Matrix::zero(self.size().0, rhs.size().1);

        for i in 0..self.size().0 {
            for j in 0..rhs.size().1 {
                let mut it = self.row_iter(i);
                while let Some((k,val)) = it.next() {
                    m[(i,j)] += val.clone() * rhs[(k,j)].clone();
                }
            }
        }
        m
    }
}

impl<T: PID, S: PID> std::ops::Mul<SparseMatrix<T>> for Matrix<S> 
    where T: std::ops::Mul<S>, <T as std::ops::Mul<S>>::Output: PID,
{
    type Output = Matrix<<T as std::ops::Mul<S>>::Output>;
    fn mul(self, rhs: SparseMatrix<T>) -> Self::Output {
        if self.size().1 != rhs.size().0 {
            panic!("size not compatible");
        }

        (rhs.transpose() * self.transpose()).transpose()
    }
}


impl<T: PID + Display> Display for SparseMatrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.size().0 {
            for (j, val) in self.row_iter(i) {
                write!(f, "({i}, {j}): {val}\n")?;
            }
        }
        write!(f, "")
    }
}

#[cfg(test)]
mod test {
    use crate::lin_alg::Matrix;

    use super::SparseMatrix;
    #[test]
    fn init_test() {
        let m = SparseMatrix::from(((2, 2), |(i,j)| if (i, j) == (0,1) {Some((2) as i64)} else {None} ));
        // row iter test
        assert_eq!(m.row_iter(1).next(), None);
        let mut iter = m.row_iter(0);
        assert_eq!(iter.next(), Some((1, &2))); 
        assert_eq!(iter.next(), None);

        // col iter test
        assert_eq!(m.col_iter(0).next(), None);
        let mut iter = m.col_iter(1);
        assert_eq!(iter.next(), Some((0, &2))); 
        assert_eq!(iter.next(), None);

        // transepose 'm' and run the same test
        let m = m.transpose();
        // col iter test
        assert_eq!(m.col_iter(1).next(), None);
        let mut iter = m.col_iter(0);
        assert_eq!(iter.next(), Some((1, &2))); 
        assert_eq!(iter.next(), None);

        // row iter test
        assert_eq!(m.row_iter(0).next(), None);
        let mut iter = m.row_iter(1);
        assert_eq!(iter.next(), Some((0, &2))); 
        assert_eq!(iter.next(), None);
        

        // Another test with big matrix
        let n = 100;
        let f = |(i,j)| {
            if i*n+j % 71 == 0 {
                Some((i*n+j) as i64)
            } else {
                None
            }
        };
        let m = SparseMatrix::from(((n, n), f));

        let mut m_ans = Matrix::zero(n, n);
        for i in 0..n {
            for j in 0..n {
                m_ans[(i,j)] = match f((i,j)) {
                    Some(x) => x,
                    None => 0,
                }
            }
        }

        assert!( m == m_ans )
    }

    #[test]
    fn eq() {
        let n = 100;
        let f = |(i,j)| {
            if i*n+j % 71 == 0 {
                Some((i*n+j) as i64)
            } else {
                None
            }
        };
        let m = SparseMatrix::from(((n, n), f));

        let mut m_ans = Matrix::zero(n, n);
        for i in 0..n {
            for j in 0..n {
                m_ans[(i,j)] = match f((i,j)) {
                    Some(x) => x,
                    None => 0,
                }
            }
        }

        assert!( m == m_ans );

        m_ans[(0,71)] = 1;
        assert!( m != m_ans );
    }

    #[test]
    fn mul() {
        // test 1
        let n = 2;
        let f = |(i,j)| {
            Some((i*n+j) as i64)
        };
        let m_sparse = SparseMatrix::from(((n, n), f));

        let mut m_dense = Matrix::zero(n, n);
        for i in 0..n {
            for j in 0..n {
                m_dense[(i,j)] = match f((i,j)) {
                    Some(x) => x,
                    None => 0,
                }
            }
        }

        assert!(m_sparse == m_dense);

        let prod = m_sparse.clone() * m_sparse;
        let ans = m_dense.clone() * m_dense;

        assert_eq!(prod, ans);

        // test 2
        let n = 100;
        let f = |(i,j)| {
            if i*n+j % 11 == 0 {
                Some((i*n+j) as i64 / 103 - 59)
            } else {
                None
            }
        };
        let m_sparse = SparseMatrix::from(((n, n), f));

        let mut m_dense = Matrix::zero(n, n);
        for i in 0..n {
            for j in 0..n {
                m_dense[(i,j)] = match f((i,j)) {
                    Some(x) => x,
                    None => 0,
                }
            }
        }

        println!("m_dense = {}", m_dense);

        let p1 = m_sparse.clone() * m_sparse.clone();
        let p2 = m_dense.clone() * m_dense.clone();
        let p4 = m_sparse.clone() * m_dense.clone();
        let p3 = m_dense * m_sparse;

        assert_eq!(p1, p2);
        assert_eq!(p2, p3);
        assert_eq!(p3, p4);
        
    }
}
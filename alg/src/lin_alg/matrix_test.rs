#[cfg(test)]
mod init_test {
    use crate::lin_alg::Matrix;
    use crate::matrix;
    #[test]
    fn mat_init() {
        let m1 = Matrix::<i64>::zero(5,4);
        let m2 = matrix!(i64; 
            [[0,0,0,0],
             [0,0,0,0],
             [0,0,0,0],
             [0,0,0,0],
             [0,0,0,0]]
        );
        assert_eq!(m1, m2);

        let m1 = Matrix::<i64>::identity(5);
        let m2 = matrix!(i64; 
            [[1,0,0,0,0],
             [0,1,0,0,0],
             [0,0,1,0,0],
             [0,0,0,1,0],
             [0,0,0,0,1]]
        );
        assert_eq!(m1, m2);
    }

    #[test]
    fn sub_mut() {
        let mut m = matrix!(i64;
            [[1, 0, 1],
             [2, 1, -1],
             [2, -3, -1]]
        );

        let mut sub = m.sub_mut(1.., 1..);
        sub[(0,1)] = 10;
        let sub_sub = sub.sub_mut(.., 0);
        assert_eq!( m, 
            matrix!(i64;
                [[1,  0,  1],
                [2,  1, 10],
                [2, -3, -1]]
            )
        );
        sub_sub.write(
            matrix!(i64;
                [[4, 5]]
            ).transpose()
        );
        assert_eq!( m, 
            matrix!(i64;
                [[1, 0,  1],
                [2, 4, 10],
                [2, 5, -1]]
            )
        );
    }
}

#[cfg(test)]
mod basic_functionality_test {
    use crate::lin_alg::Matrix;
    use crate::matrix;
    #[test] 
    fn index() {
        let mut m = matrix!(i64;
            [[1, 0, 1],
             [2, 1, -1],
             [2, -3, -1]]
        );
        m[(0,1)] = 3;
        assert_eq!(m[(0,1)], 3);
        assert_eq!(m.sub(0, ..)[(0,1)], 3);
        assert_eq!(m.sub(0, 1..)[(0,0)], 3);
    }

    #[test]
    fn transpose() {
        let m = matrix!(i64;
            [[1, 0, 1],
             [2, 1, -1],
             [2, -3, -1]]
        );
        let m_t = matrix!(i64;
            [[1, 2, 2],
             [0, 1, -3],
             [1, -1, -1]]
        );
        assert_eq!(m.clone().transpose(), m_t);

        let m_sub_t = matrix!(i64;
            [[2, 2],
             [1, -3]]
        );
        assert_eq!(m.sub(1.., ..2).transpose().as_matrix(), m_sub_t);
    }

    fn det2<T: crate::commutative::PID>(m: &Matrix<T>) -> T {
        m[(0,0)]*m[(1,1)] - m[(0,1)]*m[(1,0)]
    }

    fn det3<T: crate::commutative::PID>(m: &Matrix<T>) -> T {
        m[(0,0)]*m[(1,1)]*m[(2,2)]
         + m[(0,1)]*m[(1,2)]*m[(2,0)]
         + m[(0,2)]*m[(1,0)]*m[(2,1)]
         - m[(0,2)]*m[(1,1)]*m[(2,0)]
         - m[(0,1)]*m[(1,0)]*m[(2,2)]
         - m[(0,0)]*m[(1,2)]*m[(2,1)]
    }

    #[test]
    fn det() {
        // testing for 2 by 2 matrices
        for _ in 0..10 {
            let m = Matrix::<f64>::random(2, 2);
            assert!( (m.clone().det() - det2(&m)).abs() / det2(&m).abs() < 0.0001, "m={m:?}" );
        }

        // testing for 3 by 3 matrices
        for _ in 0..10 {
            let m = Matrix::<f64>::random(3, 3);
            assert!( (m.clone().det() - det3(&m)).abs() / det3(&m).abs() < 0.0001, "m={m:?}" );
        }
    }

    #[test]
    fn row_operation() {
        let mut m = matrix!(i64;
            [[1,  0,  1],
             [2,  1, -1],
             [2, -3, -1]]
        );

        m.swap_rows(0, 2);
        assert_eq!(m,
            matrix!(i64;
                [[2, -3, -1],
                 [2,  1, -1],
                 [1,  0,  1]]
            )
        );

        m.sub_mut(.., 0..2).swap_rows(0, 2);
        assert_eq!(m,
            matrix!(i64;
                [[1,  0, -1],
                 [2,  1, -1],
                 [2, -3, 1]]
            )
        );

        m.row_operation(2, 0, 3);
        assert_eq!(m,
            matrix!(i64;
                [[1,  0, -1],
                 [2,  1, -1],
                 [5, -3, -2]]
            )
        );

        m.sub_mut(.., 1..).row_operation(1, 2, -2);
        assert_eq!(m,
            matrix!(i64;
                [[1,  0, -1],
                 [2,  7,  3],
                 [5, -3, -2]]
            )
        );

        m.sub_mut(0, ..2).write( m.sub(0, ..2)*3 + m.sub(1, ..2)*(-2) );
        assert_eq!(m,
            matrix!(i64;
                [[-1,  -14, -1],
                 [2,  7,  3],
                 [5, -3, -2]]
            )
        );
    }

    #[test]
    fn col_operation() {
        let mut m = matrix!(i64;
            [[1,  0,  1],
             [2,  1, -1],
             [2, -3,  1]]
        );

        m.swap_cols(0, 2);
        assert_eq!(m,
            matrix!(i64;
                [[ 1,  0, 1],
                 [-1,  1, 2],
                 [ 1, -3, 2]]
            )
        );

        m.sub_mut(0..2, ..).swap_cols(0, 2);
        assert_eq!(m,
            matrix!(i64;
                [[ 1,  0,  1],
                 [ 2,  1, -1],
                 [ 1, -3,  2]]
            )
        );

        m.col_operation(2, 0, 3);
        assert_eq!(m,
            matrix!(i64;
                [[ 1,  0,  4],
                 [ 2,  1,  5],
                 [ 1, -3,  5]]
            )
        );

        m.sub_mut(1.., ..).col_operation(1, 2, -2);
        assert_eq!(m,
            matrix!(i64;
                [[ 1,   0,  4],
                 [ 2,  -9,  5],
                 [ 1, -13,  5]]
            )
        );
    }
}

#[cfg(test)]
mod arithmetic_test {
    use crate::lin_alg::Matrix;
    use crate::{matrix, rational};
    use crate::commutative::Rational;
    #[test]
    fn addition_test() {
        let m1 = matrix!(i64;
            [[1, 1],
             [2, -1]]
        );
        let m2 = matrix!(i64;
            [[1, 0],
             [2, -1]]
        );
        let m = m1 + m2;
        let ans = matrix!(i64;
            [[2, 1],
             [4, -2]]
        );
        assert_eq!(m ,ans);

        let m1 = matrix!(Rational;
            [[rational!(-3), rational!(1)],
             [rational!(5), rational!(0)]]
        );
        let m2 = matrix!(Rational;
            [[rational!(0; 1), rational!(1; 5)],
             [rational!(1; 1), rational!(3; 5)]]
        );

        assert_eq!(m1*m2, Matrix::identity(2));
    }

    #[test]
    fn multiplication_test() {
        let m1 = matrix!(i64;
            [[1, 0, 1],
             [2, 1, -1]]
        );
        let m2 = matrix!(i64;
            [[1, 0],
             [1, 0],
             [2, -1]]
        );
        let m = m1 * m2;
        let ans = matrix!(i64;
            [[3, -1],
             [1, 1]]
        );
        assert_eq!(m ,ans);

        let m1 = matrix!(Rational;
            [[rational!(-3), rational!(1)],
             [rational!(5), rational!(0)]]
        );
        let m2 = matrix!(Rational;
            [[rational!(0), rational!(1; 5)],
             [rational!(1), rational!(3; 5)]]
        );

        assert_eq!(m1*m2, Matrix::identity(2));
    }
}

#[cfg(test)]
mod operation_test {
    use crate::lin_alg::Matrix;
    use crate::matrix;
    use crate::commutative::Rational;
    #[test]
    fn row_col_operation_test() {
        let mut m = matrix!(i64;
            [[1, 0, 1],
             [2, 3, 1],
             [-1, 2, -1]]
        );
        m.row_operation(0, 2, 2);
        let m1 = matrix!(i64;
            [[-1, 4, -1],
             [2, 3, 1],
             [-1, 2, -1]]
        );
        assert_eq!(m, m1);
        m.col_operation(1, 0, -2);
        let m2 = matrix!(i64;
            [[-1, 6, -1],
             [2, -1, 1],
             [-1, 4, -1]]
        );
        assert_eq!(m, m2);
    }

    #[test]
    fn rank_as_linear_map_test() {
        use crate::rational;

        let m = Matrix::<Rational>::identity(100);
        let r = m.rank_as_linear_map();
        assert_eq!(r, 100);

        let m = matrix!(Rational;
            [[rational!(2), rational!(-2), rational!(2)],
             [rational!(2), rational!(-1), rational!(1)],
             [rational!(-1), rational!(4), rational!(-1)]]
        );
        let r = m.rank_as_linear_map();
        assert_eq!(r, 3);

        let m = matrix!(Rational;
            [[rational!(-1), rational!(0), rational!(-1)],
             [rational!(1), rational!(-1), rational!(0)],
             [rational!(0), rational!(1), rational!(1)]]
        );
        let r = m.rank_as_linear_map();
        assert_eq!(r, 2);

        let m = matrix!(Rational;
            [[rational!(0), rational!(0), rational!(-1)],
             [rational!(0),rational!(-2), rational!(0)],
             [rational!(0), rational!(0), rational!(0)]]
        );
        let r = m.rank_as_linear_map();
        assert_eq!(r, 2);

        let m = matrix!(Rational;
            [[rational!(0), rational!(0), rational!(-1)],
             [rational!(0), rational!(1), rational!(0)],
             [rational!(1), rational!(0), rational!(0)]]
        );
        let r = m.rank_as_linear_map();
        assert_eq!(r, 3);

        let m = matrix!(Rational;
            [[rational!(0), rational!(0), rational!(0)],
             [rational!(0), rational!(0), rational!(0)],
             [rational!(0), rational!(0), rational!(-1)]]
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

        m[(25, 25)] = rational!(1; 25);
        let r = m.rank_as_linear_map();
        assert_eq!(r, 1);
    }

    macro_rules! assert_eq_up_to_mul_by_unit {
        ($m: expr, $n: expr) => {
            assert_eq!( $m.size, $n.size);

            for i in 0..$m.size.0 {
                for j in 0..$m.size.1 {
                    if i==j {
                        assert!( $m[(i,j)] == $n[(i,j)] || $m[(i,j)] == -$n[(i,j)], "left={:?}right = {:?}", $m, $n );
                    } else {
                        assert!( $m[(i,j)] == 0 && $n[(i,j)] == 0, "left={:?}right = {:?}", $m, $n );
                    }
                }
            }
        }
    }

    #[test]
    fn smith_normal_form_test() {
        // test 1
        let m = matrix!(i64;
            [[2, 4, 4],
             [-6, 6, 12],
             [10, 4, 16]]
        );
        let (r_inv, r, n, c, c_inv) = m.clone().smith_normal_form();
        let ans = matrix!(i64;
            [[2, 0, 0],
             [0, 2, 0],
             [0, 0, 156]]
        );
        assert_eq_up_to_mul_by_unit!(n, ans);
        assert_eq!(r.clone() * r_inv.clone(), Matrix::identity(3));
        assert_eq!(c.clone() * c_inv.clone(), Matrix::identity(3));
        assert_eq!(r.clone() * m.clone() * c.clone(), n.clone());
        assert_eq!(r_inv.clone() * n * c_inv.clone(), m.clone());
        assert_eq!(r_inv * r, Matrix::identity(m.size.0));
        assert_eq!(c_inv * c, Matrix::identity(m.size.1));

        // test 2
        let size = (2, 3);
        let mut  m = Matrix::<i64>::zero(size.0,size.1);
        m[(0, 0)] = 1;
        for i in 0..size.0 {
            for j in 0..size.1 {
                let mut n = Matrix::<i64>::zero(size.0,size.1);
                n[(i, j)] = 1;
                let (_,_, n, _,_) = n.smith_normal_form();
                assert_eq!(m, n);
            }
        }

        // test 3
        let size = (3, 2);
        let mut  m = Matrix::<i64>::zero(size.0,size.1);
        m[(0, 0)] = 1;
        for i in 0..size.0 {
            for j in 0..size.1 {
                let mut n = Matrix::<i64>::zero(size.0,size.1);
                n[(i, j)] = 1;
                let (_,_, n, _,_) = n.smith_normal_form();
                assert_eq!(m, n);
            }
        }

        // test 4
        let m = Matrix::<i64>::zero(1,1);
        let (_,_, n, _,_) = m.clone().smith_normal_form();
        assert_eq_up_to_mul_by_unit!(m, n);
        let m = Matrix::<i64>::identity(1);
        let (_,_, n, _,_) = m.clone().smith_normal_form();
        assert_eq!(m, n);
        let m = Matrix::<i64>::identity(15);
        let (_,_, n, _,_) = m.clone().smith_normal_form();
        assert_eq!(m, n);
        let m = Matrix::<i64>::zero(5, 8);
        let (_,_, n, _,_) = m.clone().smith_normal_form();
        assert_eq!(m, n);
        let m = Matrix::<i64>::zero(8, 5);
        let (_,_, n, _,_) = m.clone().smith_normal_form();
        assert_eq!(m, n);

        // test 5
        let m = matrix!(i64;
            [[2, 4, 4],
            [-6, 6, 12]]
        );
        let (_,_, m, _,_) = m.smith_normal_form();
        let ans = matrix!(i64;
            [[2, 0, 0],
            [0, 6, 0]]
        );
        assert_eq_up_to_mul_by_unit!(m, ans);

        let m = matrix!(i64;
            [[6, -6],
             [4, 6],
             [4, 12]]
        );
        let (_,_, m, _,_) = m.smith_normal_form();
        let ans = matrix!(i64;
            [[2, 0],
             [0, 6],
             [0, 0]]
        );
        assert_eq_up_to_mul_by_unit!(m, ans);


        // test 6
        let m = matrix!(i64;
            [[ 2, 4,  9, 5, 4, 24],
             [-6, 6, 2,  0, 16,  2],
             [6, 3, 3, 9, 18, 9],
             [ 7, 8, 7,  5, 9,  1]]
        );
        let (r_inv,r, n, c, c_inv) = m.clone().smith_normal_form();
        assert_eq!(r.clone() * r_inv.clone(), Matrix::identity(4));
        assert_eq!(c.clone() * c_inv.clone(), Matrix::identity(6));
        assert_eq!(r.clone() * m.clone() * c.clone(), n.clone());
        assert_eq!(r_inv.clone() * n * c_inv.clone(), m.clone());
        assert_eq!(r_inv * r, Matrix::identity(m.size.0));
        assert_eq!(c_inv * c, Matrix::identity(m.size.1));
    }
}


#[cfg(test)]
mod numerical_functionality_test {
    use crate::matrix;
    use crate::lin_alg::Matrix;
    #[test]
    fn norms() {
        let v = matrix!(f64;
            [[-2.0],
             [8.0],
             [-13.0],
             [7.0]]
        );
        assert!( (v.one_norm() - 30.0).abs() < 0.0001, "v.one_norm()={}", v.one_norm() );
        assert!( (v.two_norm() - 16.91153).abs() < 0.0001, "v.two_norm()={}", v.two_norm() );
        assert!( (v.inf_norm() - 13.0).abs() < 0.0001, "v.inf_norm()={}", v.inf_norm() );
        assert!( (v.two_norm() - v.frobenius_norm()).abs() < 0.0001, "v.frobenius_norm()={}", v.frobenius_norm() );


        let m = matrix!(f64;
            [[ 4.0,  9.0, -5.0,  4.0,  24.0],
             [ 6.0,  2.0, 0.0,  16.0,  -2.0],
             [ 6.0, -3.0, 9.0, -18.0,  9.0],
             [ 7.0, -7.0, 5.0,  -9.0,  1.0]]
        );

        assert!( (m.one_norm() - (4.0+16.0+18.0+9.0)).abs() < 0.0001, "m.one_norm()={}", m.one_norm() );
        assert!( (m.inf_norm() - (4.0+9.0+5.0+4.0+24.0)).abs() < 0.0001, "m.inf_norm()={}", m.inf_norm() );
    }

    #[test]
    fn random_orthogonal() {
        for i in 1..10 {
            let m = Matrix::random_orthogonal(i);
            let out = m.sub(..,..).transpose() * m;
            assert!( (out - Matrix::identity(i)).frobenius_norm() / (i as f64).sqrt() < 0.0001 );
        }
    }

    #[test]
    fn solve() {
        // test 1
        let a = matrix!{ f64;
            [[1.0,  2.0],
             [4.0, -2.0]]
        };
        let b = matrix!{ f64;
            [[2.0], [3.0]]
        };

        let x = a.solve(b);
        let ans = matrix!{ f64;
            [[1.0], [0.5]]
        };

        assert!((x.clone() - ans.clone()).two_norm() / ans.two_norm() * 100.0 < 0.0001, 
            "The error is larger than 0.0001%. The output is {}",
            x
        );

        // test 2
        for i in 1..20 {
            let a = Matrix::random(i, i);
            let b = Matrix::random(i, 1);

            let x = a.clone().solve(b.clone());

            assert!( (a * x.clone() - b.clone()).two_norm() / b.two_norm() * 100.0 < 0.0001,
                "The error is larger than 0.0001%. The output is {}",
                x
            );
        }
    }

    #[test]
    fn householder_vec() {
        let u = Matrix::random(10, 1);
        let (v, b) = u.clone().householder_vec();
        let u = u.sub(.., ..);
        let v = v.sub(.., ..);
        let mut y = u - v * (v.transpose() * u) * b;
        y[(0,0)] -= u.as_matrix().two_norm();
        assert!( y.two_norm() < 0.0001, "y has to be zero, but y={:?}", y );

        let u = Matrix::random(10, 1);
        let (v, b) = u.clone().householder_vec();
        let v = v.sub(.., ..);
        let p = Matrix::identity(10) - v * v.transpose()* b;
        assert!( (p.clone().det() + 1.0)<0.00001, "p must have det 1, but det(p)={:?}", p.det() );
    }

    #[test]
    fn hessenberg_form() {
        // test 1
        let mut m = matrix!(f64;
            [[-2.0, -4.0, 2.0],
             [-2.0,  1.0, 2.0],
             [ 4.0,  2.0, 5.0]]
        );
        m.hessenberg_form();
        let ans = matrix!(f64;
            [[-2.0, 8.0/5_f64.sqrt(), -6.0/5_f64.sqrt()],
            [2.0*5_f64.sqrt(),  13.0/5.0, 14.0/5_f64],
            [ 0.0,  14.0/5_f64, 17.0/5_f64]]
        );
        assert!( (m.clone() - ans.clone()).frobenius_norm() < 0.00001, "m has to be {ans:.3} but m={m:.3}");

        // test 2
        let n = 10;
        let mut m = Matrix::random(n, n);
        let k = m.clone();
        m.hessenberg_form();
        let mut zero = 0.0;
        for i in 2..n {
            for j in 0..i-1 {
                zero += m[(i,j)].abs();
            }
        }
        assert!(zero<0.0001, "m is not Hessenberg; m={}", m);
        assert!((m.frobenius_norm() - k.frobenius_norm()).abs() < 0.00001, "The norm does not match.");
    }

    #[test]
    fn francis_step() {
        let n = 10;
        let mut m = Matrix::random(n, n);
        for i in 2..n {
            for j in 0..i-2 {
                m[(i,j)] = 0.0;
            }
        }
        let k = m.clone();
        m.francis_step();
        assert!((m.frobenius_norm() - k.frobenius_norm()).abs() < 0.00001, 
            "The norm does not match. left={}, right={}",
            m.frobenius_norm(), k.frobenius_norm()
        );
    }


    #[test]
    fn spectrum() {
        let m = matrix!(f64;
            [[3.0, 1.0],
             [2.0, 2.0]]
        );
        let spectrum = m.spectrum();
        let ans = matrix!(f64; [[1.0], [4.0]]);
        assert!( (spectrum.clone()-ans).two_norm() < 0.00001, 
            "spectrum must be [1, 4].transpose(), but spectrum = {:.2}", 
            spectrum
        );

        let m = matrix!(f64;
            [[-0.2, -0.4, 0.2],
             [-0.2,  0.1, 0.2],
             [ 0.4,  0.2, 0.5]]
        );
        let spectrum = m.spectrum();
        let ans = matrix!(f64; [[-0.5], [0.3], [0.6]]);
        let error = (spectrum.clone()-ans.clone()).two_norm() / ans.clone().two_norm() * 100.0;
        assert!( error < 0.0001, 
            "spectrum must be [-0.5, 0.3, 0.6].transpose(), but spectrum = {:.4} and the error = {}%.", 
            spectrum,
            error
        );
    }

    #[test] 
    fn spectrum_with_invariant_space() {
        let m = matrix!(f64;
            [[-0.2, -0.4, 0.2],
             [-0.2,  0.1, 0.2],
             [ 0.4,  0.2, 0.5]]
        );

        let (spectrum, eigen_vecs) = m.clone().spectrum_with_invariant_space();

        let n = spectrum.size.0;
        let m = m.sub(..,..);

        let eigen_vecs = eigen_vecs.sub(..,..);
        let spectrum_on_diagonal = {
            let mut a = Matrix::zero(n,n);
            (0..n).for_each(|i| a[(i,i)] = spectrum[(i,0)] );
            a
        };
        let error = {
            let a = m * eigen_vecs; 
            let b = eigen_vecs * spectrum_on_diagonal;
            let b = b.sub(..,..);
            (a-b).frobenius_norm() / b.frobenius_norm() * 100.0
        };

        assert!( error < 0.0001,
            "Failed on Eigenvectors. Error is {:.5}%, which is greater than 0.0001%.",
            error
        );

        // testing with larger matries
        let n = 10;
        let diagonal = {
            let mut m = Matrix::zero(n,n);
            (0..n).for_each( |i| { m[(i,i)] = (i+1) as f64; });
            m
        };
        let diagonal = diagonal.sub(..,..);

        let shuffle = Matrix::random_orthogonal(n);
        let shuffle = shuffle.sub(.., ..);
        let sample = shuffle.transpose() * diagonal * shuffle;

        let (spectrum, eigen_vecs) = sample.clone().spectrum_with_invariant_space();

        // evaluation of spectrum
        let error = {
            let mut clone = spectrum.clone();
            (0..n).for_each( |i| clone[(i,0)] -= diagonal[(i,i)] );
            clone.two_norm() / diagonal.frobenius_norm() * 100.0
        };

        assert!( error < 0.0001,
            "error is {:.5}%, which is greater than 0.0001%. the spectrum = {:.3}.",
            error, spectrum
        );

        // evaluation of eigenvectors
        let spectrum_on_diagonal = {
            let mut a = Matrix::zero(n,n);
            (0..n).for_each(|i| a[(i,i)] = spectrum[(i,0)] );
            a
        };

        let error = {
            let a = sample * eigen_vecs.clone();
            let b = eigen_vecs * spectrum_on_diagonal;
            (a-b.clone()).frobenius_norm() / b.frobenius_norm() * 100.0
        };

        assert!( error < 0.0001,
            "Failed on Eigenvectors. Error is {:.5}%, which is greater than 0.0001%.",
            error
        );
    }
}
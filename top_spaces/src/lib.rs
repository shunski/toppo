use algebra::commutative::PID;
use algebra::module::matrix::Matrix;

pub trait Space<Coeff: PID + std::fmt::Debug> {
    type Output;
    fn homology(&self, _: Coeff) -> Self::Output;
    fn get_boundary_map (&self) -> Vec<Matrix<Coeff>>;
}

pub mod complex;
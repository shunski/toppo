use algebra::commutative::PID;
use algebra::module::matrix::Matrix;

pub trait Space<Coeff: PID + std::fmt::Debug> {
    type Output;
    fn homology(&self) -> Self::Output;
    fn get_boundary_map (&self) -> Vec<Matrix<Coeff>>;
}

#[macro_export]
macro_rules! H {
    ($space:expr; $coeff:ty) => {
        <crate::simplicial::SimplicialCplx<String> as Space<$coeff>>::homology(&$space)
    };
}

pub mod simplicial;

pub mod cubical;

pub mod graph;
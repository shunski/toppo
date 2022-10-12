use algebra::commutative::{PID, Field};
use algebra::module::matrix::Matrix;

pub trait Space<Coeff: Field> {
    fn homology(&self) -> Vec<usize> {
        let boundary_maps: Vec<Matrix<Coeff>> = self.get_boundary_map();
        if boundary_maps.is_empty() { return Vec::new(); }
        let im_and_ker: Vec<_> = boundary_maps
            .into_iter()
            .map(|boundary_map| {
                let im: usize = boundary_map.rank_as_linear_map();
                let ker: usize = boundary_map.size().1 - im;
                (im, ker)
            })
            .collect();

        // println!("im_and_ker: {:?}", im_and_ker);
        let mut homology = Vec::new();
        let mut kernel = im_and_ker[0].1;
        for &(im, ker) in im_and_ker.iter().skip(1) {
            homology.push(kernel - im);
            kernel = ker;
        }
        homology.push(im_and_ker.last().unwrap().1);
        homology
    }
    fn get_boundary_map (&self) -> Vec<Matrix<Coeff>>;
}

#[macro_export]
macro_rules! H {
    ($space:expr; $coeff:ty) => {
        <_ as Space<$coeff>>::homology(&$space)
    };
}

pub mod simplicial;

pub mod graph;

pub mod cubical;
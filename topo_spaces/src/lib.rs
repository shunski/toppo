use algebra::lin_alg::Module;
use algebra::lin_alg::matrix::Matrix;

// pub trait HomologyWithFieldCoefficients <Coeff: Field> {
//     fn homology(&self) -> Vec<usize> {
//         let boundary_maps: Vec<Matrix<Coeff>> = self.get_boundary_map();
//         if boundary_maps.is_empty() { return Vec::new(); }
//         let im_and_ker: Vec<_> = boundary_maps
//             .into_iter()
//             .map(|boundary_map| {
//                 let im: usize = boundary_map.rank_as_linear_map();
//                 let ker: usize = boundary_map.size().1 - im;
//                 (im, ker)
//             })
//             .collect();

//         // println!("im_and_ker: {:?}", im_and_ker);
//         let mut homology = Vec::new();
//         let mut kernel = im_and_ker[0].1;
//         for &(im, ker) in im_and_ker.iter().skip(1) {
//             homology.push(kernel - im);
//             kernel = ker;
//         }
//         homology.push(im_and_ker.last().unwrap().1);
//         homology
//     }
//     fn get_boundary_map (&self) -> Vec<Matrix<Coeff>>;
// }


#[allow(unused)]
pub struct IntegralHomology {
    homology: Vec<Module<i64>>,
}


impl IntegralHomology {
    pub fn zero() -> IntegralHomology {
        IntegralHomology { homology: Vec::new() }
    }

    pub fn new(cplx: &impl Complex, max_dim: Option<usize>) -> IntegralHomology {
        let boundary_maps: Vec<Matrix<i64>> = cplx.boundary_map();
        if boundary_maps.is_empty() { return Self::zero(); }

        let mut rank_of_higher_map = 0;
        let mut elem_divisors_of_higher_map: Vec<(i64, usize)> = Vec::new();
        let mut homology: Vec<_> = boundary_maps
            .into_iter()
            .rev()
            .map(|boundary_map| {
                // get the rank of the map, dimention of tha chain ccomplex
                let (_, boundary_map, _) = boundary_map.smith_normal_form();
                let shorter_side = std::cmp::min(boundary_map.size().0, boundary_map.size().1);
                let rank_of_lower_map = match (0..shorter_side).find(|&i| boundary_map.get(i,i)==0){
                    Some(i) => i,
                    None => shorter_side,
                };
                let dim = boundary_map.size().1;

                // compute the homology
                let homoology_at_this_dim = Module::new(dim-rank_of_lower_map-rank_of_higher_map, elem_divisors_of_higher_map.clone());
                
                // update the variable for the computation of the lower homology
                elem_divisors_of_higher_map = Vec::new();
                (0..rank_of_lower_map).for_each(|i| {
                    let divisor = boundary_map.get(i,i);
                    if let Some((div, mul)) = elem_divisors_of_higher_map.last_mut()  {
                        if *div == divisor {
                            *mul+=1;
                        } else {
                            elem_divisors_of_higher_map.push((boundary_map.get(i,i), 1)); 
                        }
                    } else {
                        elem_divisors_of_higher_map.push((boundary_map.get(i,i), 1)); 
                    }
                });
                rank_of_higher_map = rank_of_lower_map;

                // return
                homoology_at_this_dim
            })
            .collect();

        // reverse
        homology.reverse();

        // cut off the unnecessary higher dimension
        if let Some(max_dim) = max_dim {
            while homology.len()-1 > max_dim { homology.pop(); }
        }
        while let Some(module) = homology.last() {
            if module.is_zero() {
                homology.pop();
            } else {
                break;
            }
        }

        // return
        IntegralHomology { homology: homology }
    }
}

pub trait Complex {
    fn boundary_map (&self) -> Vec<Matrix<i64>>;
}

use std::fmt::Display;
impl Display for IntegralHomology {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, module) in self.homology.iter().enumerate() {
            write!(f, "dim {}: {},\n", i, module)?;
        }
        write!(f, "")
    }
}


#[macro_export]
macro_rules! H {
    ($complex:expr) => {
        IntegralHomology::new(&$complex, None).to_string()
    };
    ($complex:expr, max_dim=) => {
        IntegralHomology::new(&$complex, Some(max_dim)).to_string()
    };
}

pub mod simplicial;

pub mod graph;

pub mod cubical;
use alg::commutative::Field;
use alg::lin_alg::{Module, FormalSum, Matrix};

pub trait HomologyWithFieldCoefficients <Coeff: Field> {
    fn homology(&self) -> Vec<usize> {
        let boundary_maps: Vec<Matrix<Coeff>> = self.get_boundary_map();
        if boundary_maps.is_empty() { return Vec::new(); }
        let im_and_ker: Vec<_> = boundary_maps
            .into_iter()
            .map(|boundary_map| {
                let im: usize = boundary_map.rank_as_linear_map();
                let ker: usize = boundary_map.size.1 - im;
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


#[allow(unused)]
pub struct IntegralHomology<T: Complex> {
    homology: Vec<Module<i128>>,
    cycles: Vec<Vec<FormalSum<T::Cell>>>,
    boundaries: Vec<Vec<FormalSum<T::Cell>>>,
}


impl<T: Complex> IntegralHomology<T> 
    where i128: std::ops::Mul<T::Cell, Output = FormalSum<T::Cell>> + std::ops::Mul<FormalSum<T::Cell>, Output = FormalSum<T::Cell>>
{
    pub fn new(cplx: &T, max_dim: Option<usize>) -> Self {
        let boundary_maps = cplx.boundary_map();
        assert!( !boundary_maps.is_empty() );

        let mut rank_of_higher_map = 0;
        let mut elem_divisors_of_higher_map: Vec<(i128, usize)> = Vec::new();
        let mut row_op_of_higher_map = Matrix::<i128>::identity( boundary_maps.last().unwrap().0.size.1 );
        let mut output: Vec<_> = boundary_maps
            .into_iter()
            .rev()
            .map(|(boundary_map, basis)| {
                // get the rank of the map, dimention of tha chain complex
                // println!("boundary map: {:?}", boundary_map );
                // println!("boundary map: {:?}", boundary_map );
                // println!("boundary map: {:?}", boundary_map.iter().filter(|&&i| i!=0).count() );
                // println!("boundary map: {:?}", boundary_map );
                let (r_op_inv, _, boundary_map, c_op, _) = boundary_map.smith_normal_form();

                let shorter_side = std::cmp::min(boundary_map.size.0, boundary_map.size.1);
                let rank_of_lower_map = match (0..shorter_side).find(|&i| boundary_map[( i,i )]==0) {
                    Some(i) => i,
                    None => shorter_side,
                };
                let dim = boundary_map.size.1;

                // compute the homology
                let homology_at_this_dim = Module::new(dim-rank_of_lower_map-rank_of_higher_map, elem_divisors_of_higher_map.clone());
                // compute cycles
                let mut cycles_at_this_dim = c_op.transpose() * basis.clone();
                cycles_at_this_dim  = cycles_at_this_dim.into_iter().enumerate().filter(|&(i, _)| i >= rank_of_lower_map ).map(|(_, x)| x ).collect();
                // compute boundaries
                let mut boundaries_at_this_dim = row_op_of_higher_map.clone().transpose() * basis;
                boundaries_at_this_dim = boundaries_at_this_dim
                    .into_iter()
                    .enumerate()
                    .filter(|&(i, _)| i < rank_of_higher_map )
                    .map(|(mut i, x)| { 
                        let mut idx = 0;
                        while elem_divisors_of_higher_map[idx].1 <= i {
                            i -= elem_divisors_of_higher_map[idx].1;
                            idx+=1;
                        }
                        elem_divisors_of_higher_map[idx].0 * x
                    } )
                    .collect();
                
                // update the variables for the computation of the lower homology
                elem_divisors_of_higher_map = Vec::new();
                (0..rank_of_lower_map).for_each(|i| {
                    let divisor = boundary_map[( i,i )];
                    if let Some((div, mul)) = elem_divisors_of_higher_map.last_mut()  {
                        if *div == divisor {
                            *mul+=1;
                        } else {
                            elem_divisors_of_higher_map.push((boundary_map[( i,i )], 1)); 
                        }
                    } else {
                        elem_divisors_of_higher_map.push((boundary_map[( i,i )], 1)); 
                    }
                });
                rank_of_higher_map = rank_of_lower_map;
                row_op_of_higher_map = r_op_inv;

                // return
                ( homology_at_this_dim, cycles_at_this_dim, boundaries_at_this_dim )
            })
            .collect();

        // reverse
        output.reverse();

        // cut off the unnecessary higher dimension
        if let Some(max_dim) = max_dim {
            while output.len()-1 > max_dim { output.pop(); }
        }
        while let Some(( module, _, _  )) = output.last() {
            if module.is_zero() {
                output.pop();
            } else {
                break;
            }
        }

        // unpack
        let (mut homology, mut cycles, mut boundaries) = (Vec::new(), Vec::new(), Vec::new());
        output.into_iter().for_each(|(h, c, b)| {
            homology.push( h );
            cycles.push( c );
            boundaries.push( b );
        } );

        // return
        IntegralHomology { homology: homology, cycles: cycles, boundaries: boundaries }
    }
}


impl<T: Complex> IntegralHomology<T> 
    where i128: std::ops::Mul<T::Cell, Output = FormalSum<T::Cell>>, T::Cell: Display
{
    pub fn get_description_of_cycles(&self) -> String {
        let mut out = String::new();
        for (i, i_cycles) in self.cycles.iter().enumerate() {
            out += &std::fmt::format( format_args!("{}-cycles: \n", i) );
            i_cycles.iter().for_each( |i_cycle| 
                out += &std::fmt::format( format_args!("\t{}, \n", i_cycle) )
            );
        };
        out
    }
    

    pub fn get_description_of_boundaries(&self) -> String {
        let mut out = String::new();
        for (i, i_boundaries) in self.boundaries.iter().enumerate() {
            out += &std::fmt::format( format_args!("{}-boundaries: \n", i) );
            i_boundaries.iter().for_each( |i_boundary| 
                out += &std::fmt::format( format_args!("\t{}, \n", i_boundary) )
            );
        };
        out
    }
}


pub trait Complex {
    type Cell: Clone + PartialEq;
    fn boundary_map (&self) -> Vec< (Matrix<i128>, Vec<Self::Cell>) >;
}

use std::fmt::Display;
impl<T: Complex> Display for IntegralHomology<T> {
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
        IntegralHomology::new(&$complex, None)
    };
    ($complex:expr, max_dim=) => {
        IntegralHomology::new(&$complex, Some(max_dim))
    };
}

pub mod simplicial;

pub mod graph;

pub mod cubical;
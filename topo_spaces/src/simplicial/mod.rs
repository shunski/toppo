mod mod_impl;
mod mod_test;

pub trait VertexLabel: std::clone::Clone + std::cmp::Eq {}
  
macro_rules! vertex_label_impl {
    ($($t:ty)*) => ($(
        impl VertexLabel for $t {}
    )*)
}

vertex_label_impl!{ String }

#[allow(unused)]
#[derive(Clone, PartialEq, Debug)] 
pub struct Simplex <T: VertexLabel> 
{
    pub dim: usize,
    pub vertices_labels: Vec<T>,
}

use alg::lin_alg::FormalSum;
use alg::formal_sum_impl;
formal_sum_impl!{ Simplex<String> }

#[macro_export]
macro_rules! simplex {
    ( $( $vertex:expr ), * ) => {
        {
            let mut v = Vec::new();
            $(
                v.push($vertex.to_string());
            )*
            Simplex::from(v)
        }
    };
}

#[allow(unused)]
#[derive(Clone, Debug, Eq)] 
struct SimplicialCell {
    vertices: Vec<u64>,
    size: usize, // size of the complex
    dim: Option<usize>, // dim=Null iff the simplicial complex in empty
}

pub struct SimplicialCplx <T: VertexLabel> 
{
    size: usize, // number of vertices
    maximal_cells: Vec<SimplicialCell>,
    vertices_labeling: Vec<T>,
}

#[macro_export]
macro_rules! simplicial_cplx {
    ( $( {$( $vertex:expr ), *}), * ) => {
        {
            let mut simplices = Vec::new();
            $(
                let mut vertices = Vec::new();
                $(
                    vertices.push($vertex.to_string());
                )*
                simplices.push( Simplex::from(vertices) );
            )*
            let mut cplx = SimplicialCplx::new();
            simplices.into_iter().for_each(|simplex| cplx.attach(simplex));
            cplx
        }
    };
}

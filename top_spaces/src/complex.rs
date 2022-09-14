use crate::Space;
use util::bitwise::count_ones;
use algebra::commutative::Field;

pub trait VertexLabel: std::clone::Clone + std::cmp::Eq {}

macro_rules! vertex_label_impl {
    ($($t:ty)*) => ($(
        impl VertexLabel for $t {}
    )*)
}

vertex_label_impl!{ String }

#[allow(unused)]
#[derive(Clone)] 
pub struct Simplex <T: VertexLabel> 
{
    dim: usize,
    vertices_labels: Vec<T>,
}

#[allow(unused)]
impl<T: VertexLabel> Simplex<T> 
{
    pub fn new_from(labels: Vec<T>) -> Simplex<T> {
        Simplex{
            dim: labels.len() - 1,
            vertices_labels: Vec::from_iter(labels.into_iter()),
        }
    }
}


#[allow(unused)]
#[derive(Clone, Debug, Eq)] 
struct SimplicialCell {
    vertices: Vec<u64>,
    size: usize, // size of the complex
    dim: usize,
}


#[allow(unused)]
impl SimplicialCell {
    fn new(s: usize, mut raw_vertices: Vec<usize>) -> SimplicialCell {
        raw_vertices.sort();
        {
            // "SimplicialCell" is used only in the implementaion of "SimplicialCplx,"
            // and hence the interface can be pretty restrictive.
            // These are the prohibited inputs
            if s == 0 || raw_vertices.is_empty() { panic!();}

            if s < raw_vertices[0] { panic!(); } 

            if raw_vertices.iter().enumerate().any(|(i, &val)| {
                if i+1 < raw_vertices.len() {
                    raw_vertices[i] == raw_vertices[i+1]
                } else {
                    false
                }
            }) { panic!(); }
        }

        let mut cell = SimplicialCell {
            vertices: vec!(0; s / 64 + 1),
            size: s,
            dim: 0,
        };

        raw_vertices.iter().for_each(|&x| cell.cone(x));

        cell.dim = SimplicialCell::dim(&cell.vertices);

        cell
    }

    #[inline]
    fn cone(&mut self, vertex: usize) {
        if vertex >= self.size {
            let new_size = vertex+1;
            self.resize(new_size);
        }
        self.vertices[vertex/64] |= 1 << (vertex % 64);
        self.dim += 1;
    }

    fn join(&mut self, cell: &SimplicialCell) {
        self.vertices = self.vertices
            .iter()
            .zip(cell.vertices.iter()).map(|(&x, &y)| x | y)
            .collect();
            
        if cell.vertices.len() > self.vertices.len() { 
            (self.vertices.len()..cell.vertices.len()).
            for_each(|i| self.vertices.push(cell.vertices[i]));
        }
    }

    fn intersection(&mut self, cell: &SimplicialCell) -> SimplicialCell {
        let vertices = self.vertices
            .iter()
            .zip(cell.vertices.iter()).map(|(&x, &y)| x & y)
            .collect();
        
        SimplicialCell {
            dim: SimplicialCell::dim(&vertices),
            vertices: vertices,
            size: self.size, 
        }
    }

    fn dim(vertices: &Vec<u64>) -> usize {
        vertices.iter()
        .map( |&vertex| {
            count_ones(vertex)
        } )
        .fold(0, |acc, x| acc+x)
        -1
    }

    #[inline] 
    fn resize(&mut self, new_size: usize) {
        while new_size > 64 * self.vertices.len() {
            self.vertices.push(0);
        }
        self.size = new_size;
    }

    fn is_connected_with(&self, cell: &SimplicialCell) -> bool {
        self.vertices.iter().zip(cell.vertices.iter()).any(|(&x, &y)| x & y != 0)
    }

    fn contains(&self, cell: &SimplicialCell) -> bool {
        self.vertices.iter().zip(cell.vertices.iter()).all(|(&x, &y)| x | y == x)
    }

    fn get_face(&self) -> Vec<SimplicialCell> {
        let n = self.dim-1;

        let mut n_sub_cells = Vec::new();
        n_sub_cells.reserve(n+2);
        let mut idx = 0;
        for j in 0..n+2 {
            let mut sub_cell = self.clone();
            sub_cell.dim -= 1;
            idx = (idx..).find(|&new_idx| (sub_cell.vertices[new_idx/64]>>(new_idx%64)) % 2 == 1).unwrap();
            sub_cell.vertices[idx/64] &= !(1 << (idx%64));
            n_sub_cells.push(sub_cell);
            idx+=1;
        };
        n_sub_cells
    }
}

impl PartialEq for SimplicialCell {
    fn eq(&self, rhs: &Self) -> bool {
        if self.vertices.len() == rhs.vertices.len() {
            return self.vertices == rhs.vertices
        } else {
            self.vertices.iter().zip(rhs.vertices.iter()).all(|(&x, &y)| x==y)
        }
    }
}

use std::cmp::Ordering;
impl Ord for SimplicialCell {
    fn cmp(&self, rhs: &Self) -> Ordering {
        if self.dim > rhs.dim {
            return Ordering::Greater
        } else if self.dim < rhs.dim {
            return Ordering::Less
        }

        // if self.dim == rhs.dim
        let result = self.vertices.iter()
            .zip(rhs.vertices.iter())
            .map(|(x, y)| x.cmp(y))
            .find(|&cmp| cmp != Ordering::Equal);
        match result {
            Some(cmp) => cmp,
            None => Ordering::Equal,
        }
    }
}

impl PartialOrd for SimplicialCell {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}

#[cfg(test)]
mod simplicial_cell_test {
    use super::SimplicialCell;

    #[test]
    fn init_test() {
        let vertices = vec![10,6,2,3];
        let mut cell = SimplicialCell::new(30, vertices);

        let mut answer = 0;
        answer |= 1 << 2;
        answer |= 1 << 3;
        answer |= 1 << 6;
        answer |= 1 << 10;
        assert_eq!(cell.vertices[0], answer);
        assert_eq!(cell.size, 30);
        assert_eq!(cell.dim, 3);

        cell.cone(150);
        let mut answer2 = 0;
        answer2 |= 1 << 150 % 64;
        assert_eq!(cell.vertices[0], answer);
        assert_eq!(cell.vertices[1], 0);
        assert_eq!(cell.vertices[2], answer2);
        assert_eq!(cell.size, 151);
        assert_eq!(cell.dim, 4);
    }

    #[test]
    fn functionality_test(){
        let vertices = vec![1,2,3,10];
        let mut cell1 = SimplicialCell::new(100, vertices);

        let vertices = vec![9, 49, 99];
        let mut cell2 = SimplicialCell::new(100, vertices);

        // testing is_connected_with()
        assert_eq!(cell1.is_connected_with(&cell2), false);
        cell1.cone(99);
        assert_eq!(cell1.is_connected_with(&cell2), true);

        // testing contains()
        assert_eq!(cell1.contains(&cell2), false);
        cell1.cone(9);
        cell1.cone(49);
        assert_eq!(cell1.contains(&cell2), true);
        cell2.cone(33);
        assert_eq!(cell1.contains(&cell2), false);
    }

    #[test]
    fn sub_cell_test() {
        let vertices = vec![9, 49, 99, 200, 3000];
        let cell1 = SimplicialCell::new(4000, vertices);
        let sub_cells = cell1.get_face();
        let ans = vec![
            SimplicialCell::new(4000, vec![49, 99, 200, 3000]),
            SimplicialCell::new(4000, vec![9, 99, 200, 3000]),
            SimplicialCell::new(4000, vec![9, 49, 200, 3000]),
            SimplicialCell::new(4000, vec![9, 49, 99, 3000]),
            SimplicialCell::new(4000, vec![9, 49, 99, 200])
        ];
        assert_eq!(sub_cells, ans);
    }
}


pub struct SimplicialCplx <T: VertexLabel> 
{
    size: usize, // number of vertices
    maximal_cells: Vec<SimplicialCell>,
    vertices_labeling: Vec<T>,
}

#[allow(unused)]
impl<T: VertexLabel> SimplicialCplx<T> 
{
    pub fn new() -> SimplicialCplx<T> {
        SimplicialCplx{
            size: 0,
            maximal_cells: Vec::new(),
            vertices_labeling: Vec::new(),
        }
    }

    pub fn new_from_zero_skeleton(vertices_labels: Vec<T>) -> SimplicialCplx<T> {
        SimplicialCplx {
            size: vertices_labels.len(),
            maximal_cells: 
                (0..vertices_labels.len())
                .map(|x| SimplicialCell::new(1, vec![x]))
                .collect(),
            vertices_labeling: vertices_labels,
        }
    }

    pub fn new_from(simplex: Simplex<T>) -> SimplicialCplx<T> 
    {
        let mut cplx = SimplicialCplx::new();
        cplx.attach(simplex);
        cplx
    }

    pub fn attach(&mut self, simplex: Simplex<T>){
        // add vertices to the complex if necessary
        for v_label in simplex.vertices_labels.iter() {
            if !self.vertices_labeling.contains(&v_label) {
                self.size += 1;
                self.vertices_labeling.push(v_label.clone());
                self.maximal_cells.push(SimplicialCell::new(self.size, vec![self.size-1]));
            }
        }

        // build the ordered simplex using Simplex data
        let mut new_cell = SimplicialCell::new( self.size, 
            simplex
            .vertices_labels
            .into_iter()
            .map( |label1| {
                self.vertices_labeling.iter().position(|label2|label1 == *label2).unwrap()
            })
            .collect()
        );

        // if the new simplex is contained in the complex then do nothing
        if self.maximal_cells.iter()
            .any(|maximal_cell| maximal_cell.contains(&new_cell)) {
            return;
        }

        // Otherwise we will add this simplex to the complex
        // First remove all the simplices that are contained in the new simplex, 
        // then modify the connected components (some components may be glued by adding),
        // and finally add the new simplex.
        self.maximal_cells.retain(|maximal_cell| !new_cell.contains(maximal_cell));

        self.maximal_cells.push(new_cell);
    }

    pub fn get_labeled_zero_skeleton(&self) -> &Vec<T> {
        return &self.vertices_labeling;
    }

    pub fn dim(&self) -> usize {
        self.maximal_cells
            .iter()
            .map(|cell| cell.dim)
            .max()
            .unwrap()
    }

    pub fn is_connected(&self) -> bool {
        // if the complex is empty, then it is vacuously connected.
        if self.maximal_cells.is_empty() {
            return true;
        }

        // Otherwise, we take the naive approach to detect connectedness
        let mut join = self.maximal_cells[0].clone();
        let mut prev_dim = join.dim;
        loop {
            // iterate each cell and join it if connected.
            let cells_connected_to_join: Vec<_> = self.maximal_cells.iter()
                .filter(|&cell| cell.is_connected_with(&join))
                .collect();
            cells_connected_to_join.iter().for_each(|cell| { join.join(cell); });

            // if the join eventually gets the join of all points, then the complex is connected
            if join.dim == self.size-1 {
                return true;
            }

            // if the join comes here and if the join stoped growing, then the complex is not connected.
            if join.dim == prev_dim {
                return false;
            }

            prev_dim = join.dim;
        }
    }

    fn get_i_cells(&self, i: usize, i_plus_1_cells: Option<&Vec<SimplicialCell>>) -> Vec<SimplicialCell> {
        let mut i_cells:Vec<_> = self.maximal_cells
            .iter()
            .filter(|&maximal_cell|maximal_cell.dim == i)
            .map(|maximal_cell| maximal_cell.clone())
            .collect();
        if i_plus_1_cells == None {return i_cells;}
        let i_plus_1_cells = i_plus_1_cells.unwrap();
        i_cells.reserve((i+2)*i_plus_1_cells.len());
        i_plus_1_cells.iter().for_each( |i_plus_1_cell| {
            let mut i_sub_cells = i_plus_1_cell.get_face();
            i_cells.append(&mut i_sub_cells);
        });

        i_cells.sort();
        i_cells.dedup();

        i_cells
    }
}


#[cfg(test)]
mod simplicial_cplx_tests {
    #[test]
    fn object_instantiation() {
        use crate::complex::*;

        let mut cplx = SimplicialCplx::new_from_zero_skeleton(vec!["v0".to_string(), "v1".to_string()]);
        assert_eq!(cplx.is_connected(), false);

        let mut labels = vec!["v1".to_string(), "v2".to_string(), "v3".to_string()];
        let simplex = Simplex::new_from(labels.clone());
        assert_eq!(simplex.dim, 2);

        cplx.attach(simplex.clone());
        labels.push("v0".to_string());
        assert_eq!(
            cplx.get_labeled_zero_skeleton().len(), 
            labels.len()
        );

        cplx.attach(simplex.clone());
        assert_eq!(
            cplx.get_labeled_zero_skeleton().len(), 
            labels.len()
        );

        let mut labels = vec!["v0".to_string(), "v1".to_string()];
        let simplex = Simplex::new_from(labels.clone());
        let cplx = SimplicialCplx::new_from(simplex);
        assert_eq!(cplx.dim(), 1);

        labels.pop();
        labels.push("v2".to_string());
        labels.push("v3".to_string());
        labels.push("v4".to_string());
        let simplex = Simplex::new_from(labels.clone());
        let cplx = SimplicialCplx::new_from(simplex);
        assert_eq!(cplx.dim(), 3);
        assert_eq!(cplx.is_connected(), true);
    }
}


use algebra::module::matrix::Matrix;


impl<T: VertexLabel, Coeff: Field + std::fmt::Debug + std::fmt::Display> Space<Coeff> for SimplicialCplx<T> 
{
    type Output = Vec<usize>; 
    fn homology(&self, _: Coeff ) -> Self::Output {
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

        let mut homology = Vec::new();
        homology.reserve( self.dim()+1 );
        let mut kernel = im_and_ker[0].1;
        for &(im, ker) in im_and_ker.iter().skip(1) {
            homology.push(kernel - im);
            kernel = ker;
        }
        homology.push(im_and_ker.last().unwrap().1);
        homology
    }

    fn get_boundary_map(&self) -> Vec<Matrix<Coeff>> {
        let mut boundary_maps: Vec<Matrix<Coeff>> = Vec::new();
        if self.maximal_cells.is_empty() {
            return boundary_maps;
        }

        let mut join = self.maximal_cells[0].clone();
        self.maximal_cells.iter().skip(1).for_each(|maximal_cell| join.join(maximal_cell) );
        boundary_maps.resize(self.dim()+1, Matrix::zero(1,1));
        let mut n_cells = self.get_i_cells(self.dim(), None);
        // iterating over dimensions
        for n in (1..=self.dim()).rev() {
            let n_minus_1_cells = self.get_i_cells(n-1, Some(&n_cells));

            let mut boundary_map = Matrix::zero(n_minus_1_cells.len(), n_cells.len());

            // iterate over n-cells
            for (j, n_cell) in n_cells.iter().enumerate() {
                let boundary = n_cell.get_face();
                let mut sign = Coeff::one();
                // iterate over (n-1)-cells
                n_minus_1_cells.iter()
                    .enumerate()
                    .filter( |(_, n_minus_1_cell)| boundary.contains(n_minus_1_cell))
                    .for_each(|(i, _)| {
                        boundary_map.write(i, j, sign);
                        sign = sign * (Coeff::zero() - Coeff::one());
                } );
            }
            n_cells = n_minus_1_cells;
            boundary_maps[n] = boundary_map;
        }

        // for dimension zero
        boundary_maps[0] = Matrix::zero(1, n_cells.len());

        println!("boundary_maps: {boundary_maps:?}");

        boundary_maps
    }


}#[cfg(test)]
mod simpilicial_cplx_space_tests {
    use crate::complex::*;
    use algebra::commutative::rational::*;
    use algebra::commutative::Zero;

    #[test]
    fn homology_test() {
        let mut circle = SimplicialCplx::new();
        let edge1 = Simplex::new_from( vec!["v1".to_string(), "v2".to_string()] );
        let edge2 = Simplex::new_from( vec!["v2".to_string(), "v3".to_string()] );
        let edge3 = Simplex::new_from( vec!["v3".to_string(), "v1".to_string()] );
        circle.attach(edge1.clone());
        circle.attach(edge2.clone());
        circle.attach(edge3.clone());

        let homology = circle.homology(Rational::zero());
        assert_eq!(homology, vec![1,1]);

        let mut sphere = SimplicialCplx::new();
        let faces = 
            vec![vec!["v0".to_string(), "v1".to_string(), "v2".to_string()],
                 vec!["v0".to_string(), "v2".to_string(), "v3".to_string()],
                 vec!["v0".to_string(), "v3".to_string(), "v4".to_string()],
                 vec!["v0".to_string(), "v4".to_string(), "v1".to_string()],
                 vec!["v5".to_string(), "v1".to_string(), "v2".to_string()],
                 vec!["v5".to_string(), "v2".to_string(), "v3".to_string()],
                 vec!["v5".to_string(), "v3".to_string(), "v4".to_string()],
                 vec!["v5".to_string(), "v4".to_string(), "v1".to_string()]];
        let faces: Vec<_> = faces.into_iter().map(|v| Simplex::new_from(v)).collect();
        faces.iter().for_each(|face| sphere.attach(face.clone()));

        let homology = sphere.homology(Rational::zero());
        assert_eq!(homology, vec![1,0,1]);

        let mut torus = SimplicialCplx::new();
        let faces = 
            vec![vec!["v0".to_string(), "v1".to_string(), "v3".to_string()],
                 vec!["v1".to_string(), "v3".to_string(), "v5".to_string()],
                 vec!["v1".to_string(), "v2".to_string(), "v5".to_string()],
                 vec!["v0".to_string(), "v2".to_string(), "v7".to_string()],
                 vec!["v0".to_string(), "v3".to_string(), "v7".to_string()],
                 vec!["v2".to_string(), "v5".to_string(), "v7".to_string()],
                 vec!["v3".to_string(), "v4".to_string(), "v5".to_string()],
                 vec!["v4".to_string(), "v5".to_string(), "v6".to_string()],
                 vec!["v5".to_string(), "v6".to_string(), "v7".to_string()],
                 vec!["v6".to_string(), "v7".to_string(), "v8".to_string()],
                 vec!["v3".to_string(), "v7".to_string(), "v8".to_string()],
                 vec!["v3".to_string(), "v4".to_string(), "v8".to_string()],
                 vec!["v2".to_string(), "v4".to_string(), "v8".to_string()],
                 vec!["v0".to_string(), "v2".to_string(), "v4".to_string()],
                 vec!["v1".to_string(), "v2".to_string(), "v8".to_string()],
                 vec!["v1".to_string(), "v6".to_string(), "v8".to_string()],
                 vec!["v0".to_string(), "v1".to_string(), "v6".to_string()],
                 vec!["v0".to_string(), "v4".to_string(), "v6".to_string()]];
        let faces: Vec<_> = faces.into_iter().map(|v| Simplex::new_from(v)).collect();
        faces.iter().for_each(|face| torus.attach(face.clone()));

        let homology = torus.homology(Rational::zero());
        assert_eq!(homology, vec![1,2,1]);
    }
}
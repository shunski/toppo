use crate::simplicial::{Simplex, SimplicialCell, SimplicialCplx, VertexLabel};
use crate::Complex;
use util::bitwise::count_ones;

#[allow(unused)]
impl SimplicialCell {
    pub fn new(s: usize, mut raw_vertices: Vec<usize>) -> SimplicialCell {
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
            dim: Some(0),
        };

        raw_vertices.iter().for_each(|&x| cell.cone(x));

        cell.dim = SimplicialCell::dim(&cell.vertices);

        cell
    }

    #[inline]
    pub fn cone(&mut self, vertex: usize) {
        if vertex >= self.size {
            let new_size = vertex+1;
            self.resize(new_size);
        }
        self.vertices[vertex/64] |= 1 << (vertex % 64);
        self.dim = if let Some(d)=self.dim {
            Some(d+1)
        } else {
            Some(0)
        }
    }

    pub fn join(&mut self, cell: &SimplicialCell) {
        self.vertices = self.vertices
            .iter()
            .zip(cell.vertices.iter()).map(|(&x, &y)| x | y)
            .collect();
            
        if cell.vertices.len() > self.vertices.len() { 
            (self.vertices.len()..cell.vertices.len()).
                for_each(|i| self.vertices.push(cell.vertices[i]));
        }

        self.dim = SimplicialCell::dim(&self.vertices);
    }

    pub fn intersection(&self, cell: &SimplicialCell) -> SimplicialCell {
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

    pub fn dim(vertices: &Vec<u64>) -> Option<usize> {
        let dim_plus_1 = 
            vertices.iter()
            .map( |&vertex| {
                count_ones(vertex)
            } )
            .sum::<usize>();
        if dim_plus_1 > 0 {
            Some( dim_plus_1-1 )
        } else {
            None
        }
    }

    #[inline] 
    pub fn resize(&mut self, new_size: usize) {
        while new_size > 64 * self.vertices.len() {
            self.vertices.push(0);
        }
        self.size = new_size;
    }

    pub fn is_connected_with(&self, cell: &SimplicialCell) -> bool {
        self.vertices.iter().zip(cell.vertices.iter()).any(|(&x, &y)| x & y != 0)
    }

    pub fn is_empty(&self) -> bool {
        self.vertices.iter().all(|&x| x==0)
    }

    pub fn contains(&self, cell: &SimplicialCell) -> bool {
        self.vertices.iter().zip(cell.vertices.iter()).all(|(&x, &y)| x | y == x)
    }

    pub fn get_face(&self) -> Vec<SimplicialCell> {
        let n = if let Some(d) = self.dim {
            d-1
        } else {
            panic!();
        };

        let mut n_sub_cells = Vec::new();
        n_sub_cells.reserve(n+2);
        let mut idx = 0;
        for j in 0..n+2 {
            let mut sub_cell = self.clone();
            sub_cell.dim = if let Some(d) = sub_cell.dim {
                Some(d-1)
            } else {
                panic!();
            };
            idx = (idx..).find(|&new_idx| (sub_cell.vertices[new_idx/64]>>(new_idx%64)) % 2 == 1).unwrap();
            sub_cell.vertices[idx/64] &= !(1 << (idx%64));
            n_sub_cells.push(sub_cell);
            idx+=1;
        };
        n_sub_cells
    }
}

#[allow(unused)]
impl<T: VertexLabel> Simplex<T> 
{
    pub fn from(labels: Vec<T>) -> Simplex<T> {
        Simplex{
            dim: labels.len() - 1,
            vertices_labels: Vec::from_iter(labels.into_iter()),
        }
    }
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

    fn get_simplex(&self, cell: &SimplicialCell) -> Simplex<T> {
        let vertices_labels = (0..cell.size)
            .filter(|&i| cell.contains( &SimplicialCell::new( i+1, vec![i] )))
            .map(|i| self.vertices_labeling[i].clone())
            .collect::<Vec<_>>();
        Simplex { dim: cell.dim.unwrap(), vertices_labels: vertices_labels }
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
            .map(|cell| cell.dim.unwrap_or(0))
            .max()
            .unwrap()
    }

    pub fn is_connected(&self) -> bool {
        let join = self.maximal_cells.first();
        let mut join = match join {
            Some(x) => x.clone(),
            None => return true,
            // if the complex is empty, then it is vacuously connected.
        };
        let mut prev_dim = join.dim;
        loop {
            // iterate each cell and join it if connected.
            let cells_connected_to_join: Vec<_> = self.maximal_cells.iter()
                .filter(|&cell| cell.is_connected_with(&join))
                .collect();
            println!("cells_connected_to_join.len(): {}", cells_connected_to_join.len());
            cells_connected_to_join.iter().for_each(|cell| { join.join(cell); });

            // if the join eventually gets the join of all points, then the complex is connected
            if join.dim == Some(self.size-1) {
                return true;
            }

            // if the join comes here and if the join stoped growing, then the complex is not connected.
            if join.dim == prev_dim {
                return false;
            }

            prev_dim = join.dim;
        }
    }

    #[inline]
    pub fn n_vertices(&self) -> usize {
        self.vertices_labeling.len()
    }

    pub fn n_cells_of_dim(&self, dim: usize) -> usize {
        use util::TupleIterable;

        // closure that computes the binomial
        let binomial = |x: usize, y: usize| { 
            if x < y {
                0
            } else {
                (x-y+1..=x).product::<usize>() / (1..=y).product::<usize>()
            }
        };

        // Perform the principle of inclusion and exclusion
        // Itertate over the number of subsets
        let mut sign = true;
        let mut out: i64 = 0;
        let mut intersections = self.maximal_cells.clone();
        for i in 1..=self.maximal_cells.len() {
            // take all possible intersections of i maximal cells
            let mut intersections = Vec::new();
            let largest_cell = SimplicialCell{
                vertices: vec![!0_u64; self.size/64+1],
                size: self.size,
                dim: Some(self.size),
            };
            for cells in self.maximal_cells.tuple_iter(i) {
                let intersection = cells.iter()
                    .fold(largest_cell.clone(), |accum, &cell| accum.intersection(cell));
                if !intersection.is_empty() {
                    intersections.push(intersection);
                }
            }
            if intersections.is_empty() {
                break;
            }


            let sum = intersections.iter()
                .map( |cell| binomial(cell.dim.unwrap()+1, dim+1))
                .sum::<usize>();
            if sign {
                out += sum as i64;
            } else {
                out -= sum as i64;
            }

            // sign must alternate
            sign = !sign;
        }
        out as usize
    }

    fn get_i_cells(&self, i: usize, i_plus_1_cells: Option<&Vec<SimplicialCell>>) -> Vec<SimplicialCell> {
        let mut i_cells:Vec<_> = self.maximal_cells
            .iter()
            .filter(|&maximal_cell|maximal_cell.dim.unwrap() == i)
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


use alg::lin_alg::Matrix;
impl<T: VertexLabel> Complex for SimplicialCplx<T> 
{
    type Cell = Simplex<T>;
    fn boundary_map(&self) -> Vec<( Matrix<i128>, Vec<Self::Cell> )> {
        let mut boundary_maps: Vec<Matrix<i128>> = vec![Matrix::zero(1,1); self.dim()+1];
        let mut basis: Vec<Vec<Simplex<T>>> = vec![Vec::new(); self.dim()+1];

        let mut n_cells = self.get_i_cells(self.dim(), None);
        // iterating over dimensions
        for n in (1..=self.dim()).rev() {
            basis[n] = n_cells.iter().map(|cell| self.get_simplex(cell)).collect();

            let n_minus_1_cells = self.get_i_cells(n-1, Some(&n_cells));

            let mut boundary_map = Matrix::zero(n_minus_1_cells.len(), n_cells.len());

            // iterate over n-cells
            for (j, n_cell) in n_cells.iter().enumerate() {
                let boundary = n_cell.get_face();
                let mut sign = 1;
                // iterate over (n-1)-cells
                n_minus_1_cells.iter()
                    .enumerate()
                    .filter( |(_, n_minus_1_cell)| boundary.contains(n_minus_1_cell))
                    .for_each(|(i, _)| {
                        boundary_map[( i, j )] = sign;
                        sign = -sign;
                } );
            }
            n_cells = n_minus_1_cells;
            boundary_maps[n] = boundary_map;
        }

        // for dimension zero
        boundary_maps[0] = Matrix::zero(1, n_cells.len());
        basis[0] = n_cells.iter().map(|cell| self.get_simplex(cell)).collect();

        let boundary_maps = boundary_maps.into_iter().zip( basis.into_iter() ).collect::<Vec<_>>();

        boundary_maps
    }
}
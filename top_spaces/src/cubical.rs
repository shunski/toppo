use crate::Space;
use crate::graph::*;
use algebra::commutative::Field;
use algebra::module::matrix::Matrix;

#[allow(unused)]
#[derive(Clone)]
struct Branch {
    parent_idx: Option<usize>,
    child_indeces: Vec<usize>,
    val: CubeFactor,
    curr_dim: usize,
}

#[allow(unused)]
impl Branch {
    fn new_parent(val: CubeFactor) -> Branch {
        let curr_dim = match val {
            CubeFactor::Edge(_, _) => 1,
            CubeFactor::Vertex(_) => 0,
        };
        Branch { 
            parent_idx: None, 
            child_indeces: Vec::new(), 
            val: val,
            curr_dim: curr_dim,
        }
    }

    fn new(val: CubeFactor, parent_idx: usize, cplx: &CubicalCplx) -> Branch {
        let mut curr_dim = cplx.data[parent_idx].dim();
        curr_dim += match val {
            CubeFactor::Edge(_, _) => 1,
            CubeFactor::Vertex(_) => 0,
        };
        Branch { 
            parent_idx: Some(parent_idx), 
            child_indeces: Vec::new(), 
            val: val,
            curr_dim: curr_dim,
        }
    }

    pub fn dim(&self) -> usize {
        self.curr_dim
    }

    fn is_parent(&self) -> bool {
        return self.parent_idx == None
    }

    fn get_parent_idx(&self) -> Option<usize> {
        self.parent_idx
    }
}

#[allow(unused)]
pub struct CubicalCplx {
    factor_pool: SimpleGraph,
    data: Vec<Branch>,
    first_branch_indeces: Vec<usize>,
    final_branch_indeces: Vec<usize>,
    dim: usize,
    ranks_of_chain_cplx: Vec<usize>,
}

#[allow(unused)]
impl CubicalCplx {
    pub fn udc(graph: &Graph, n_marked_points: usize) -> CubicalCplx {
        let factor_pool = SimpleGraph::subdivided_graph_from(graph, n_marked_points);

        // if graph is empty or consisting of only one vertex, then we panic. TODO: handle these boundary cases.
        if factor_pool.n_vertices() < 2 { panic!() };

        let boundary_vertex = CubeFactor::Vertex( factor_pool.n_vertices() - (n_marked_points - 1) );
        let data: Vec<_> = factor_pool
            .factor_iter()
            .filter(|&factor| !factor.is_in_maximal_tree() || factor.is_disjoint_from(boundary_vertex))
            .filter(|&factor| !factor.is_in_maximal_tree() || factor < boundary_vertex)
            .map(|factor| Branch::new_parent(factor))
            .collect();

        let first_branch_indeces: Vec<_> = (0..data.len()).collect();

        let mut cplx = CubicalCplx {
            factor_pool: factor_pool,
            data: data,
            first_branch_indeces: first_branch_indeces.clone(),
            final_branch_indeces: Vec::new(),
            dim: 0,
            ranks_of_chain_cplx: vec![0; n_marked_points+1],
        };
        
        first_branch_indeces.iter().for_each(|&idx| cplx.expand(idx, n_marked_points-1));

        cplx
    }

    fn expand(&mut self, branch_idx: usize, n_rem_factors: usize) {
        // base case
        if n_rem_factors == 0 { 
            self.final_branch_indeces.push(branch_idx);
            self.ranks_of_chain_cplx[self.data[branch_idx].dim()] += 1;
            if self.data[branch_idx].curr_dim > self.dim { self.dim = self.data[branch_idx].curr_dim; }
            return; 
        }

        // general case
        let boundary_vertex = CubeFactor::Vertex(self.factor_pool.n_vertices() - (n_rem_factors-1));
        // find all the branches that we need to add
        let new_branches: Vec<_> = self.factor_pool
            .factor_iter()
            .filter(|&factor| !factor.is_in_maximal_tree() || factor.is_disjoint_from(boundary_vertex))
            .filter(|&factor| !factor.is_in_maximal_tree() || factor < boundary_vertex)
            .filter(|&factor| self.factor_iter(branch_idx).all(|x| x.is_disjoint_from(factor)))
            .filter(|&factor| factor > self.data[branch_idx].val)
            .map(|factor| Branch::new(factor, branch_idx, self ))
            .collect();

        // for each of the new branches, ...
        new_branches.into_iter().for_each( |new_branch| {
            let new_branch_idx = self.data.len();
            self.data[branch_idx].child_indeces.push(new_branch_idx);
            self.data.push(new_branch);
            self.expand(new_branch_idx, n_rem_factors-1);
        });
    }


    fn get_cell_from_idx(&self, idx: usize) -> Option<Vec<CubeFactor>> {
        if self.final_branch_indeces.contains(&idx) {
            let mut out: Vec<_> = self.factor_iter(idx).collect();
            out.reverse();
            Some(out)
        } else {
            None
        }
    }

    fn get_boundaries_of(&self, branch_idx: usize) -> Vec<(bool, usize)> {
        let mut boundaries = Vec::new();
        let cell = self.get_cell_from_idx(branch_idx).unwrap();

        let mut sign = true;
        for (i, factor) in cell.iter().enumerate().filter(|&(_, factor)| factor.is_edge()) {
            let (mut first_face, mut second_face) = (cell.clone(), cell.clone());
            first_face[i] = CubeFactor::Vertex(factor.edge().0);
            second_face[i] = CubeFactor::Vertex(factor.edge().1);
            boundaries.push((sign, first_face)); boundaries.push((!sign, second_face));

            sign = !sign;
        }

        boundaries.iter_mut().for_each(|(bool, cell)| cell.sort());
        let boundaries = boundaries.into_iter().map(|(sign, cell)| (sign, self.get_index_of(&cell).unwrap())).collect();

        boundaries
    }

    fn get_index_of(&self, cell: &[CubeFactor]) -> Option<usize> {
        let mut branch_idx = match self.first_branch_indeces.iter().find(|&&idx| self.data[idx].val == cell[0]) {
            Some(&idx) => idx,
            None => { return None; },
        };
        for i in 1..cell.len() {
            branch_idx = match self.data[branch_idx].child_indeces.iter().find(|&&idx| self.data[idx].val == cell[i] ) {
                Some(&idx) => idx,
                None => { return None; },
            };
        };
        Some(branch_idx)
    }

    fn contains(&self, cell: &[CubeFactor]) -> bool {
        match self.get_index_of(cell) {
            Some(_) => true,
            None => false,
        }
    }

    pub fn n_cells(&self) -> usize {
        self.final_branch_indeces.len()
    }
}


struct CubicalFactorIter<'a> {
    cplx: &'a CubicalCplx,
    curr_idx: usize,
    is_done: bool,
}

impl<'a> Iterator for CubicalFactorIter<'a> {
    type Item = CubeFactor;

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_done { return None; }
        let out = self.cplx.data[self.curr_idx].val;
        let parent = self.cplx.data[self.curr_idx].get_parent_idx();

        match parent {
            Some(idx) => self.curr_idx = idx,
            None => self.is_done = true,
        }

        Some(out)
    }
}

impl<'a> CubicalFactorIter<'a> {
    fn new(cplx: &'a CubicalCplx, idx: usize) -> CubicalFactorIter<'a> {
        CubicalFactorIter {
            cplx: cplx,
            curr_idx: idx,
            is_done: false,
        }
    }
}


impl<'a> CubicalCplx {
    fn factor_iter( &'a self, idx: usize ) -> CubicalFactorIter<'a> {
        CubicalFactorIter::new(self, idx)
    }
}

use std::fmt;
impl fmt::Display for CubicalCplx {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for i in 0..=self.dim {
            write!(f, "{}-cells:", i)?;
            for &cell_idx in self.final_branch_indeces.iter().filter(|&&idx| self.data[idx].dim() == i) {
                // printing each cell
                write!(f, "\n\t{{")?;
                for &factor in self.get_cell_from_idx(cell_idx).unwrap().iter() {
                    // printping each factor
                    write!(f, "{:?}", factor)?;
                    // match factor {
                    //     CubeFactor::Edge(v, w) => write!(f, "[{}, {}]", self.factor_pool[v], self.factor_pool[w])?,
                    //     CubeFactor::Vertex(v) => write!(f, "{}", self.factor_pool[v])?,
                    // };
                    write!(f, ", ")?;
                }
                write!(f, "}}")?;
                if cell_idx != *self.final_branch_indeces.last().unwrap() {
                    write!(f, ", ")?;
                };
            }
            write!(f, "\n")?;
        };
        write!(f, "")
    }
}


#[cfg(test)]
mod cubical_cells {
    use crate::graph;
    use crate::graph::*;
    use crate::cubical::CubicalCplx;
    #[test]
    fn dc_init_test() {
        // test 1
        let graph = graph!{
            ["v0", "v1"],
            ["v1", "v1"]
        };
        let cplx = CubicalCplx::udc(&graph, 2);
        println!("{}", cplx);
        assert_eq!(cplx.ranks_of_chain_cplx, vec![6, 8, 1]);
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(1)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(2)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(3)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(1), CubeFactor::Vertex(2)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(1), CubeFactor::Vertex(3)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(2), CubeFactor::Vertex(3)]));

        // test 2
        let cplx = CubicalCplx::udc(&graph, 3);
        println!("{}", cplx);
        assert_eq!(cplx.ranks_of_chain_cplx, vec![20, 36, 16, 2]);
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(1), CubeFactor::Vertex(2)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(1), CubeFactor::Vertex(3)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(1), CubeFactor::Vertex(4)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(1), CubeFactor::Vertex(5)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(2), CubeFactor::Vertex(3)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(2), CubeFactor::Vertex(4)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(2), CubeFactor::Vertex(5)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(3), CubeFactor::Vertex(4)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(3), CubeFactor::Vertex(5)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(4), CubeFactor::Vertex(5)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(1), CubeFactor::Vertex(2), CubeFactor::Vertex(3)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(1), CubeFactor::Vertex(2), CubeFactor::Vertex(4)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(1), CubeFactor::Vertex(2), CubeFactor::Vertex(5)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(1), CubeFactor::Vertex(3), CubeFactor::Vertex(4)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(1), CubeFactor::Vertex(3), CubeFactor::Vertex(5)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(1), CubeFactor::Vertex(4), CubeFactor::Vertex(5)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(2), CubeFactor::Vertex(3), CubeFactor::Vertex(4)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(2), CubeFactor::Vertex(3), CubeFactor::Vertex(5)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(2), CubeFactor::Vertex(4), CubeFactor::Vertex(5)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(3), CubeFactor::Vertex(4), CubeFactor::Vertex(5)]));


        // test 3
        let graph = graph!{
            ["v0", "v1"],
            ["v0", "v2"],
            ["v0", "v3"]
        };
        let cplx = CubicalCplx::udc(&graph, 3);
        println!("{}", cplx);
        assert_eq!(cplx.ranks_of_chain_cplx, vec![(7*6*5)/(3*2*1), (5*4)/(2*1) * 6, 4*3 + 2*3 + 3*3, 4]);
        // assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(1), CubeFactor::Vertex(6)]));
    }
}

impl<Coeff: Field> Space<Coeff> for CubicalCplx {
    fn get_boundary_map (&self) -> Vec<Matrix<Coeff>> {
        let mut numbering = vec![0; self.data.len()];
        let mut counter = vec![0; self.dim+1];
        self.final_branch_indeces.iter().for_each(|&idx| {
            numbering[idx] = counter[self.data[idx].dim()];
            counter[self.data[idx].dim()] += 1;
        });

        // debug_assert!();

        // instantiating the boundary maps
        let mut boundary_maps = vec![Matrix::zero(1, self.ranks_of_chain_cplx[0])];
        (1..=self.dim).for_each(|i| boundary_maps.push(Matrix::zero(self.ranks_of_chain_cplx[i-1], self.ranks_of_chain_cplx[i])));

        // filling up the boundary maps
        for &idx in self.final_branch_indeces.iter().filter(|&&idx| self.data[idx].dim() > 0 ) {
            let boundary_indeces = self.get_boundaries_of(idx);
            let dim = self.data[idx].dim();
            boundary_indeces.iter().for_each(|&(sign, b_idx)| {
                let val = boundary_maps[dim].get(numbering[b_idx], numbering[idx]);
                boundary_maps[dim].write(numbering[b_idx],numbering[idx], if sign {val+Coeff::one()} else {val+Coeff::negative_one()})
            });
        }

        boundary_maps
    }
}


#[cfg(test)]
mod homology_test {
    use algebra::commutative::rational::Rational;
    use crate::{graph, Space};
    use crate::cubical::CubicalCplx;
    #[test]
    fn homology_test() {
        // test 1
        let graph = graph!{
            ["v0", "v1"],
            ["v0", "v2"],
            ["v0", "v3"]
        };
        let cplx = CubicalCplx::udc(&graph, 2);
        assert_eq!( H!(cplx; Rational), vec![1, 1] );

        // test 2
        let cplx = CubicalCplx::udc(&graph, 3);
        assert_eq!( H!(cplx; Rational), vec![1,3,0,0] );

        // test 3
        let graph = graph!{
            ["v0", "v1"],
            ["v1", "v1"]
        };
        let cplx = CubicalCplx::udc(&graph, 5);
        assert_eq!( H!(cplx; Rational), vec![1, 5, 0, 0, 0, 0] );

        // test 4
        let graph = graph!{
            ["v1", "v2"],
            ["v1", "v3"],
            ["v1", "v4"],
            ["v2", "v3"],
            ["v2", "v4"],
            ["v3", "v4"]
        };
        let cplx = CubicalCplx::udc(&graph, 2);
        assert_eq!( H!(cplx; Rational), vec![1, 4, 0] );
        let cplx = CubicalCplx::udc(&graph, 3);
        assert_eq!( H!(cplx; Rational), vec![1, 4, 3, 0] );        

        // test 5
        let graph = graph!{
            ["v1", "v2"],
            ["v1", "v3"],
            ["v1", "v4"],
            ["v1", "v5"],
            ["v2", "v3"],
            ["v2", "v4"],
            ["v2", "v5"],
            ["v3", "v4"],
            ["v3", "v5"],
            ["v4", "v5"]
        };
        let cplx = CubicalCplx::udc(&graph, 2);
        assert_eq!( H!(cplx; Rational), vec![1, 6, 0] );
        let cplx = CubicalCplx::udc(&graph, 3);
        assert_eq!( H!(cplx; Rational), vec![1, 6, 30, 0] );
    }
}
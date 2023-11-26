use crate::Complex;
use crate::graph::*;
use alg::formal_sum_impl;
use alg::lin_alg::Vector;
use alg::lin_alg::Matrix;
use alg::non_commutative::permutation::ConstPermutation;
use util::LabeledSet;

#[allow(unused)]
#[derive(Clone, Debug)]
struct Branch {
    parent_idx: Option<usize>,
    child_indeces: Vec<usize>,
    val: CubeFactor,
}

#[allow(unused)]
impl Branch {
    fn new_parent(val: CubeFactor) -> Branch {
        Branch { 
            parent_idx: None, 
            child_indeces: Vec::new(), 
            val: val,
        }
    }

    fn new(val: CubeFactor, parent_idx: usize ) -> Branch {
        Branch { 
            parent_idx: Some(parent_idx), 
            child_indeces: Vec::new(), 
            val: val,
        }
    }

    fn is_parent(&self) -> bool {
        return self.parent_idx == None
    }

    fn get_parent_idx(&self) -> Option<usize> {
        self.parent_idx
    }
}


impl std::ops::Deref for Branch {
    type Target = CubeFactor;

    fn deref(&self) -> &Self::Target {
        &self.val
    }
}

impl std::ops::DerefMut for Branch {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.val
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
    pub fn udc(graph: &Graph<i64>, n_marked_points: usize) -> CubicalCplx {
        let factor_pool = SimpleGraph::subdivided_graph_from(graph, n_marked_points);

        // if graph is empty or consisting of only one vertex, then we panic. TODO: handle these boundary cases.
        if factor_pool.n_vertices() < 2 { panic!() };

        let boundary_vertex = CubeFactor::Vertex( factor_pool.n_vertices() - (n_marked_points - 1) );
        let data: Vec<_> = factor_pool
            .iter()
            .filter(|&factor| (factor.is_edge() && !factor.is_order_respecting()) || factor.is_disjoint_from(boundary_vertex))
            .filter(|&factor| (factor.is_edge() && !factor.is_order_respecting()) || factor < boundary_vertex)
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
            let dim_of_cell = self.dim_of(branch_idx);
            self.ranks_of_chain_cplx[dim_of_cell] += 1;
            if dim_of_cell > self.dim { self.dim = dim_of_cell; }
            return; 
        }

        // general case
        let boundary_vertex = CubeFactor::Vertex(self.factor_pool.n_vertices() - (n_rem_factors-1));
        // find all the branches that we need to add
        let new_branches: Vec<_> = self.factor_pool
            .iter()
            .filter(|&factor| (factor.is_edge() && !factor.is_order_respecting()) || factor.is_disjoint_from(boundary_vertex))
            .filter(|&factor| (factor.is_edge() && !factor.is_order_respecting()) || factor < boundary_vertex)
            .filter(|&factor| self.factor_iter(branch_idx).all(|x| x.is_disjoint_from(factor)))
            .filter(|&factor| factor > self.data[branch_idx].val)
            .map(|factor| Branch::new(factor, branch_idx))
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
            first_face[i] = CubeFactor::Vertex(factor.edge()[0]);
            second_face[i] = CubeFactor::Vertex(factor.edge()[1]);
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

    pub fn dim_of(&self, branch_idx: usize) -> usize {
        let dim_of_branch = match self.data[branch_idx].val {
            CubeFactor::Edge(_) => 1,
            CubeFactor::Vertex(_) => 0,
        };

        // return the recursive call
        dim_of_branch + if let Some(p_idx) = self.data[branch_idx].parent_idx {
            self.dim_of(p_idx)
        } else {
            0
        }
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


impl<'a> CubicalCplx {
    fn factor_iter( &'a self, idx: usize ) -> CubicalFactorIter<'a> {
        CubicalFactorIter {
            cplx: self,
            curr_idx: idx,
            is_done: false,
        }
    }
}

use std::collections::VecDeque;
use std::fmt;
impl fmt::Display for CubicalCplx {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for i in 0..=self.dim {
            write!(f, "{}-cells:", i)?;
            for &cell_idx in self.final_branch_indeces.iter().filter(|&&idx| self.dim_of(idx) == i) {
                // printing each cell
                write!(f, "\n\t{{")?;
                for &factor in self.get_cell_from_idx(cell_idx).unwrap().iter() {
                    // printping each factor
                    write!(f, "{}", factor)?;
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
        assert_eq!(cplx.ranks_of_chain_cplx, vec![10, 15, 4]);
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(1)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(2)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(3)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(0), CubeFactor::Vertex(4)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(1), CubeFactor::Vertex(2)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(1), CubeFactor::Vertex(3)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(1), CubeFactor::Vertex(4)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(2), CubeFactor::Vertex(3)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(2), CubeFactor::Vertex(4)]));
        assert!(cplx.contains(&vec![CubeFactor::Vertex(3), CubeFactor::Vertex(4)]));

        // test 2
        let cplx = CubicalCplx::udc(&graph, 3);
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

impl Complex for CubicalCplx {
    type Cell = Cube;
    fn boundary_map (&self) -> Vec< (Matrix<i64>, Vec<Self::Cell>) > {
        let basis: Vec<Vec<Self::Cell>> = vec![Vec::new(); self.dim+1];

        let mut numbering = vec![0; self.data.len()];
        let mut counter = vec![0; self.dim+1];
        self.final_branch_indeces.iter().for_each(|&idx| {
            numbering[idx] = counter[self.dim_of(idx)];
            counter[self.dim_of(idx)] += 1;
        });

        // instantiating the boundary maps
        let mut boundary_maps = vec![Matrix::zero(1, self.ranks_of_chain_cplx[0])];
        (1..=self.dim).for_each(|i| boundary_maps.push(Matrix::zero(self.ranks_of_chain_cplx[i-1], self.ranks_of_chain_cplx[i])));

        // filling up the boundary maps
        for &idx in self.final_branch_indeces.iter().filter(|&&idx| self.dim_of(idx) > 0 ) {
            let boundary_indeces = self.get_boundaries_of(idx);
            let dim = self.dim_of(idx);
            boundary_indeces.iter().for_each(|&(sign, b_idx)| {
                let val = boundary_maps[dim][( numbering[b_idx], numbering[idx] )];
                boundary_maps[dim][( numbering[b_idx], numbering[idx] )] = if sign {val+1} else {val-1}
            });
        }

        let boundary_maps = boundary_maps.into_iter().zip(basis.into_iter()).collect::<Vec<_>>();
        boundary_maps
    }
}


// #[cfg(test)]
// mod homology_test {
//     use crate::graph;
//     use crate::cubical::CubicalCplx;
//     use crate::IntegralHomology;
    
    // #[test]
    // fn homology_test() {
    //     // test 1
    //     let graph = graph!{
    //         ["v0", "v1"],
    //         ["v0", "v2"],
    //         ["v0", "v3"]
    //     };
    //     let cplx = CubicalCplx::udc(&graph, 2);
    //     assert_eq!( H!(cplx), "dim 0: Z,\ndim 1: Z,\n" );

    //     // test 2
    //     let cplx = CubicalCplx::udc(&graph, 3);
    //     assert_eq!( H!(cplx), "dim 0: Z,\ndim 1: Z^3,\n" );

        // // test 3
        // let graph = graph!{
        //     ["v0", "v1"],
        //     ["v1", "v1"]
        // };
        // let cplx = CubicalCplx::udc(&graph, 3);
        // assert_eq!( H!(cplx), "dim 0: Z,\ndim 1: Z^3,\n" );
        // let cplx = CubicalCplx::udc(&graph, 5);
        // assert_eq!( H!(cplx), "dim 0: Z,\ndim 1: Z^5,\n" );

        // // test 4
        // let graph = graph!{
        //     ["v1", "v2"],
        //     ["v1", "v3"],
        //     ["v1", "v4"],
        //     ["v2", "v3"],
        //     ["v2", "v4"],
        //     ["v3", "v4"]
        // };
        // let cplx = CubicalCplx::udc(&graph, 2);
        // assert_eq!( H!(cplx), "dim 0: Z,\ndim 1: Z^4,\n" );
        // let cplx = CubicalCplx::udc(&graph, 3);
        // assert_eq!( H!(cplx), "dim 0: Z,\ndim 1: Z^4,\ndim 2: Z^3,\n" );

        // // test 5
        // let graph = graph!{
        //     ["v1", "v2"],
        //     ["v1", "v3"],
        //     ["v1", "v4"],
        //     ["v1", "v5"],
        //     ["v2", "v3"],
        //     ["v2", "v4"],
        //     ["v2", "v5"],
        //     ["v3", "v4"],
        //     ["v3", "v5"],
        //     ["v4", "v5"]
        // };
        // let cplx = CubicalCplx::udc(&graph, 2);
        // assert_eq!( H!(cplx), "dim 0: Z,\ndim 1: Z^6 + Z_2,\n" );
        // let cplx = CubicalCplx::udc(&graph, 3);
        // assert_eq!( H!(cplx), "dim 0: Z,\ndim 1: Z^6 + Z_2,\ndim 2: Z^30,\n" );
//     }
// }


#[derive(Debug, PartialEq)]
enum MorseProp {
    Critical,
    Collapsible(usize),
    Redundant(usize),
}


#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Cube (Vec<CubeFactor>);

formal_sum_impl!{ Cube }

impl std::ops::Deref for Cube {
    type Target = Vec<CubeFactor>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for Cube {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<Vec<CubeFactor>> for Cube {
    fn from(value: Vec<CubeFactor>) -> Self {
        Self(value)
    }
}

impl Cube {
    fn from_index( idx: usize, cplx: &UdnMorseCplx ) -> Cube {
        let mut factors: Vec<_> = cplx.factor_iter(idx).collect();
        factors.reverse();
        Cube( factors )
    }

    pub fn flow(mut self, graph: &RawSimpleGraph) -> Vector<i64, Cube> {
        match self.property(graph) {
            MorseProp::Collapsible(_) => Vector::zero(),
            MorseProp::Critical => 1*self,
            MorseProp::Redundant(i) => {
                // do the special reduction if applicable: if an edge starting from unblocked vertex does not disrespect any vertex, 
                // then we can flow to the terminal vertex of the edge
                let flow = self[i].flow(graph);
                let mut vertices_in_cube = Vec::new(); 
                self.iter()
                    .filter(|&&f| f != flow.initial() ) // remove the unblocked vertex itself
                    .for_each(|&f| 
                        if f.is_edge() {
                            vertices_in_cube.push( f.initial() );
                            vertices_in_cube.push( f.terminal() );
                        } else { // if f is a vertex
                            vertices_in_cube.push( f );
                        }
                    );
                if vertices_in_cube.iter().all( |&v| v < flow.terminal() || flow.initial() < v ) {
                    // apply the reduction
                    // let sign = if self.iter().take_while(|&&f| f.is_disjoint_from(flow) ).filter(|&&f| f.is_edge()).count() % 2 == 0 {
                    //     1
                    // } else {
                    //     -1
                    // };
                    self[i] = flow.terminal();
                    self.sort();
                    return self.flow(graph);
                }

                // else do the computation following the definition of relative flow
                let mut cpy = self.clone();
                cpy[i] = self[i].flow(&graph);
                cpy.sort();
                let mut boundary = cpy.get_boundary();
                // self[i] = self[i].initial();
                boundary = if boundary[&self] == 1 {
                    - boundary + self
                } else if boundary[&self] == -1 {
                    boundary + self
                } else {
                    panic!("Implementation of boundary map is wrong. ")
                };
                boundary.into_iter().map(|(i, cube)| cube.flow(graph) * i).sum::<Vector<i64,_>>()
            },
        }
    }

    pub fn get_boundary(&self) -> Vector<i64, Cube> {
        let mut boundary = Vector::zero();
        let mut sign = -1;
        for (i, &edge) in self.0.iter().enumerate().filter( |&(_, &f)| f.is_edge() ) {
            let mut b1 = self.clone(); b1[i] = edge.terminal(); b1.sort();
            let mut b2 = self.clone(); b2[i] = edge.initial();  b2.sort();

            boundary = boundary + ( b2 - b1 ) * sign;
            sign *= -1;
        }

        boundary
    }

    fn property(&self, graph: &RawSimpleGraph) -> MorseProp {
        let f = self.0.iter()
            .enumerate()
            .filter(|&(_, &f)| 
                ( f.is_vertex() && !f.is_blocked(&self, graph) ) || 
                ( f.is_edge() && f.is_in_maximal_tree(graph) && ( f.is_order_respecting() || self.0.iter().all(|&other| !other.is_vertex() || !other.is_disrespected_by(f)) ))
            )
            .next();
        match f {
            Some((i, factor)) => if factor.is_vertex() { MorseProp::Redundant(i) } else { MorseProp::Collapsible(i) },
            None => MorseProp::Critical,
        }
    }

    // This function panics if the cube is not an edge, that is, if the cube is not one dimensional
    pub fn get_ends(self) -> [Cube; 2] {
        let mut out = [Self(Vec::new()), Self(Vec::new())];

        let mut edge_count = 0;

        for f in self.0 {
            match f {
                CubeFactor::Vertex(_) => {
                    out[0].push(f); out[1].push(f); 
                },
                CubeFactor::Edge([v,w]) => {
                    out[0].push(CubeFactor::Vertex(v)); 
                    out[1].push(CubeFactor::Vertex(w));
                    edge_count += 1;
                }
            }
        }

        if edge_count != 1 {
            panic!("input to this function must be an edge");
        }

        out
    }
}

impl std::fmt::Display for Cube {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut iter = self.iter();
        write!(f, "{{{}", iter.next().unwrap())?;
        while let Some(factor) = iter.next() {
            write!(f, ", {}", factor)?;
        }
        write!(f, "}}")
    }
}


#[cfg(test)]
mod cube_test {
    use crate::cubical::*;
    use crate::graph::graphs_for_tests;

    #[test]
    fn property_test() {
        let graph = graphs_for_tests::k_5_5_with_4_marked_points();
        assert_eq!( MorseProp::Critical, Cube( vec![ CubeFactor::Vertex(0), CubeFactor::Vertex(1), CubeFactor::Vertex(2), CubeFactor::Vertex(3)] ).property(&graph) );
        assert_eq!( MorseProp::Critical, Cube( vec![ CubeFactor::Edge([6, 9]), CubeFactor::Vertex(7), CubeFactor::Vertex(19), CubeFactor::Edge([18, 23])] ).property(&graph) );
        assert_eq!( MorseProp::Critical, Cube( vec![ CubeFactor::Edge([6, 9]), CubeFactor::Vertex(7), CubeFactor::Vertex(10), CubeFactor::Vertex(11)] ).property(&graph) );
        assert_eq!( MorseProp::Critical, Cube( vec![ CubeFactor::Edge([0, 24]), CubeFactor::Vertex(1), CubeFactor::Edge([6, 9]), CubeFactor::Vertex(7)] ).property(&graph) );

        assert_eq!( MorseProp::Collapsible(1), Cube( vec![ CubeFactor::Vertex(0), CubeFactor::Edge([6, 9]), CubeFactor::Vertex(10), CubeFactor::Vertex(11)] ).property(&graph) );
        assert_eq!( MorseProp::Collapsible(2), Cube( vec![ CubeFactor::Vertex(0), CubeFactor::Edge([3, 22]), CubeFactor::Edge([6, 9]), CubeFactor::Vertex(11)] ).property(&graph) );

        assert_eq!( MorseProp::Redundant(0), Cube( vec![ CubeFactor::Vertex(1), CubeFactor::Edge([6, 9]), CubeFactor::Vertex(10), CubeFactor::Vertex(11)] ).property(&graph) );
        assert_eq!( MorseProp::Redundant(0), Cube( vec![ CubeFactor::Vertex(1), CubeFactor::Edge([3, 22]), CubeFactor::Edge([6, 9]), CubeFactor::Vertex(11)] ).property(&graph) );
        assert_eq!( MorseProp::Redundant(1), Cube( vec![ CubeFactor::Edge([3, 22]), CubeFactor::Vertex(5), CubeFactor::Edge([6, 9]), CubeFactor::Vertex(11)] ).property(&graph) );
        assert_eq!( MorseProp::Redundant(2), Cube( vec![ CubeFactor::Edge([6, 9]), CubeFactor::Vertex(7), CubeFactor::Vertex(12), CubeFactor::Vertex(13)] ).property(&graph) );
    }

    #[test]
    fn cubical_boundary_test(){
        // test 1
        let boundary = Cube( vec![CubeFactor::Vertex(6), CubeFactor::Vertex(7), CubeFactor::Vertex(9), CubeFactor::Vertex(24)] ).get_boundary();
        let ans = Vector::zero();
        assert_eq!(boundary, ans);

        // test 2
        let boundary = Cube( vec![CubeFactor::Edge([6, 9]), CubeFactor::Vertex(7), CubeFactor::Vertex(10), CubeFactor::Vertex(11)] ).get_boundary();
        let boundary: Vector<i64,_> = boundary.into_iter().map( |(i, cube)| i*cube ).sum();
        let ans = Cube( vec![CubeFactor::Vertex(6), CubeFactor::Vertex(7), CubeFactor::Vertex(10), CubeFactor::Vertex(11)] )
            - Cube( vec![CubeFactor::Vertex(7), CubeFactor::Vertex(9), CubeFactor::Vertex(10), CubeFactor::Vertex(11)] );
        assert_eq!(boundary, ans);

        // test 3 (a non order-respecting edge)
        let boundary = Cube( vec![ CubeFactor::Edge([0, 8]), CubeFactor::Vertex(1), CubeFactor::Edge([6, 9]), CubeFactor::Vertex(7) ] ).get_boundary();
        let boundary: Vector<i64,_> = boundary.into_iter().map( |(i, cube)| i*cube ).sum();
        let ans = Cube( vec![ CubeFactor::Vertex(0), CubeFactor::Vertex(1), CubeFactor::Edge([6,9]), CubeFactor::Vertex(7) ] )
            - Cube( vec![ CubeFactor::Vertex(1), CubeFactor::Edge([6,9]), CubeFactor::Vertex(7), CubeFactor::Vertex(8) ] ) 
            - Cube( vec![ CubeFactor::Edge([0, 8]), CubeFactor::Vertex(1), CubeFactor::Vertex(6), CubeFactor::Vertex(7) ] ) 
            + Cube( vec![ CubeFactor::Edge([0, 8]), CubeFactor::Vertex(1), CubeFactor::Vertex(7), CubeFactor::Vertex(9) ] );
        assert_eq!(boundary, ans);
    }

    #[test]
    fn flow_test() {
        let graph = graphs_for_tests::k_5_5_with_4_marked_points();

        // test 1
        let critical_cell = Cube( vec![CubeFactor::Vertex(6), CubeFactor::Vertex(7), CubeFactor::Vertex(9), CubeFactor::Vertex(24)] ).flow( &graph );
        let ans = 1*Cube( vec![CubeFactor::Vertex(0), CubeFactor::Vertex(1), CubeFactor::Vertex(2), CubeFactor::Vertex(3)] );
        assert_eq!(critical_cell, ans);

        // test 2
        let critical_cell = Cube( vec![CubeFactor::Edge([6, 9]), CubeFactor::Vertex(7), CubeFactor::Vertex(19), CubeFactor::Vertex(23)] ).flow( &graph );
        let ans = 1*Cube( vec![CubeFactor::Edge([6,9]), CubeFactor::Vertex(7), CubeFactor::Vertex(10), CubeFactor::Vertex(11)] );
        assert_eq!(critical_cell, ans);

        // test 3 (a non order-respecting edge)
        let critical_cell = Cube( vec![ CubeFactor::Edge([0, 8]), CubeFactor::Vertex(1), CubeFactor::Vertex(2), CubeFactor::Vertex(9) ] ).flow( &graph );
        let critical_cell: Vector<i64,_> = critical_cell.into_iter().map( |(i, cube)| i*cube ).sum();
        let ans = Cube( vec![ CubeFactor::Edge([0, 8]), CubeFactor::Vertex(1), CubeFactor::Vertex(2), CubeFactor::Vertex(3) ] )
            + Cube( vec![ CubeFactor::Vertex(0), CubeFactor::Vertex(1), CubeFactor::Edge([6,9]), CubeFactor::Vertex(7) ] );
        assert_eq!(critical_cell, ans);
    }

    #[test]
    fn boundary_test() {
        let graph = graphs_for_tests::k_5_5_with_4_marked_points();

        // test 1
        let two_cell = Cube( vec![CubeFactor::Edge([3, 13]), CubeFactor::Edge([6, 9]), CubeFactor::Vertex(7), CubeFactor::Vertex(10)] );
        let boundary = two_cell.get_boundary();
        let boundary = boundary.into_iter().map(|(i, x)| x.flow(&graph) * i ).sum::<Vector<i64,_>>();
        let boundary = boundary.into_iter().map(|(i, x)| x.get_boundary() * i ).sum::<Vector<i64,_>>();
        let mut boundary = boundary.into_iter().map(|(i, x)| x.flow(&graph) * i ).sum::<Vector<i64,_>>();
        boundary.iter_mut().for_each(|(_, cube)| cube.sort() );
        // let boundary: Vector<i64,_> = boundary[];
        assert_eq!(boundary, Vector::zero() );

        // test 2
        let two_cell = Cube( vec![ CubeFactor::Edge([0, 24]), CubeFactor::Edge([3, 13]), CubeFactor::Edge([6, 9]), CubeFactor::Vertex(7) ] );
        let boundary = two_cell.get_boundary();
        let boundary = boundary.into_iter().map(|(i, x)| x.flow(&graph) * i ).sum::<Vector<i64,_>>();
        let boundary = boundary.into_iter().map(|(i, x)| x.get_boundary() * i ).sum::<Vector<i64,_>>();
        let boundary = boundary.into_iter().map(|(i, x)| x.flow(&graph) * i ).sum::<Vector<i64,_>>();
        assert_eq!(boundary, Vector::zero() );


        use crate::graph;
        let cube = graph! {
            ["v0", "v1"],
            ["v1", "v2"],
            ["v2", "v3"],
            ["v3", "v0"],
            ["v4", "v5"],
            ["v5", "v6"],
            ["v6", "v7"],
            ["v7", "v4"],
            ["v0", "v4"],
            ["v1", "v5"],
            ["v2", "v6"],
            ["v3", "v7"]
        };
        let graph = SimpleGraph::subdivided_graph_from(&cube, 4);

        // test 3
        // {[0 11], [27 30], 28, 29} - {[0 11], [22 25], 23, 26}
        let two_cell = Cube( vec![CubeFactor::Edge([0, 11]), CubeFactor::Edge([27, 30]), CubeFactor::Vertex(28), CubeFactor::Vertex(29)] )
            - Cube( vec![CubeFactor::Edge([0, 11]), CubeFactor::Edge([22, 25]), CubeFactor::Vertex(23), CubeFactor::Vertex(26)] );
        let mut boundary = two_cell.into_iter().map(|(i, x)| x.get_boundary() * i ).sum::<Vector<i64,_>>();
        boundary = boundary.into_iter().map(|(i, x)| x.flow(&graph) * i ).sum::<Vector<i64,_>>();
        boundary.iter_mut().for_each(|(_, cube)| cube.sort() );
        assert_eq!(boundary, Vector::zero() );

        // test 4
        // - {[0 11], 1, [17 20], 18} - {[0 11], [17 20], 18, 21} - {[9 12], 10, [17 20], 18}
        let two_cell = -Cube( vec![CubeFactor::Edge([0, 11]), CubeFactor::Vertex(1), CubeFactor::Edge([17, 20]), CubeFactor::Vertex(18)] )
            - Cube( vec![CubeFactor::Edge([0, 11]), CubeFactor::Edge([17, 20]), CubeFactor::Vertex(18), CubeFactor::Vertex(21)] )
            - Cube( vec![CubeFactor::Edge([9, 12]), CubeFactor::Vertex(10), CubeFactor::Edge([17, 20]), CubeFactor::Vertex(18)] );
        let mut boundary = two_cell.into_iter().map(|(i, x)| x.get_boundary() * i ).sum::<Vector<i64,_>>();
        boundary = boundary.into_iter().map(|(i, x)| x.flow(&graph) * i ).sum::<Vector<i64,_>>();
        boundary.iter_mut().for_each(|(_, cube)| cube.sort() );
        assert_eq!(boundary, Vector::zero() );

        // test 6
        // {[3 24], 4, [14 31], 15}
        // let two_cell = 
        //     Cube( vec![CubeFactor::Edge(0, 19), CubeFactor::Edge(6, 29), CubeFactor::Vertex(7), CubeFactor::Edge(14, 31)] }
        //     + Cube( vec![CubeFactor::Edge(3, 24), CubeFactor::Edge(6, 29), CubeFactor::Vertex(7), CubeFactor::Edge(14, 31)] }
        //     - Cube( vec![CubeFactor::Edge(3, 24), CubeFactor::Vertex(4), CubeFactor::Edge(6, 29), CubeFactor::Edge(14, 31)] }
        //     - Cube( vec![CubeFactor::Edge(0, 19), CubeFactor::Vertex(1), CubeFactor::Edge(6, 29), CubeFactor::Edge(14, 31)] }
        //     + Cube( vec![CubeFactor::Edge(0, 11), CubeFactor::Edge(6, 29), CubeFactor::Edge(17, 20), CubeFactor::Vertex(18)] }
        //     + Cube( vec![CubeFactor::Edge(0, 11), CubeFactor::Edge(6, 29), CubeFactor::Edge(22, 25), CubeFactor::Vertex(23)] }
        //     - Cube( vec![CubeFactor::Edge(0, 11), CubeFactor::Edge(14, 31), CubeFactor::Edge(17, 20), CubeFactor::Vertex(18)] }
        //     - Cube( vec![CubeFactor::Edge(0, 11), CubeFactor::Edge(14, 31), CubeFactor::Edge(22, 25), CubeFactor::Vertex(23)] }
        //     - Cube( vec![CubeFactor::Edge(0, 19), CubeFactor::Edge(9, 12), CubeFactor::Vertex(10), CubeFactor::Edge(14, 31)] }
        //     + Cube( vec![CubeFactor::Edge(3, 24), CubeFactor::Edge(9, 12), CubeFactor::Vertex(10), CubeFactor::Edge(14, 31)] };
        // let mut boundary = two_cell.into_iter().map(|(i, x)| i * x.get_boundary() ).sum::<Vector<i64,_>>();
        // boundary = boundary.into_iter().map(|(i, x)| i * x.flow(&graph) ).sum::<Vector<i64,_>>();
        // println!("{}", boundary);
        // assert_eq!(boundary, Vector::zero() );

        let two_cell = 
            1*Cube( vec![CubeFactor::Edge([9, 12]), CubeFactor::Vertex(10), CubeFactor::Edge([22, 25]), CubeFactor::Vertex(23)] );
        let mut boundary = two_cell.into_iter().map(|(i, x)| x.get_boundary() * i ).sum::<Vector<i64,_>>();
        boundary = boundary.into_iter().map(|(i, x)| x.flow(&graph) * i ).sum::<Vector<i64,_>>();
        println!("{}", boundary);
        assert_eq!(boundary, Vector::zero() );
    }
}


pub struct UdnMorseCplx {
    pub factor_pool: SimpleGraph, // the base graph. subdivided in 'fn new(..)'
    data: Vec<Branch>,
    first_branch_indeces: Vec<usize>,
    final_branch_indeces: Vec<usize>,
    dim: usize,
    ranks_of_chain_cplx: Vec<usize>,
}

#[allow(unused)]
impl UdnMorseCplx {
    pub fn new(graph: &Graph<i64>, n_marked_points: usize, upper_dim: Option<usize>) -> Self {
        let factor_pool = SimpleGraph::subdivided_graph_from(graph, n_marked_points);
        let upper_dim = if let Some(d) = upper_dim {
            d
        } else {
            usize::MAX
        };
        // println!("{:?}", *factor_pool);
        Self::new_from_simple_graph(factor_pool, n_marked_points, upper_dim)
    }

    fn new_from_simple_graph(factor_pool: SimpleGraph, n_marked_points: usize, upper_dim: usize) -> Self {
        // if graph is empty or consisting of only one vertex, then we panic. TODO: handle these boundary cases.
        if factor_pool.n_vertices() < 2 { panic!() };

        let data = Vec::new();
        // Create the initial possible factors, namely, edges, 
        // - initializing possible factors
        let mut possible_factors = LabeledSet::new(
            factor_pool.iter().collect::<Vec<_>>()
        );
        // - adding initial possible factors
        // -- the only possibility for initial 0-dim'l factor is the vertex 0 
        possible_factors.add( CubeFactor::Vertex(0) );
        // -- initial 1-dim'l factors are the non-order-respecting edges
        factor_pool
            .edge_iter()
            .filter( |&factor| !factor.is_order_respecting() )
            .for_each( |factor| possible_factors.add(factor) );

        
        let first_branch_indeces: Vec<_> = (0..data.len()).collect();

        let mut cplx = UdnMorseCplx {
            factor_pool: factor_pool,
            data: data,
            first_branch_indeces: first_branch_indeces.clone(),
            final_branch_indeces: Vec::new(),
            dim: 0,
            ranks_of_chain_cplx: vec![0; n_marked_points+1],
        };

        cplx.expand(0, n_marked_points, possible_factors, upper_dim);

        cplx
    }

    fn expand(&mut self, branch_idx: usize, n_rem_factors: usize, mut possible_factors: LabeledSet<CubeFactor>, upper_dim: usize) -> LabeledSet<CubeFactor> {
        // if !self.data.is_empty() {
        //     println!( "so far chosen: {}", Cube::from_index(branch_idx, self) );
        // }
        // println!("possible_factors: {}", possible_factors);

        // possible factors are provided by upstairs, so just collect them and make new_branches
        let new_branches = if self.data.len()==0 {
            // if 'self.data.len()' is zero, then this is the initial call of this function
            possible_factors.iter().map(|&factor| Branch::new_parent(factor) ).collect::<Vec<_>>()
        } else {
            // if this is not the initial call, then we need to take care of the upper bound and lower bound when collecting factors
            let unpaid_edges = {
                let vertices:Vec<_> = self.factor_iter(branch_idx).filter(|&factor| factor.is_vertex() ).collect();
                self.factor_iter(branch_idx)
                    .filter(|&factor| factor.is_edge() ) // look only at the edges in the cell
                    .filter(|&factor| factor.is_in_maximal_tree(&self.factor_pool) )
                    .filter(|&factor| !vertices.iter().any(  |&v| v.is_disrespected_by(factor) ) ) // remove all the edges that are "paid back", i.e., the ones with a disrespected vertex
                    .collect::<Vec<_>>()
            };
            let upper_bound = {
                // upper bound of new factor can be determined by the minimum of the two numbers
                // 1. the uppper bound determined by the unpaid edges.
                let upper_bound1 = if let Some(&factor) = unpaid_edges.first() {
                    factor.initial()
                } else {
                    CubeFactor::Vertex( usize::MAX )
                };

                // 2. the uppper bound determined by the number of remaining factors. This guerantees that the branch has children in downstairs
                let upper_bound2 = {
                    let mut upper_bound2 = CubeFactor::Vertex(self.factor_pool.n_vertices());
                    (0..n_rem_factors-1).for_each(|_| {
                        upper_bound2 = CubeFactor::Vertex( upper_bound2.vertex()-1 ); 
                        while self.factor_iter(branch_idx).any(|factor| !factor.is_disjoint_from( upper_bound2 )) {
                            upper_bound2 = CubeFactor::Vertex( upper_bound2.vertex()-1 ); 
                        }
                    });
                    upper_bound2
                };

                std::cmp::min(upper_bound1, upper_bound2)
            };

            let mut new_branches = possible_factors.iter()
                .filter(|&&factor| self.factor_iter(branch_idx).all(|prev_factor| factor.is_disjoint_from(prev_factor) ))
                .filter(|&&factor| factor > *self.data[branch_idx] ) // ignore lower bounds
                .filter(|&&factor| factor < upper_bound || ( factor.is_edge() && factor.terminal() < upper_bound ) ) // ignore upper bounds 
                .filter(|&&factor| {
                    if factor.is_edge() && factor.is_in_maximal_tree(&self.factor_pool) {
                        let vertices_candidates: Vec<_> = (factor.terminal().vertex()+1..factor.initial().vertex()).map(|v| CubeFactor::Vertex(v)).filter(|&v| self.factor_iter(branch_idx).all(|f| f.is_disjoint_from(v))).collect();
                        !vertices_candidates.is_empty()
                    } else {
                        true
                    }
                }) 
                .map(|&factor| Branch::new(factor, branch_idx ))
                .collect::<Vec<_>>();

            // further, if the 'n_rem_factor' equals the number of debts, then we need to pay back right now, so remove all the factors that do not pay for unpaid edges.
            let n_debts = unpaid_edges.len();
            if n_rem_factors == n_debts {
                new_branches.retain(|b| (**b).is_vertex() );
            }

            // Similarly if the 'n_rem_factors' is 1, then we cannot add a non-order-respecting edge that is in the maximal tree anymore.
            if n_rem_factors == 1 {
                new_branches.retain(|b| (**b).is_vertex() || !(**b).is_in_maximal_tree(&self.factor_pool) );
            }

            // Finally, if the dimention of the cell is already the upper dim, then we cannot add edges. Hence, we will remove them.
            if self.dim_of(branch_idx) == upper_dim {
                new_branches.retain(|b| (**b).is_vertex() );
            }
            
            new_branches
        };
        let new_branch_indeces: Vec<_> = (0..new_branches.len()).map( |i| i+self.data.len()).collect();
        if new_branch_indeces.is_empty() { return possible_factors; }


        // put new indeces to the appropriate place. again, if 'self.data.len()' is zero, then this is the initial call of this function
        if self.data.len()==0 { 
            self.first_branch_indeces = new_branch_indeces.clone();
        } else {
            self.data[branch_idx].child_indeces = new_branch_indeces.clone();
        };

        self.data.append(&mut new_branches.clone());

        // base case
        if n_rem_factors == 1 {
            self.final_branch_indeces.append( &mut new_branch_indeces.clone() );
            // modify the rank of the complex
            let dim_of_cells: Vec<_> = new_branch_indeces.iter().map(|&idx| self.dim_of(idx)).collect();
            dim_of_cells.iter().for_each( |&dim| self.ranks_of_chain_cplx[dim] += 1 );
            // mofify the dimention of the complex if necessary
            let max_dim_of_cells = dim_of_cells.into_iter().max().unwrap();
            if max_dim_of_cells > self.dim { self.dim = max_dim_of_cells; }
            return possible_factors;
        }

        // for each branch...
        for (idx, branch) in new_branch_indeces.into_iter().zip( new_branches.iter() ) {
            let curr_factor = **branch;
            // record the factors added to 'possible_factors' in order to recover it later.
            let mut factors_added = Vec::new(); 

            // modify 'possible_factors' for the factors downstairs.
            if curr_factor.is_vertex() { 
                let new_vertex = curr_factor.vertex()+1;
                if new_vertex < self.factor_pool.n_vertices() { factors_added.push( CubeFactor::Vertex( curr_factor.vertex()+1 ) ); }
            } else if curr_factor.is_in_maximal_tree(&self.factor_pool) {
                // Coming here means that 'curr_factor' is an edge in the maximal tree.
                let [_, j] = curr_factor.edge(); // get two ends of the edge
                // add vertices (no edges will be added)
                self.factor_pool
                    .vertex_iter()
                    .filter(|&v| v > curr_factor.terminal())
                    .filter(|&v| v.flow(&self.factor_pool).terminal() == curr_factor.terminal() )
                    .for_each(|factor| factors_added.push(factor) ); 
                if j+1<self.factor_pool.n_vertices() { factors_added.push( CubeFactor::Vertex(j+1) )}; // add a vertex blocked by the initial vertex of the edge
            } else {
                // Coming here means that 'curr_factor' is an edge not in the maximal tree.
                let [_, j] = curr_factor.edge(); // get two ends of the edge
                // add vertices (no edges will be added)
                self.factor_pool
                    .vertex_iter()
                    .filter(|&v| v > curr_factor.terminal())
                    .filter( |&v| v.flow(&self.factor_pool).terminal() == curr_factor.terminal() )
                    .for_each(|factor| factors_added.push(factor) ); 
            }

            factors_added.iter().for_each( |&factor| possible_factors.add( factor ));

            // go downstairs
            possible_factors = self.expand(idx, n_rem_factors-1, possible_factors, upper_dim);

            // recover the 'possible_factors' for other.
            possible_factors.add( curr_factor );
            factors_added.iter().for_each(|&factor| possible_factors.remove( factor ));
        }

        possible_factors
    }

    pub fn dim_of(&self, branch_idx: usize) -> usize {
        let dim_of_branch = match self.data[branch_idx].val {
            CubeFactor::Edge(_) => 1,
            CubeFactor::Vertex(_) => 0,
        };

        // return the recursive call
        dim_of_branch + if let Some(p_idx) = self.data[branch_idx].parent_idx {
            self.dim_of(p_idx)
        } else {
            0
        }
    }

    fn get_cell_from_idx(&self, idx: usize) -> Vec<CubeFactor> {
        let mut out: Vec<_> = self.factor_iter(idx).collect();
        out.reverse();
        out
    }

    fn get_index_of(&self, cell: &[CubeFactor]) -> usize {
        let mut branch_idx = match self.first_branch_indeces.iter().find(|&&idx| self.data[idx].val == cell[0]) {
            Some(&idx) => idx,
            None => { panic!("cell {:?} not found in the complex.", cell); },
        };
        for i in 1..cell.len() {
            branch_idx = match self.data[branch_idx].child_indeces.iter().find(|&&idx| self.data[idx].val == cell[i] ) {
                Some(&idx) => idx,
                None => { panic!("cell {:?} not found in the complex.", cell); },
            };
        };
        branch_idx
    }

    fn get_boundary_of(&self, branch_idx: usize) -> Vector<i64, Cube> {
        let boundary = Cube::from_index(branch_idx, self).get_boundary();
        let mut boundary = boundary.into_iter().map(|(i, cube)| cube.flow(&self.factor_pool) * i ).sum::<Vector<i64,_>>();
        boundary.iter_mut().for_each( |(_, cube)| cube.sort() );
        let boundary = boundary.into_iter().map(|(i, cube)| i*cube ).sum::<Vector<i64,_>>();

        boundary
    }

    pub fn get_cells_of_dim(&self, d: usize) -> Vec<Cube> {
        self.final_branch_indeces.iter()
            .filter(|&&i| self.dim_of(i) == d )
            .map(|&i| Cube::from_index(i, &self) )
            .collect()
    }

    // pub fn get_edge_path_of_edge_in_maximal_tree(&self, c: Cube) -> Vec<Vec<[usize; 2]>> {
    //     assert!(
    //         c.0.iter().filter(|f| f.is_edge()).count() == 1,
    //         "The input must be a critical edge, but it is a vertex."
    //     );

    //     // get two ends of c
    //     let mut non_order_respecting_edge = c.0.iter().find(|&&f| f.is_edge() ).unwrap().edge();
    //     non_order_respecting_edge.swap(0, 1); // Now the second element is smaller than the first.

    //     let mut path = VecDeque::from([vec![non_order_respecting_edge]]);

    //     let first_snake_len = c.iter()
    //         .filter(|f| f.is_vertex())
    //         .filter(|f| (non_order_respecting_edge[1]..non_order_respecting_edge[0]).contains(&f.vertex()))
    //         .count();

    //     let first_snake_head = c.iter()
    //         .filter(|f| f.is_vertex())
    //         .map(|f| f.vertex() )
    //         .filter(|v| (non_order_respecting_edge[1]..non_order_respecting_edge[0]).contains(&v))
    //         .max()
    //         .unwrap();

    //     let second_snake_len = c.iter()
    //         .filter(|f| f.is_vertex())
    //         .filter(|f| f.vertex() > non_order_respecting_edge[0])
    //         .count();
    //     let second_snake_head = non_order_respecting_edge[0] + second_snake_len;

    //     let stationary_vertices_len = c.len()-first_snake_len - second_snake_len - 1;

    //     // add some motions to the non-order-respecting motion
    //     path.front_mut().unwrap().append(
    //         &mut (non_order_respecting_edge[0]..non_order_respecting_edge[0]+second_snake_len).map(|i| [i+1, i] ).collect::<Vec<_>>()
    //     );

    //     // generating motions before the non-order-respecting motion
    //     let mut outward_path = VecDeque::new();
    //     if first_snake_len != 0 {
    //         outward_path = self.generate_path::<false>(first_snake_head, first_snake_len, stationary_vertices_len);
    //     }
    //     let mut path_of_second_snake = self.generate_path::<false>(second_snake_head, second_snake_len+1, first_snake_len + stationary_vertices_len);
    //     for (p1, p2) in outward_path.iter_mut().rev().skip(first_snake_len).zip(path_of_second_snake.iter_mut().rev()) {
    //         p1.append(p2);
    //     }
    //     path.front_mut().unwrap().append(&mut outward_path.pop_back().unwrap());
    //     outward_path.append(&mut path);
    //     path = outward_path;

    //     // generating motions after the non-order-respecting motion
    //     let mut return_path = self.generate_path::<true>(first_snake_head, first_snake_len+1, stationary_vertices_len);
    //     if second_snake_len != 0 { 
    //         path_of_second_snake = self.generate_path::<true>(second_snake_head-1, second_snake_len, first_snake_len + stationary_vertices_len + 1);
    //     }
    //     for (p1, p2) in return_path.iter_mut().skip(first_snake_len).zip(path_of_second_snake.iter_mut()) {
    //         p1.append(p2);
    //     }
    //     path.append(&mut return_path);

    //     path
    // }

    // pub fn get_edge_path_of_edge_not_in_maximal_tree(&self, c: Cube) -> Vec<Vec<[usize; 2]>> {
    //     assert!(
    //         c.0.iter().filter(|f| f.is_edge()).count() == 1,
    //         "The input must be a critical edge, but it is a vertex."
    //     );

    //     let edge = c.0.iter().find(|&&f| f.is_edge() ).unwrap().edge();
    //     let path = VecDeque::new();
    //     path
    // }

    fn generate_path<const FORWARD: bool>(&self, snake_head: usize, snake_len: usize, flow_terminal: usize) -> VecDeque<Vec<[usize; 2]>>{

        assert!(snake_len!=0);
        let mut path = VecDeque::new();
        let mut path_of_last_vertex = Vec::new();
        let mut curr = snake_head;
        while curr > flow_terminal {
            // flow to the terminal
            let terminal = (0..curr).rev()
                .find(|&i| self.factor_pool.contains(i, curr) && self.factor_pool.maximal_tree_contains([i, curr]) )
                .unwrap();

            if FORWARD {
                path_of_last_vertex.push( [curr, terminal] );
            } else {
                path_of_last_vertex.push( [terminal, curr] );
            }
            curr = terminal;
        }

        (0..path_of_last_vertex.len() - snake_len + 1)
            .map(|i| (i..i+snake_len).map(|j| path_of_last_vertex[j]).collect::<Vec<_>>() )
            .for_each(|edge| if FORWARD {
                    path.push_back(edge)
                } else {
                    path.push_front(edge)
                } );

        path
    }
}

pub struct UdnMorseFactorIter<'a> {
    cplx: &'a UdnMorseCplx,
    curr_idx: usize,
    is_done: bool,
}

impl<'a> Iterator for UdnMorseFactorIter<'a> {
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

impl<'a> UdnMorseCplx {
    pub fn factor_iter( &'a self, branch_idx: usize ) -> UdnMorseFactorIter<'a> {
        UdnMorseFactorIter {
            cplx: self,
            curr_idx: branch_idx,
            is_done: false,
        }
    }
}


impl fmt::Display for UdnMorseCplx {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for i in 0..=self.dim {
            write!(f, "{}-cells: {} cells in total", i, self.ranks_of_chain_cplx[i])?;
            for &cell_idx in self.final_branch_indeces.iter().filter(|&&idx| self.dim_of(idx) == i) {
                // printing each cell
                write!(f, "\n\t{{")?;
                for &factor in self.get_cell_from_idx(cell_idx).iter() {
                    // printping each factor
                    write!(f, "{}", factor)?;
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


#[derive(Debug, PartialEq, Eq, Clone)]
pub struct EdgePath (Vec<Vec<[usize; 2]>>);

impl EdgePath {
    pub fn trivial_path() -> Self {
        Self(Vec::new())
    }
    pub fn reduce_to_geodesic(self) -> Self {
        let mut out = vec![Vec::new()];

        let commute = |v: &[[usize; 2]], w: &[usize; 2]| -> bool {
            for i in v {
                if i[0]==w[1] || i[1]==w[0]|| i[1]==w[1] || i[0]==w[0] {
                    return false;
                }
            }
            true
        };

        let composable = |v: &[[usize; 2]], w: &[usize; 2]| -> bool {
            for i in v {
                if i[1]==w[0] || i[0]==w[1] || i[1]==w[1] || i[0]==w[0] {
                    return false;
                }
            }
            true
        };

        let cancelable = |v: &[usize; 2], w: &[usize; 2]| -> bool {
            if v[1]==w[0] && v[0]==w[1] {
                true
            } else {
                false
            }
        };

        let mut some_canceled = false;

        for step in self.0.into_iter().map(|f| f.into_iter()).flatten() {
            if let Some(i) = (0..out.len()).rev().find(|&i| !commute(&out[i], &step)) {
                if composable(&out[i], & step) {
                    out[i].push(step);
                } else if let Some(idx) = (0..out[i].len()).find(|&j| cancelable(&out[i][j], &step)) {
                    some_canceled = true;

                    out[i].remove(idx);
                    if out[i].is_empty() {
                        out.remove(i);
                    }
                } else if i+1 == out.len() {
                    out.push(vec![step]);
                } else {
                    out[i+1].push(step);
                }
            } else {
                out[0].push(step);
            };
        }

        if some_canceled {
            Self(out).reduce_to_geodesic()
        } else {
            Self(out)
        }
    }

    pub fn composed_with(mut self, mut other: Self) -> Self {
        self.0.append(&mut other.0);
        self
    }

    pub fn evaluate_permutation<const N: usize>(&self) -> ConstPermutation<N> {
        let mut out = {
            let mut out = [0; N];
            (0..N).for_each(|i| out[i]=i );
            out
        };

        for step in self.iter() {
            let mut next = out;
            for edge in step.iter() {
                let update_idx = (0..N).find(|&i| out[i] == edge[0]).unwrap_or_else(|| panic!("path={self:?}, \n step = {step:?},\n edge={edge:?}") );
                next[update_idx] = edge[1];
            }
            out = next;
        }

        ConstPermutation::from(out)
    }

    pub fn inverse(mut self) -> Self {
        // reverse each step
        for edge in self.0.iter_mut().map(|step| step.iter_mut()).flatten() {
            edge.swap(0, 1);
        }

        for step in self.0.iter_mut() {
            step.reverse();
        }

        let path = Self( self.0.into_iter().rev().collect::<Vec<_>>() );
        path.reduce_to_geodesic()
    }
}

impl std::ops::Deref for EdgePath {
    type Target = Vec<Vec<[usize; 2]>>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<Vec<Vec<[usize; 2]>>> for EdgePath {
    fn from(value: Vec<Vec<[usize; 2]>>) -> Self {
        EdgePath(value)
    }
}


#[cfg(test)]
mod udn_morse_cplx_test {
    use crate::cubical::*;
    use crate::graph;

    #[test]
    fn init_test1() {
        let graph = graphs_for_tests::k_5_5_with_4_marked_points();

        // create a morse complex with 4 marked points
        let complex = UdnMorseCplx::new_from_simple_graph(graph, 4, usize::MAX);

        // check the morse complex
        // - the 0-cell
        use crate::graph::CubeFactor::Vertex;
        let zero_cells = complex.final_branch_indeces.iter().filter(|&&idx| complex.dim_of(idx)==0 ).map(|&idx| complex.get_cell_from_idx(idx)).collect::<Vec<_>>();
        assert_eq!( zero_cells, vec![ vec![Vertex(0), Vertex(1), Vertex(2), Vertex(3)]] );

        // println!("{}", complex);
        assert!( complex.final_branch_indeces.iter().all(|&idx| Cube::from_index(idx, &complex).property(&complex.factor_pool)==MorseProp::Critical) );
    }


    #[test]
    fn init_test2() {
        let graph = graphs_for_tests::skeleton_of_cube_with_4_marked_points();

        // create a morse complex with 4 marked points
        let complex = UdnMorseCplx::new_from_simple_graph(graph, 4, usize::MAX);

        // check the morse complex
        // - the 0-cell
        use crate::graph::CubeFactor::Vertex;
        let zero_cells = complex.final_branch_indeces.iter().filter(|&&idx| complex.dim_of(idx)==0 ).map(|&idx| complex.get_cell_from_idx(idx)).collect::<Vec<_>>();
        assert_eq!( zero_cells, vec![ vec![Vertex(0), Vertex(1), Vertex(2), Vertex(3)]] );

        // println!("{}", complex);
        assert!( complex.final_branch_indeces.iter().all(|&idx| Cube::from_index(idx, &complex).property(&complex.factor_pool)==MorseProp::Critical) );

        assert_eq!( complex.ranks_of_chain_cplx, vec![1, 33, 109, 58, 2]);
        println!("{}", complex);
    }


    #[test]
    fn init_test3() {
        let graph = graph!{
            ["v1", "v2"],
            ["v1", "v3"],
            ["v1", "v4"],
            ["v2", "v3"],
            ["v2", "v4"],
            ["v3", "v4"]
        };
        let complex = UdnMorseCplx::new(&graph, 2, None);

        use crate::graph::CubeFactor::{ Vertex, Edge };

        // 0-cells
        let zero_cells = complex.final_branch_indeces.iter().filter(|&&idx| complex.dim_of(idx)==0 ).map(|&idx| complex.get_cell_from_idx(idx)).collect::<Vec<_>>();
        assert_eq!( zero_cells, vec![ vec![Vertex(0), Vertex(1)]] );

        // 1-cells
        let one_cells = complex.final_branch_indeces.iter().filter(|&&idx| complex.dim_of(idx)==1 ).map(|&idx| complex.get_cell_from_idx(idx)).collect::<Vec<_>>();
        assert_eq!( one_cells, vec![ vec![Vertex(0), Edge([2,8])], vec![Edge([0,5]), Vertex(1)], vec![Edge([4,6]), Vertex(5)], vec![Edge([2,8]), Vertex(3)], vec![Edge([7,9]), Vertex(8)], vec![Edge([0,9]), Vertex(1)] ] );


        assert!( complex.final_branch_indeces.iter().all(|&idx| Cube::from_index(idx, &complex).property(&complex.factor_pool)==MorseProp::Critical) );
        println!("{}", complex);
    }

    #[test]
    fn init_test4() {
        let graph = graph!{
            ["v0", "v1"],
            ["v1", "v1"]
        };
        let complex = UdnMorseCplx::new(&graph, 3, None);

        use crate::graph::CubeFactor::{ Vertex, Edge };

        // 0-cells
        let zero_cells = complex.final_branch_indeces.iter().filter(|&&idx| complex.dim_of(idx)==0 ).map(|&idx| complex.get_cell_from_idx(idx)).collect::<Vec<_>>();
        assert_eq!( zero_cells, vec![ vec![Vertex(0), Vertex(1), Vertex(2)]] );

        // 1-cells
        let one_cells = complex.final_branch_indeces.iter().filter(|&&idx| complex.dim_of(idx)==1 ).map(|&idx| complex.get_cell_from_idx(idx)).collect::<Vec<_>>();
        assert_eq!( one_cells, vec![ vec![Vertex(0), Vertex(1), Edge([2,5])], vec![Vertex(0), Edge([2,5]), Vertex(3)], vec![Edge([2,5]), Vertex(3), Vertex(4)] ] );

        // 2-cells
        let two_cells = complex.final_branch_indeces.iter().filter(|&&idx| complex.dim_of(idx)==2 ).map(|&idx| complex.get_cell_from_idx(idx)).collect::<Vec<_>>();
        assert!( two_cells.is_empty() );


        assert!( complex.final_branch_indeces.iter().all(|&idx| Cube::from_index(idx, &complex).property(&complex.factor_pool)==MorseProp::Critical) );
        println!("{}", complex);
    }
}



impl Complex for UdnMorseCplx {
    type Cell = Cube;
    fn boundary_map (&self) -> Vec< (Matrix<i64>, Vec<Self::Cell>) > {
        let mut numbering = vec![0; self.data.len()];
        let mut counter = vec![0; self.dim+1];
        self.final_branch_indeces.iter().for_each(|&idx| {
            numbering[idx] = counter[self.dim_of(idx)];
            counter[self.dim_of(idx)] += 1;
        });

        // instantiating the boundary maps
        let mut boundary_maps = vec![Matrix::zero(1, self.ranks_of_chain_cplx[0])];
        (1..=self.dim).for_each(|i| boundary_maps.push(Matrix::zero(self.ranks_of_chain_cplx[i-1], self.ranks_of_chain_cplx[i])));
        let mut basis: Vec<Vec<Self::Cell>> = vec![Vec::new(); self.dim+1];
        let zero_cell_idx = *self.final_branch_indeces.iter().find(|&&i| self.dim_of(i) == 0).unwrap();
        basis[0].push( Cube::from_index(zero_cell_idx, self) );

        // filling up the boundary maps
        for &i in self.final_branch_indeces.iter().filter(|&&i| self.dim_of(i) > 0 ) {
            let boundary = self.get_boundary_of(i);
            let dim = self.dim_of(i);
            // encode the boundary information
            boundary.into_iter().for_each(|(coeff, cube)| {
                let coeff = coeff as i64;
                let j = self.get_index_of(&cube);
                let val = boundary_maps[dim][( numbering[j], numbering[i] )];
                boundary_maps[dim][( numbering[j],numbering[i] )] = val+coeff;
            });

            // add the cell to the basis
            basis[dim].push( Cube::from_index(i, self) );
        }

        let boundary_maps = boundary_maps.into_iter().zip(basis.into_iter()).collect::<Vec<_>>();
        debug_assert!(boundary_maps.iter().all(|(map, basis)| map.size().1 == basis.len() ), "{:?}:  {:?}", 
            boundary_maps.iter().map(|(map, _)| map.size()).collect::<Vec<_>>(),
            boundary_maps.iter().map(|(_, basis)| basis.len()).collect::<Vec<_>>()
        );
        boundary_maps
    }
}

#[cfg(test)]
mod udn_morse_homology_test {
    use crate::graph;
    use crate::cubical::UdnMorseCplx;
    use crate::IntegralHomology;
    
    // #[test]
    fn homology_test() {
        // test 1
        let graph = graph!{
            ["v0", "v1"],
            ["v0", "v2"],
            ["v0", "v3"]
        };
        let cplx = UdnMorseCplx::new(&graph, 2, None);
        assert_eq!( H!(cplx).to_string(), "dim 0: Z,\ndim 1: Z,\n" );

        // test 2
        let cplx = UdnMorseCplx::new(&graph, 3, None);
        assert_eq!( H!(cplx).to_string(), "dim 0: Z,\ndim 1: Z^3,\n" );

        // test 3
        let graph = graph!{
            ["v0", "v1"],
            ["v1", "v1"]
        };
        let cplx = UdnMorseCplx::new(&graph, 5, None);
        assert_eq!( H!(cplx).to_string(), "dim 0: Z,\ndim 1: Z^5,\n" );

        // test 4
        let graph = graph!{
            ["v1", "v2"],
            ["v1", "v3"],
            ["v1", "v4"],
            ["v2", "v3"],
            ["v2", "v4"],
            ["v3", "v4"]
        };
        let cplx = UdnMorseCplx::new(&graph, 2, None);
        assert_eq!( H!(cplx).to_string(), "dim 0: Z,\ndim 1: Z^4,\n" );
        let cplx = UdnMorseCplx::new(&graph, 3, None);
        println!( "{:?}", cplx.ranks_of_chain_cplx );
        println!( "{}", cplx );
        assert_eq!( H!(cplx).to_string(), "dim 0: Z,\ndim 1: Z^4,\ndim 2: Z^3,\n" );

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
        let cplx = UdnMorseCplx::new(&graph, 2, None);
        println!("{cplx}");
        assert_eq!( H!(cplx).to_string(), "dim 0: Z,\ndim 1: Z^6 + Z_2,\n" );
        let cplx = UdnMorseCplx::new(&graph, 3, None);
        assert_eq!( H!(cplx).to_string(), "dim 0: Z,\ndim 1: Z^6 + Z_2,\ndim 2: Z^30,\n" );

        use super::graphs_for_tests;
        let graph = graphs_for_tests::k_5_5_with_4_marked_points();
        let complex = UdnMorseCplx::new_from_simple_graph(graph, 3, usize::MAX);
        assert_eq!( H!(complex).to_string(), "dim 0: Z,\ndim 1: Z^6 + Z_2,\ndim 2: Z^30,\n" );
    }
}
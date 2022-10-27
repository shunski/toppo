use crate::Complex;
use crate::graph::*;
use algebra::lin_alg::matrix::Matrix;
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
    pub fn udc(graph: &Graph, n_marked_points: usize) -> CubicalCplx {
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

    pub fn dim_of(&self, branch_idx: usize) -> usize {
        let dim_of_branch = match self.data[branch_idx].val {
            CubeFactor::Edge(_,_) => 1,
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
    fn boundary_map (&self) -> Vec<Matrix<i64>> {
        let mut numbering = vec![0; self.data.len()];
        let mut counter = vec![0; self.dim+1];
        self.final_branch_indeces.iter().for_each(|&idx| {
            numbering[idx] = counter[self.dim_of(idx)];
            counter[self.dim_of(idx)] += 1;
        });

        // debug_assert!();

        // instantiating the boundary maps
        let mut boundary_maps = vec![Matrix::zero(1, self.ranks_of_chain_cplx[0])];
        (1..=self.dim).for_each(|i| boundary_maps.push(Matrix::zero(self.ranks_of_chain_cplx[i-1], self.ranks_of_chain_cplx[i])));

        // filling up the boundary maps
        for &idx in self.final_branch_indeces.iter().filter(|&&idx| self.dim_of(idx) > 0 ) {
            let boundary_indeces = self.get_boundaries_of(idx);
            let dim = self.dim_of(idx);
            boundary_indeces.iter().for_each(|&(sign, b_idx)| {
                let val = boundary_maps[dim].get(numbering[b_idx], numbering[idx]);
                boundary_maps[dim].write(numbering[b_idx],numbering[idx], if sign {val+1} else {val-1})
            });
        }

        boundary_maps
    }
}


#[cfg(test)]
mod homology_test {
    use crate::graph;
    use crate::cubical::CubicalCplx;
    use crate::IntegralHomology;
    
    #[test]
    fn homology_test() {
        // test 1
        let graph = graph!{
            ["v0", "v1"],
            ["v0", "v2"],
            ["v0", "v3"]
        };
        let cplx = CubicalCplx::udc(&graph, 2);
        assert_eq!( H!(cplx), "dim 0: Z,\ndim 1: Z,\n" );

        // test 2
        let cplx = CubicalCplx::udc(&graph, 3);
        assert_eq!( H!(cplx), "dim 0: Z,\ndim 1: Z^3,\n" );

        // // test 3
        // let graph = graph!{
        //     ["v0", "v1"],
        //     ["v1", "v1"]
        // };
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
    }
}


struct UdnMorseCplx {
    factor_pool: SimpleGraph,
    data: Vec<Branch>,
    first_branch_indeces: Vec<usize>,
    final_branch_indeces: Vec<usize>,
    dim: usize,
    ranks_of_chain_cplx: Vec<usize>,
}

#[allow(unused)]
impl UdnMorseCplx {
    pub fn new(graph: &Graph, n_marked_points: usize) -> Self {
        let factor_pool = SimpleGraph::subdivided_graph_from(graph, n_marked_points);
        Self::new_from_simple_graph(factor_pool, n_marked_points)
    }

    fn new_from_simple_graph(factor_pool: SimpleGraph, n_marked_points: usize) -> Self {
        // if graph is empty or consisting of only one vertex, then we panic. TODO: handle these boundary cases.
        if factor_pool.n_vertices() < 2 { panic!() };

        let data = Vec::new();
        // Create the initial possible factors, namely, edges, 
        // - initializing possible factors
        let mut possible_factors = LabeledSet::new(
            factor_pool.iter().collect::<Vec<_>>()
        );
        factor_pool
            .iter()
            .filter(|factor| factor.is_vertex() || factor.is_in_maximal_tree( &factor_pool ) )
            .for_each(|factor| possible_factors.add(factor) );
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

        cplx.expand(0, n_marked_points, possible_factors);

        cplx
    }

    fn expand(&mut self, branch_idx: usize, n_rem_factors: usize, mut possible_factors: LabeledSet<CubeFactor>) -> LabeledSet<CubeFactor> {
        println!("{}",n_rem_factors);
        // possible factors are provided by the people upstairs, so just collect them
        let new_branches = if self.data.len()==0 {
            // if 'self.data.len()' is zero, then this is the initial call of this function
            possible_factors.iter().map(|&factor| Branch::new_parent(factor) ).collect::<Vec<_>>()
        } else {
            // if this is not the initial call, then we need to take care of the upper bound and lower bound when collecting factors
            let unpaid_edges = {
                let vertices:Vec<_> = self.factor_iter(branch_idx).filter(|&factor| factor.is_vertex() ).map(|factor| factor.vertex()).collect();
                self.factor_iter(branch_idx)
                    .filter(|&factor| factor.is_edge() )
                    .filter(|&factor| factor.is_in_maximal_tree(&self.factor_pool) )
                    .map(|factor| factor.edge()) // now we have the iterator over the non-order-respecting edges in the maximal tree
                    .filter(|&(i, _)| !vertices.iter().any(  |&v| self.factor_pool.contains(i, v) ) ) // remove all the edges that are "paid back", i.e., the ones with a disrespected vertex
                    .collect::<Vec<_>>()
            };
            let upper_bound = if let Some(&(_, j)) = unpaid_edges.first() {
                CubeFactor::Vertex( j )
            } else {
                CubeFactor::Vertex( usize::MAX )
            };
            let mut new_branches = possible_factors.iter()
                .skip_while(|&&factor| factor > self.factor_iter(branch_idx).next().unwrap() ) // lower bound when collecting factors
                .take_while(|&&factor| factor < upper_bound ) // upper bound when collecting factors
                .map(|&factor| Branch::new(factor, branch_idx ))
                .collect::<Vec<_>>();

            // further, if the 'n_rem_factor' equals the number of debts, then we need to pay back right now, so remove all the cells that do not pay for unpaid edges.
            let n_debts = unpaid_edges.len();
            if n_rem_factors == n_debts {
                new_branches.retain(|b| b.val.is_vertex() );
            }
            new_branches
        };
        let new_branch_indeces: Vec<_> = (0..new_branches.len()).map( |i| i+self.data.len()).collect();


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
                let new_factor = CubeFactor::Vertex( curr_factor.vertex()+1 );
                possible_factors.add( new_factor );
                factors_added.push( new_factor );
            } else if curr_factor.is_in_maximal_tree(&self.factor_pool) {
                // Coming here means that 'curr_factor' is an edge in the maximal tree.
                let (i, j) = curr_factor.edge(); // get two ends of the edge
                // add vertices (no edges will be added)
                self.factor_pool
                    .vertex_iter()
                    .skip_while(|&v| v <= curr_factor.terminal())
                    .take_while(|&v| v < curr_factor.initial() )
                    .filter( |&k| k.is_disrespected_by(curr_factor) )
                    .for_each(|factor| factors_added.push(factor) ); 
                factors_added.push( CubeFactor::Vertex(j+1) ); // add a vertex blocked by the initial vertex of the edge
            } else {
                // Coming here means that 'curr_factor' is an edge not in the maximal tree.
                let (i, j) = curr_factor.edge(); // get two ends of the edge
                // add vertices (no edges will be added)
                self.factor_pool
                    .vertex_iter()
                    .skip_while(|&v| v <= curr_factor.terminal() )
                    .take_while(|&v| v < curr_factor.initial() )
                    .filter( |&v| v.is_disrespected_by(curr_factor) )
                    .for_each(|factor| factors_added.push(factor) ); 
            }

            factors_added.iter().for_each( |&factor| possible_factors.add( factor ));

            // go downstairs
            possible_factors = self.expand(idx, n_rem_factors-1, possible_factors);

            // recover the 'possible_factors' for other.
            possible_factors.add( curr_factor );
            factors_added.iter().for_each(|&factor| possible_factors.remove( factor ));
        }

        possible_factors
    }

    pub fn dim_of(&self, branch_idx: usize) -> usize {
        let dim_of_branch = match self.data[branch_idx].val {
            CubeFactor::Edge(_,_) => 1,
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

struct UdnMorseFactorIter<'a> {
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
    fn factor_iter( &'a self, branch_idx: usize ) -> UdnMorseFactorIter<'a> {
        UdnMorseFactorIter {
            cplx: self,
            curr_idx: branch_idx,
            is_done: false,
        }
    }
}


#[cfg(test)]
mod udn_morse_cplx_test {
    use crate::cubical::*;

    #[test]
    fn init_test() {
        // we manually create the graph of K_{5,5} with a maximal tree.
        let graph = {
            let mut raw_graph = RawSimpleGraph::new(25);
            // add deleted edges
            raw_graph.add_edge(0, 8);
            raw_graph.add_edge(0, 15);
            raw_graph.add_edge(0, 24);
            raw_graph.add_edge(3, 13);
            raw_graph.add_edge(3, 22);
            raw_graph.add_edge(6, 20);
            // add edges in the maximal tree
            raw_graph.add_edge(6, 9);
            raw_graph.add_edge(11, 14);
            raw_graph.add_edge(11, 16);
            raw_graph.add_edge(18, 21);
            raw_graph.add_edge(18, 23);
            (0..24).filter(|&i| i!=8 && i!=13 && i!=15 && i!=20 ).for_each(|i| raw_graph.add_edge(i, i+1));

            raw_graph.is_maximal_tree_chosen = true;

            SimpleGraph::without_lables_from(raw_graph)
        };

        // create a morse complex with 3 marked points
        let _morse_complex = UdnMorseCplx::new_from_simple_graph(graph, 3);

    }
}

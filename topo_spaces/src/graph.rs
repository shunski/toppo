use algebra::permutation;
use util::BinaryPool;
use algebra::commutative::rational::*;
use std::fmt;
use algebra::non_commutative::{permutation::*, Group};

#[allow(unused)]
#[derive(Clone)]
struct Edge {
    distance: f64,
}

impl Edge {
    fn new() -> Edge {
        Edge{
            distance: 1.0,
        }
    }
}


#[allow(unused)]
pub struct Graph {
    data: Vec<Vec<Vec<Edge>>>,
    vertices_labeling: Vec<String>,

    n_non_loop_edges: usize,
    n_non_loop_edges_up_to_multiplicity: usize,
    n_loops: usize,
}


#[allow(unused)]
impl Graph {
    pub fn from( edges: Vec<(String, String)> ) -> Graph {
        let mut n_non_loop_edges = edges.len();
        let mut n_non_loop_edges_up_to_multiplicity = edges.len();
        let mut n_loops = 0;

        let mut vertices_labeling = Vec::new();
        // registering vertices
        edges.iter().for_each( |(x, y)| {
            if !(vertices_labeling.contains(x)) {
                vertices_labeling.push(x.clone());
            } 
            if !(vertices_labeling.contains(y)) {
                vertices_labeling.push(y.clone());
            }
        });

        // registering edges
        let n_vertices = vertices_labeling.len();
        let mut data = vec![vec![Vec::new(); n_vertices]; n_vertices];
        edges.iter().for_each( |(v, w)| {
            let v_idx = vertices_labeling.iter().position(|x| x == v);
            let mut v_idx = match v_idx {
                Some(idx) => idx,
                None => panic!(),
            };

            let w_idx = vertices_labeling.iter().position(|x| x == w);
            let mut w_idx = match w_idx {
                Some(idx) => idx,
                None => panic!(),
            };

            // order in such a way that v_idx <= w_idx
            let (v_idx, w_idx) = if v_idx <= w_idx {
                (v_idx, w_idx)
            } else {
                (w_idx, v_idx)
            };

            data[v_idx][w_idx].push(Edge::new());

            // update some members
            if data[v_idx][w_idx].len() > 1 {
                n_non_loop_edges_up_to_multiplicity -= 1;
            } 
            if v_idx == w_idx {
                n_loops += 1;
                n_non_loop_edges -= 1;
                n_non_loop_edges_up_to_multiplicity -= 1;
            }
        });

        Graph {
            data: data,
            vertices_labeling: vertices_labeling,
            n_non_loop_edges: n_non_loop_edges,
            n_non_loop_edges_up_to_multiplicity: n_non_loop_edges_up_to_multiplicity,
            n_loops: n_loops,
        }
    }

    pub fn n_vertices(&self) -> usize { self.vertices_labeling.len() }

    pub fn n_non_loop_edges(&self) -> usize { 
        self.n_non_loop_edges
    }

    pub fn n_non_loop_edges_up_to_multiplicity(&self) -> usize {
        self.n_non_loop_edges_up_to_multiplicity
    }

    pub fn n_loops(&self) -> usize {
        self.n_loops
    }
}


#[macro_export]
macro_rules! graph {
    ( $([$v1:expr, $v2:expr]), *) => {
        {
            let mut v: Vec<(String, String)> = Vec::new();
            $(
                v.push(($v1.to_string(), $v2.to_string()));
            )*
            crate::graph::Graph::from(v)
        }
    };
}


mod graph_test {
    #[test]
    fn init_test() {
        // test 1
        let graph = graph!{
            ["v0", "v1"],
            ["v1", "v1"]
        };
        assert_eq!(graph.vertices_labeling, vec!["v0".to_string(), "v1".to_string()]);
        assert_eq!(graph.n_vertices(), 2);
        assert_eq!(graph.n_non_loop_edges(), 1);
        assert_eq!(graph.n_non_loop_edges_up_to_multiplicity(), 1);
        assert_eq!(graph.n_loops(), 1);
        assert_eq!(graph.data[0][1].len(), 1);
        assert_eq!(graph.data[1][1].len(), 1);


        // test 2
        let graph = graph!{
            ["v0", "v0"],
            ["v0", "v1"],
            ["v1", "v0"],
            ["v1", "v1"]
        };
        assert_eq!(graph.vertices_labeling, vec!["v0".to_string(), "v1".to_string()]);
        assert_eq!(graph.n_vertices(), 2);
        assert_eq!(graph.n_non_loop_edges(), 2);
        assert_eq!(graph.n_non_loop_edges_up_to_multiplicity(), 1);
        assert_eq!(graph.n_loops(), 2);
        assert_eq!(graph.data[0][0].len(), 1);
        assert_eq!(graph.data[0][1].len(), 2);
        assert_eq!(graph.data[1][1].len(), 1);


        // test 3
        let graph = graph! {
            ["v0", "v1"],
            ["v1", "v2"],
            ["v1", "v3"],
            ["v0", "v3"],
            ["v0", "v4"]
        };
        assert_eq!(graph.vertices_labeling, vec!["v0".to_string(), "v1".to_string(), "v2".to_string(), "v3".to_string(), "v4".to_string()]);
        assert_eq!(graph.data[0][1].len(), 1);
        assert_eq!(graph.data[1][2].len(), 1);
        assert_eq!(graph.data[1][3].len(), 1);
        assert_eq!(graph.data[0][3].len(), 1);
        assert_eq!(graph.data[0][4].len(), 1);
    }
}

#[allow(unused)]
#[derive(Clone)]
pub struct VertexVector {
    start: String,
    end: String,
    direction: usize,
    pos: Rational,
}

#[allow(unused)]
impl VertexVector {
    fn new(start: String, end: String, direction: usize, pos: Rational ) -> VertexVector {
        VertexVector { start: start, end: end, direction: direction, pos: pos }
    }
}

impl fmt::Display for VertexVector {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({})({}, {})", self.pos, self.start, self.end)?;
        if self.direction != 0 {
            write!(f, "{}", self.direction)?;
        };
        write!(f, "")
    }
}

#[allow(unused)]
#[derive(Clone)]
pub enum SimpleVertex {
    Essential(String),
    Redundant(VertexVector),
}

impl fmt::Display for SimpleVertex {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.clone() {
            SimpleVertex::Essential(label) => write!(f, "{}", label),
            SimpleVertex::Redundant(vec) => write!(f, "{}", vec),
        }
    }
}


#[derive(Eq, PartialEq, Debug, Clone)]
pub struct RawSimpleGraph {
    data: BinaryPool,
    n_vertices: usize,
    pub is_maximal_tree_chosen: bool,
}

impl RawSimpleGraph {
    pub fn new(n_vertices: usize) -> RawSimpleGraph {
        RawSimpleGraph {
            data: BinaryPool::new(n_vertices * (n_vertices-1)/2),
            n_vertices: n_vertices,
            is_maximal_tree_chosen: false,
        }
    }

    // function returing the raw index for self.data. The condition that i < j is required.
    fn raw_idx(i: usize, j: usize) -> usize {
        debug_assert!(i < j, "The condition i < j is required, but not satisfied: i={} but j={}", i, j);
        j*(j-1)/2 + i
    }

    pub fn add_edge(&mut self, i: usize, j: usize) {
        self.data.add(RawSimpleGraph::raw_idx(i, j));
    }

    pub fn remove_edge(&mut self, i: usize, j: usize) {
        self.data.remove(RawSimpleGraph::raw_idx(i, j));
    }

    pub fn contains(&self, i: usize, j: usize) -> bool {
        self.data.contains(RawSimpleGraph::raw_idx(i, j))
    }

    pub fn n_edges(&self) -> usize {
        self.data.n_elements()
    }

    pub fn n_vertices(&self) -> usize {
        self.n_vertices
    }

    pub fn swap(&mut self, i: usize, j: usize) {
        debug_assert!(i < j && j < self.n_vertices(), "Two condition i <= j and j < self.n_vertices() are required, but not satisfied: i={}, j={}, self.n_vertices()={}.", i, j, self.n_vertices());

        if i ==j { return; }

        let isotopy = self.contains(i,j);
        // Now i < j

        let tmp_i: Vec<_> = self.adjacent_vertices_iter(i).filter(|&v| v != j ).map(|v| if j < v {(j, v)} else {(v, j)} ).collect();
        let tmp_j: Vec<_> = self.adjacent_vertices_iter(j).filter(|&v| v != i ).map(|v| if i < v {(i, v)} else {(v, i)} ).collect();
        self.remove_edges_adjacent_to(i);
        self.remove_edges_adjacent_to(j);
        tmp_i.iter().for_each(|&(i,j)| self.add_edge(i,j));
        tmp_j.iter().for_each(|&(i,j)| self.add_edge(i,j));

        if isotopy==true { self.add_edge(i,j); };
    }

    fn remove_edges_adjacent_to(&mut self, i: usize) {
        let edges_to_delete: Vec<_> = self.adjacent_vertices_iter(i).map(|j| if i<j {(i,j)} else {(j,i)} ).collect();
        edges_to_delete.iter().for_each(|&(i, j)| self.remove_edge(i, j) );
    }

    fn degree_of(&self, idx: usize) -> usize {
        if idx >= self.n_vertices() { panic!(); };

        (0..self.n_vertices())
        .filter(|&j| j != idx )
        .map(|j| {
            if j < idx {
                if self.contains(j, idx) { 1 } else { 0 }
            } else {
                if self.contains(idx, j) { 1 } else { 0 }
            }
        })
        .sum()
    }

    pub fn is_connected(&self) -> bool {
        // we store all the vertex visited, and see if this array gets full (which means that the graph is connected)
        let mut vertices_visited = BinaryPool::new(self.n_vertices());
        self.adjacent_vertices_iter(0).for_each( |v| vertices_visited.add(v) );
        let mut new_vertices: Vec<_> = self.adjacent_vertices_iter(0).collect();

        loop {
            // base case
            if vertices_visited.is_full() { return true }

            let v = match new_vertices.last() {
                Some(&v) => v,
                None => return false,
            };
            new_vertices.pop();
            self.adjacent_vertices_iter(v).for_each( |w| {
                if !vertices_visited.contains(w) {
                    vertices_visited.add(w);
                    new_vertices.push(w);
                }
            });
        }
    }
}


pub struct AdjacentVerticesIter<'a> {
    raw_graph: &'a RawSimpleGraph,
    curr: usize,
    vertex: usize,
}

impl<'a> AdjacentVerticesIter<'a> {
    fn new(raw_graph: &'a RawSimpleGraph, vertex: usize) -> AdjacentVerticesIter<'a> {
        AdjacentVerticesIter { raw_graph: raw_graph, curr: 0, vertex: vertex }
    }
}

impl<'a> Iterator for AdjacentVerticesIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr >= self.raw_graph.n_vertices() { return None; }

        let mut out = None;

        if self.vertex < self.curr {
            if self.raw_graph.data.contains(RawSimpleGraph::raw_idx(self.vertex, self.curr)) {
                out = Some(self.curr);
            }
        } else if self.vertex > self.curr {
            if self.raw_graph.data.contains(RawSimpleGraph::raw_idx(self.curr, self.vertex)) {
                out = Some(self.curr);
            }
        }

        self.curr += 1;
        if out == None {
           self.next()
        } else {
            out
        }
    }
}


impl<'a> RawSimpleGraph {
    pub fn adjacent_vertices_iter(&'a self, idx: usize) -> AdjacentVerticesIter<'a> {
        AdjacentVerticesIter::new(self, idx)
    }
}


// Display
impl std::fmt::Display for RawSimpleGraph {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for j in 1..self.n_vertices() {
            for i in 0..j {
                write!(f, "{}", if self.contains(i, j) {1} else {0} )?;
            }
            write!(f, "\n")?;
        }
        write!(f, "")
    }
}


// symmetric action on 'RawSimpleGraph'
use algebra::non_commutative::Gset;
impl Gset<Permutation> for RawSimpleGraph {
    fn gen_action_by(&mut self, elem: Permutation) {
        // First, process the case that 'elem' is the identity
        if elem == Permutation::identity() { return; }

        // again, self is a generator, so it must be a 2-cycle in our convension.
        let unfixed_indeces: Vec<_> = elem.unfixed_indeces_iter().collect();

        debug_assert!(unfixed_indeces.len() == 2, "gen_act must receive a 2-cycle, but it receives: {}", unfixed_indeces.len());

        let (i, j) = (unfixed_indeces[0], unfixed_indeces[1]);
        self.swap(i, j);
    }
}

impl std::ops::Shl<Permutation> for RawSimpleGraph {
    type Output = RawSimpleGraph;

    fn shl(mut self, elem: Permutation) -> Self::Output {
        self.action_by(elem);
        self
    }
}

impl std::ops::ShlAssign<Permutation> for RawSimpleGraph {
    fn shl_assign(&mut self, elem: Permutation) {
        self.action_by(elem);
    }
}


// debugging purpose macro
#[allow(unused)]
macro_rules! raw_graph {
    (n_vertices=$n_vertices: expr, $( ($v: expr, $w: expr) ), + ) => {
        {
            let mut raw_graph = RawSimpleGraph::new($n_vertices);
            $(
                raw_graph.add_edge($v, $w);
            )*
            raw_graph
        }
    };
}



#[cfg(test)]
mod raw_simple_graph_test {
    use crate::graph::RawSimpleGraph;
    use algebra::{non_commutative::{permutation::*, Group}, permutation};

    #[test]
    fn iter_and_erase_test() {
        let mut raw_graph = raw_graph!{ n_vertices=10, 
            (0, 5),
            (1, 5),
            (0, 6),
            (3, 6),
            (0, 7),
            (4, 7),
            (1, 8),
            (2, 8),
            (1, 9),
            (3, 9)
        };
        assert_eq!(raw_graph.adjacent_vertices_iter(0).collect::<Vec<_>>(), vec![5,6,7]);
        raw_graph.remove_edges_adjacent_to(0);
        assert_eq!(raw_graph.adjacent_vertices_iter(0).collect::<Vec<_>>(), vec![]);
    }

    #[test]
    fn degree_of_test () {
        let raw_graph = raw_graph!{ n_vertices=10, 
            (0, 5),
            (1, 5),
            (0, 6),
            (3, 6),
            (0, 7),
            (4, 7),
            (1, 8),
            (2, 8),
            (1, 9),
            (3, 9)
        };
        assert_eq!(raw_graph.degree_of(0), 3);
        assert_eq!(raw_graph.degree_of(1), 3);
        assert_eq!(raw_graph.degree_of(2), 1);
        assert_eq!(raw_graph.degree_of(3), 2);
        assert_eq!(raw_graph.degree_of(4), 1);
        assert_eq!(raw_graph.degree_of(5), 2);
        assert_eq!(raw_graph.degree_of(6), 2);
        assert_eq!(raw_graph.degree_of(7), 2);
        assert_eq!(raw_graph.degree_of(8), 2);
        assert_eq!(raw_graph.degree_of(9), 2);
    }


    #[test]
    fn swap_test() {
        let mut raw_graph = raw_graph!{ n_vertices=10, 
            (0, 5),
            (1, 5),
            (0, 6),
            (3, 6),
            (0, 7),
            (4, 7),
            (1, 8),
            (2, 8),
            (1, 9),
            (3, 9)
        };

        raw_graph.swap(1, 5);
        let answer = raw_graph!{ n_vertices=10, 
            (0, 1),
            (0, 6),
            (0, 7),
            (1, 5),
            (2, 8),
            (3, 6),
            (3, 9),
            (4, 7),
            (5, 8),
            (5, 9)
        };
        assert_eq!(raw_graph, answer);

        raw_graph.swap(2, 5);
        raw_graph.swap(3, 8);
        raw_graph.swap(4, 5);
        raw_graph.swap(5, 9);
        raw_graph.swap(6, 8);
        raw_graph.swap(7, 8);
        let answer = raw_graph!{ n_vertices=10, 
            (0, 1),
            (0, 7),
            (0, 8),
            (1, 2),
            (2, 3),
            (2, 5),
            (3, 4),
            (5, 6),
            (6, 7),
            (8, 9)
        };
        assert_eq!(raw_graph, answer);
    }


    #[test]
    fn connectedness_test() {
        // test 1
        let raw_graph = raw_graph!{ n_vertices=10, 
            (0, 5),
            (1, 5),
            (0, 6),
            (3, 6),
            (0, 7),
            (4, 7),
            (1, 8),
            (2, 8),
            (1, 9),
            (3, 9)
        };
        assert!( raw_graph.is_connected() );

        // test 2
        let raw_graph = raw_graph!{ n_vertices=2, 
            (0, 1)
        };
        assert!( raw_graph.is_connected() );

        // test 3
        let raw_graph = raw_graph!{ n_vertices=4, 
            (0, 1),
            (2, 3)
        };
        assert!( !raw_graph.is_connected() );

        // test 4
        let raw_graph = RawSimpleGraph::new(10);
        assert!( !raw_graph.is_connected() );

        // test 5
        let raw_graph = RawSimpleGraph::new(1);
        assert!( !raw_graph.is_connected() );

        // test 6
        let raw_graph = raw_graph!{ n_vertices=10, 
            (0, 5),
            (1, 5),
            (2, 6),
            (3, 6),
            (4, 5),
            (5, 6),
            (7, 9),
            (8, 9)
        };
        assert!( !raw_graph.is_connected() );

        // test 7
        let raw_graph = raw_graph!{ n_vertices=10, 
            (2, 3),
            (3, 4),
            (2, 5),
            (5, 6),
            (6, 7),
            (0, 7),
            (0, 8),
            (8, 9)
        };
        assert!( !raw_graph.is_connected() );

    }
    
    #[test]
    fn symmetric_action_test() {
        let mut graph = raw_graph!{ n_vertices=10, 
            (0, 5),
            (1, 5),
            (0, 6),
            (3, 6),
            (0, 7),
            (4, 7),
            (1, 8),
            (2, 8),
            (1, 9),
            (3, 9)
        };
        let permutations = vec![ 
            permutation!{ (1,2,4,6) },
            permutation!{ (0,3,4,2,5) },
            permutation!{ (5,6) },
            permutation!{ (1,2,3,5,8,9) }
        ];

        // inverse test
        let graph_cpy = graph.clone();
        for p in permutations.iter() {
            graph <<= p.clone();
            graph <<= p.clone().inverse();
            assert_eq!(graph, graph_cpy);
        }

        // mutiplication test
        let mut graph1 = graph.clone();
        let mut graph2 = graph.clone();
        
        graph1 <<= permutations[0].clone();
        graph1 <<= permutations[1].clone();

        graph2 <<= permutations[1].clone() * permutations[0].clone();

        assert_eq!(graph1, graph2);

        // product test
        let mut graph1 = graph.clone();
        let mut graph2 = graph.clone();

        permutations.iter().for_each(|p| graph1 <<= p.clone());
        graph2 <<= permutations.iter().map(|p|p.clone()).product();

        assert_eq!(graph1, graph2);
    }
}


impl RawSimpleGraph {
    fn choose_maximal_tree(&mut self) -> Permutation {
        // The following permutation records the permutation of vertices done in this function.
        let mut p = Permutation::identity();

        // Step 1: choose a base point.
        p = self.choose_base_point() * p;

        // Step 2: choose a maximal tree in the depth first manner
        p = self.choose_maximal_tree_by_depth_first() * p;

        // Now the maximal tree is chosen:
        self.is_maximal_tree_chosen = true;

        // Step 3: Modify the ordering of the tree in such a way that any non-separating edge come first
        p = self.modify_ordering_of_the_branches() * p;

        p
    }

    fn choose_base_point(&mut self) -> Permutation {
        // -- if the vertex 0 is of degree one, we can just let this point to be the base point
        if self.degree_of(0) != 1 {
            // -- otherwuse, we want a vertex v such that the graph with a vertex v removed is connected.
            for i in 0..self.n_vertices() {
                // remove edges adjacent to i and see if the resulting graph is connected
                let tmp_edges: Vec<_> = self.adjacent_vertices_iter(i).map(|j| if i<j {(i, j)} else {(j, i)} ).collect();
                self.remove_edges_adjacent_to(i);
                self.add_edge( tmp_edges.first().unwrap().0, tmp_edges.first().unwrap().1 ); // add one edge back in, so that the vertex i belongs to one of the components
                let connected = self.is_connected();
                // recover the graph
                tmp_edges.iter().for_each(|&(i, j)| self.add_edge(i, j) );

                //  if connected bring the vertex to the front
                if connected { 
                    if i==0 {
                        return Permutation::identity()
                    } else {
                        self.swap(0, i);
                        // return the permutation that represents this swap
                        return permutation!{ (0, i) };
                    }
                }
            }
            panic!("good basepoint does not exist.");
        };
        Permutation::identity()
    }

    // used in 'choose_maximal_tree_by_depth_first(_)'. Find edge (j,k)
    fn fetch_right_and_swap( &mut self, i: usize, j: usize ) -> Option<Permutation> {
        let op_k = self.adjacent_vertices_iter(i).filter(|&k| k > j).next();
        match op_k {
            Some(k) => { 
                let p = permutation![(j,k)];
                *self <<= p.clone(); 
                Some( p )
            },
            None => None,
        }
    }

    fn choose_maximal_tree_by_depth_first(&mut self) -> Permutation {
        // The following 'Permutation' is the record of action done in this function, which will also be the return value of this function
        let mut p = Permutation::identity();

        for j in 1..self.n_vertices() {
            for i in (0..j).rev() {
                if self.contains(i, j) { break; }
                let result = self.fetch_right_and_swap(i, j);
                if let Some(q) = result { 
                    p = q * p;
                    break; 
                }  
            }
        }

        // return the record of action
        p
    }

    fn modify_ordering_of_the_branches(&mut self) -> Permutation {
        // The following 'Permutation' is the record of action done in this function, which will also be the return value of this function
        let mut p = Permutation::identity();
        
        for i in 0..self.n_vertices()-1 {
            if self.degree_of(i) <= 2 { continue; }

            // Now i is an essential vertex
            let mut start = i+1;
            let mut branches: Vec<_> = (i+2..self.n_vertices())
                .filter(|&j| self.contains(i,j) )
                .filter(|&j| self.maximal_tree_contains(i,j) )
                .map(|j| { 
                    let range = start..j;
                    let i_j_is_separating = (0..i).any(|k| (start..j).any(|l| self.contains(k, l)));
                    start = j; // update 'start'
                    (range, i_j_is_separating)
                })
                .collect();

            // add the last branch
            {
                let range = start..self.n_vertices();
                let i_j_is_separating = (0..i).any(|k| (start..self.n_vertices()).any(|l| self.contains(k, l)));
                branches.push((range, i_j_is_separating));
            }

            // we order the branshes in such a way that any branch with the second element 'false' comes earlier than that with 'true';
            let mut end = self.n_vertices();
            let action = branches.iter()
                .rev()
                .filter(|&&(_, i_j_is_separating)| i_j_is_separating)
                .map(|(branch, _)|{
                    let branch_size = branch.end - branch.start;
                    let out = if end == branch.end {
                        Permutation::identity()
                    } else {
                        Permutation::swap(branch.clone(), branch.end..end)
                    };
                    end -= branch_size; // update 'end'
                    out
                })
                .product::<Permutation>();

            *self <<= action.clone();
            p = action * p;
        }
        return p;
    }

    fn maximal_tree_contains(&self, i: usize, j: usize ) -> bool {
        if !self.is_maximal_tree_chosen { panic!("Trying to access the maximal tree, but the tree is not chosen yet.") };

        if !self.contains(i, j) { return false; };

        if (i+1..j).any(|k| self.contains(k, j)) {
            return false;
        };

        return true;
    }
}


#[cfg(test)]
mod maximal_graph_test {
    use super::RawSimpleGraph;
    use algebra::{non_commutative::{permutation::*, Group}, permutation};

    fn check_edges(before_choosing_maximal_tree: &RawSimpleGraph, after_choosing_maximal_tree: &RawSimpleGraph, p: &Permutation) {
        for j in 0..before_choosing_maximal_tree.n_vertices() {
            for i in 0..j {
                if !before_choosing_maximal_tree.contains(i,j) { continue; };
                let i_ = std::cmp::min( p.map(i), p.map(j) );
                let j_ = std::cmp::max( p.map(i), p.map(j) );
                if !after_choosing_maximal_tree.contains(i_, j_) {
                    panic!("'before' contains ({}, {}), but 'after' does not.", i, j);
                }
            };
        };
    }

    fn check_maximality_and_uniqueness_of_zero_cell (graph: &RawSimpleGraph) {
        let mut record_of_vertices = Vec::new();
        check_maximality_and_uniqueness_of_zero_cell_recc(graph, 0, &mut record_of_vertices);
        assert_eq!(record_of_vertices, (1..graph.n_vertices()).collect::<Vec<_>>());
    }

    fn check_maximality_and_uniqueness_of_zero_cell_recc (
        graph: &RawSimpleGraph, 
        i:usize, 
        record_of_vertices: &mut Vec<usize>
    ) {

        // base case
        if i > graph.n_vertices()-2 { return; };
        if !graph.contains(i,i+1) { return; } 

        // general case
        graph.adjacent_vertices_iter(i)
            .filter(|&j| j>i)
            .filter(|&j| graph.contains(i, j))
            .take_while(|&j| graph.maximal_tree_contains(i, j))
            .for_each(|j| {
                record_of_vertices.push(j);
                check_maximality_and_uniqueness_of_zero_cell_recc(graph, j, record_of_vertices);
            })
    }

    #[test]
    fn trivial_case() {
        let mut tree = raw_graph!{ n_vertices=10, 
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (2, 6),
            (6, 7),
            (2, 8),
            (8, 9)
        };
        let mut cpy = tree.clone();
        cpy.is_maximal_tree_chosen = true;

        let p = tree.choose_maximal_tree();
        assert_eq!(p, Permutation::identity());
        assert_eq!(tree, cpy);

        check_edges(&cpy, &tree, &p);
        check_maximality_and_uniqueness_of_zero_cell( &tree );
    }

    #[test]
    fn choose_base_point_test() {
        let mut graph = raw_graph!{ n_vertices=10, 
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (2, 5),
            (5, 6),
            (6, 7),
            (0, 7),
            (0, 8),
            (8, 9)
        };

        let p = graph.choose_base_point();
        assert_eq!(p, permutation!{ (0,1) });
    }

    #[test]
    fn choose_maximal_tree_by_depth_first_test() {
        let mut graph = raw_graph!{ n_vertices=10, 
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (2, 5),
            (5, 6),
            (6, 7),
            (0, 7),
            (0, 8),
            (8, 9)
        };

        let mut graph_cpy = graph.clone();

        let p = graph.choose_base_point();
        let p = graph.choose_maximal_tree_by_depth_first() * p;
        assert_eq!(p, permutation!{(0,1)(2,5,4,7)(3,6)} );

        graph_cpy <<= p.clone();
        assert_eq!(graph, graph_cpy);
    }

    #[test]
    fn modify_ordering_test() {
        let mut graph = raw_graph!{ n_vertices=10, 
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (2, 5),
            (5, 6),
            (6, 7),
            (0, 7),
            (0, 8),
            (8, 9)
        };

        let mut graph_cpy = graph.clone();
        let before = graph.clone();

        let p = graph.choose_base_point();
        let p = graph.choose_maximal_tree_by_depth_first() * p;
        graph.is_maximal_tree_chosen = true;
        let p = graph.modify_ordering_of_the_branches() * p;

        assert_eq!(p, permutation!{(0,1)(2, 7, 4, 9, 3, 8)(5,6)} );

        graph_cpy <<= p.clone();
        graph_cpy.is_maximal_tree_chosen = true;
        assert_eq!(graph, graph_cpy);

        // if we come here, the resulting 'graph' is endowed with the maiximal tree.
        // using this graph, we do a test for the test functions.
        check_edges(&before, &graph, &p);
        check_maximality_and_uniqueness_of_zero_cell(&graph);

    }
}



#[allow(unused)]
pub struct SimpleGraph {
    data: RawSimpleGraph,
    vertices_list: Vec<SimpleVertex>,
    n_essential_vertices: usize,
}


// non public interfaces of SimpleGraph
#[allow(unused)]
impl SimpleGraph {
    fn new() -> SimpleGraph {
        SimpleGraph {
            data: RawSimpleGraph::new(0),
            vertices_list: Vec::new(),
            n_essential_vertices: 0,
        }
    }

    pub fn from(graph: &Graph ) -> Self {
        Self::subdivided_graph_from(graph, 2)
    }

    // this function is needed by the test functions outside this module
    pub fn without_lables_from( graph: RawSimpleGraph ) -> Self {
        let n_vertices  = graph.n_vertices();
        SimpleGraph {
            data: graph,
            vertices_list: vec![SimpleVertex::Essential("".to_string()); n_vertices],
            n_essential_vertices: 0,
        }
    }

    pub fn subdivided_graph_from(graph: &Graph, n_marked_points: usize) -> Self {
        // When n_marked_points < 2, we create a graph subdivided for n_marked_points = 2 (This will guerantee that the resulting graph is simple).
        let n_marked_points = if n_marked_points < 2 { 
            2
        } else {
            n_marked_points
        };

        let n_vertices = 
            graph.n_vertices() 
            + graph.n_non_loop_edges() * std::cmp::max(n_marked_points-2, 1)
            + graph.n_loops() * n_marked_points;

        let data = RawSimpleGraph::new(n_vertices );

        // Adding vertices
        let vertices_list:Vec<_> = (0..graph.n_vertices()).map(|i| SimpleVertex::Essential(graph.vertices_labeling[i].clone())).collect();

        let mut subdivided_graph = SimpleGraph {
            data: data,
            vertices_list: vertices_list,
            n_essential_vertices: graph.n_vertices(),
        };

        // for each essential vertex...
        for i in 0..subdivided_graph.n_essential_vertices() {
            // Iterating over non-loop edges
            for j in i+1..subdivided_graph.n_essential_vertices() {
                if graph.data[i][j].len() == 0 {
                    continue;
                }

                let n_mid_vertices = std::cmp::max(n_marked_points-2, 1);


                for (k, edge) in graph.data[i][j].iter().enumerate() {
                    subdivided_graph.data.add_edge(i, subdivided_graph.n_vertices());
                    subdivided_graph.vertices_list.push(SimpleVertex::Redundant(VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[j].clone(), k, 1.over(n_marked_points-1))));

                    for l in 2..=n_mid_vertices {
                        subdivided_graph.data.add_edge(subdivided_graph.n_vertices()-1, subdivided_graph.n_vertices());
                        subdivided_graph.vertices_list.push(SimpleVertex::Redundant(VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[j].clone(), k, l.over(n_marked_points-1))));
                    }

                    subdivided_graph.data.add_edge(j, subdivided_graph.n_vertices()-1);
                }
            }

            // Finally, we process loops, i.e. edges of the form [i, i]
            let n_mid_vertices = n_marked_points;

            for (k, edge) in graph.data[i][i].iter().enumerate() {
                subdivided_graph.data.add_edge(i, subdivided_graph.n_vertices());
                subdivided_graph.vertices_list.push(SimpleVertex::Redundant(VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[i].clone(), k, 1.over(n_marked_points+1))));

                for l in 2..=n_mid_vertices {
                    subdivided_graph.data.add_edge(subdivided_graph.n_vertices()-1, subdivided_graph.n_vertices());
                    subdivided_graph.vertices_list.push(SimpleVertex::Redundant(VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[i].clone(), k, l.over(n_marked_points+1))));
                }

                subdivided_graph.data.add_edge(i, subdivided_graph.n_vertices()-1);
            }
        }

        let p = subdivided_graph.data.choose_maximal_tree();
        subdivided_graph.vertices_list <<= p;

        subdivided_graph
    }

    pub fn n_vertices(&self) -> usize{
        self.vertices_list.len()
    }

    pub fn n_essential_vertices(&self) -> usize{
        self.n_essential_vertices
    }
}

impl std::ops::Deref for SimpleGraph {
    type Target = RawSimpleGraph;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}


#[cfg(test)]
mod simple_graph_test {
    use crate::graph::RawSimpleGraph;
    use super::SimpleGraph;
    #[test]
    fn init_test() {
        let graph = graph!{
            ["v0", "v0"],
            ["v0", "v1"],
            ["v0", "v1"],
            ["v1", "v1"]
        };

        let simple_graph = SimpleGraph::from(&graph);
        assert_eq!(simple_graph.n_vertices(), 8);
    }

    #[test]
    fn subdivision_test() {
        // test 1
        let graph = graph!{
            ["v0", "v1"],
            ["v1", "v1"]
        };

        let subdivided_graph = SimpleGraph::subdivided_graph_from(&graph, 3);
        assert_eq!(subdivided_graph.n_vertices(), 6);
        let mut answer_data = RawSimpleGraph::new(6);
        {
            answer_data.add_edge(0, 1);
            answer_data.add_edge(1, 2);
            answer_data.add_edge(2, 3);
            answer_data.add_edge(3, 4);
            answer_data.add_edge(4, 5);
            answer_data.add_edge(2, 5);
            answer_data.is_maximal_tree_chosen = true;
        }
        assert_eq!(subdivided_graph.data, answer_data);
    }
}


impl fmt::Display for SimpleGraph {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "vertices: ")?;
        for v_idx in 0..self.vertices_list.len() {
            // printing each cell
            write!(f, "{}", self.vertices_list[v_idx])?;
            if v_idx != self.vertices_list.len()-1 {
                write!(f, ", ")?;
            };
        };

        write!(f, "\n")?;

        write!(f, "edges: ")?;
        for (i, (v, w)) in self.edge_iter().map(|e| e.edge()).enumerate() {
            // printing each cell
            write!(f, "[{}, {}]", self.vertices_list[v], self.vertices_list[w])?;
            if i != self.n_edges()-1 {
                write!(f, ", ")?;
            };
        };

        write!(f, "\n")
    }
}


use std::ops::Index;
impl Index<usize> for SimpleGraph {
    type Output = SimpleVertex;

    fn index(&self, idx: usize) -> &Self::Output {
        if idx >= self.vertices_list.len() { panic!(); }
        &self.vertices_list[idx]
    }
}


#[derive(Eq, PartialEq, Clone, Copy, Ord, Debug)]
pub enum CubeFactor {
    Edge(usize, usize),
    Vertex(usize)
}

impl CubeFactor {
    pub fn is_disjoint_from(self, other: Self) -> bool {
        match (self, other) {
            (Self::Edge(w,x), Self::Edge(y,z)) => w != y && w != z && x != y && x != z,
            (Self::Vertex(x), Self::Edge(y,z)) =>  x != y && x != z,
            (Self::Edge(x, y), Self::Vertex(z)) =>  x != z && y != z,
            (Self::Vertex(x), Self::Vertex(y)) =>  x != y
        }
    }

    pub fn is_in_maximal_tree(self, graph: &RawSimpleGraph) -> bool {
        match self {
            Self::Edge(x, y) => graph.maximal_tree_contains(x, y),
            Self::Vertex(_) => true,
        }
    }

    pub fn is_edge(self) -> bool {
        match self {
            Self::Edge(_, _) => true,
            Self::Vertex(_) => false,
        }
    }

    pub fn is_vertex(self) -> bool {
        !self.is_edge()
    }

    pub fn edge(self) -> (usize, usize) {
        match self {
            Self::Edge(x, y) => (x, y),
            Self::Vertex(x) => panic!("unwrap Failed: ({}) is not an edge", x),
        }
    }

    pub fn vertex(self) -> usize {
        match self {
            Self::Edge(x, y) => panic!("unwrap Failed: ({}, {}) is not a vertex", x, y),
            Self::Vertex(x) => x,
        }
    }

    pub fn is_order_respecting(self) -> bool {
        // note that 'self' must be an edge
        let (i, j) = self.edge();
        return j - i == 1
    }

    pub fn is_disrespected_by(self, other: CubeFactor) -> bool {
        // 'self' must be a vertex, and 'other' has to be an edge
        let (i, j) = other.edge();
        let k = self.vertex();

        if k<=i || k>=j { return false; }

        k - i == 1
    }

    pub fn terminal(self) -> CubeFactor {
        // note that 'self' must be an edge
        CubeFactor::Vertex( self.edge().0 )
    }

    pub fn initial(self) -> CubeFactor {
        // note that 'self' must be an edge
        CubeFactor::Vertex( self.edge().1 )
    }
}


pub struct VertexIter<'a> {
    _graph: &'a RawSimpleGraph,
    iter: std::ops::Range<usize>,
}

impl<'a> Iterator for VertexIter<'a> {
    type Item = CubeFactor;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(v) = self.iter.next() {
            Some( CubeFactor::Vertex(v) )
        } else {
            None
        }
    }
}

pub struct EdgeIter<'a> {
    graph: &'a RawSimpleGraph,
    curr: (usize, usize),
    is_done: bool,
}

impl<'a> Iterator for EdgeIter<'a> {
    type Item = CubeFactor;

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_done {
            return None;
        };

        let out = CubeFactor::Edge(self.curr.0, self.curr.1);
        
        loop {
            if self.curr.0 == 0 {
                self.curr.0 = self.curr.1;
                self.curr.1 += 1;
            } else {
                self.curr.0 -= 1;
            }

            if self.curr.1 >= self.graph.n_vertices() {
                self.is_done = true;
                break;
            } else if self.graph.contains(self.curr.0, self.curr.1) {
                break;
            };
        }

        Some(out)
    }
}

#[allow(unused)]
impl<'a> RawSimpleGraph {
    pub fn edge_iter(&'a self) -> EdgeIter<'a> {
        let edge = (0, 1);
        EdgeIter {
            graph: self,
            curr: edge,
            is_done: self.n_vertices() < 2,
        }
    }

    pub fn vertex_iter(&'a self) -> VertexIter<'a> {
        VertexIter {
            iter: (0..self.n_vertices()),
            _graph: self,
        }
    }
}


pub struct FactorIter<'a> {
    v_iter: VertexIter<'a>,
    e_iter: EdgeIter<'a>,
}

impl<'a> Iterator for FactorIter<'a> {
    type Item = CubeFactor;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(v) = self.v_iter.next() {
            Some(v)
        } else {
            self.e_iter.next()
        }
    }
}

impl<'a> RawSimpleGraph {
    pub fn iter(&'a self) -> FactorIter<'a> {
        FactorIter { 
            v_iter: self.vertex_iter(), 
            e_iter: self.edge_iter(),
        }
    }
}

use std::cmp::Ordering;
impl PartialOrd for CubeFactor {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // they are equal if they are equal
        if self == other { return Some( Ordering::Equal ) };

        // if the two cells are not dijoint, then they are incomparable.
        if !self.is_disjoint_from(*other) {
            return None;
        }

        // otherwise the factors are comparable
        match (self, other) {
            (Self::Edge(x, _), Self::Edge(y, _)) => Some(x.cmp(y)),
            (Self::Vertex(x), Self::Edge(y, _)) =>  Some(x.cmp(y)),
            (Self::Edge(x, _), Self::Vertex(y)) => Some(x.cmp(y)),
            (Self::Vertex(x), Self::Vertex(y)) =>  Some(x.cmp(y)),
        }
    } 
}


impl fmt::Display for CubeFactor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CubeFactor::Edge(x, y) => write!(f, "[{} {}]", x, y),
            CubeFactor::Vertex(x) => write!(f, "{}", x),
        }
    }
}


#[cfg(test)]
mod simple_graph_utility_test {
    use super::SimpleGraph;
    use super::CubeFactor;

    #[test]
    fn display_test() {
        let graph = graph!{
            ["v0", "v1"],
            ["v1", "v1"]
        };

        let subdivided_graph = SimpleGraph::subdivided_graph_from(&graph, 3);
        let description = subdivided_graph.to_string();
        let ans = 
            "vertices: v0, (1/2)(v0, v1), v1, (1/4)(v1, v1), (1/2)(v1, v1), (3/4)(v1, v1)\nedges: [v0, (1/2)(v0, v1)], [(1/2)(v0, v1), v1], [v1, (1/4)(v1, v1)], [(1/4)(v1, v1), (1/2)(v1, v1)], [(1/2)(v1, v1), (3/4)(v1, v1)], [v1, (3/4)(v1, v1)]\n"
            .to_string();
        assert_eq!(description, ans);
    }

    #[test]
    fn edge_iter_test() {
        // test 1
        let graph = graph!{
            ["v0", "v1"],
            ["v1", "v1"]
        };

        let subdivided_graph = SimpleGraph::subdivided_graph_from(&graph, 2);
        let edges:Vec<_> = subdivided_graph.edge_iter().collect();
        let answer = vec![
            CubeFactor::Edge(0,1),
            CubeFactor::Edge(1,2),
            CubeFactor::Edge(2,3),
            CubeFactor::Edge(3,4),
            CubeFactor::Edge(2,4)
        ];
        assert_eq!(edges, answer);

        let subdivided_graph = SimpleGraph::subdivided_graph_from(&graph, 3);
        let edges:Vec<_> = subdivided_graph.edge_iter().collect();
        let answer = vec![
            CubeFactor::Edge(0,1),
            CubeFactor::Edge(1,2),
            CubeFactor::Edge(2,3),
            CubeFactor::Edge(3,4),
            CubeFactor::Edge(4,5),
            CubeFactor::Edge(2,5)
        ];
        assert_eq!(edges, answer);
    }
}
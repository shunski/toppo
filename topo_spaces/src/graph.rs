use alg::permutation;
use util::BinaryPool;
use alg::commutative::{PID, Rational};
use alg::lin_alg::{Matrix, ConstVector};
use alg::rational;
use std::collections::HashMap;
use std::fmt;
use alg::non_commutative::{permutation::*, Group};
use crate::cubical::Cube;
use std::collections::VecDeque;
use crate::cubical::EdgePath;

#[allow(unused)]
#[derive(Clone)]
struct Edge<T: PID> {
    weight: T,
}

impl<T: PID> Edge<T> {
    fn new() -> Edge<T> {
        Edge{
            weight: T::one(),
        }
    }
}


#[allow(unused)]
pub struct Graph<T: PID> {
    data: Vec<Vec<Vec<Edge<T>>>>,
    vertices_labeling: Vec<String>,

    n_non_loop_edges: usize,
    n_non_loop_edges_up_to_multiplicity: usize,
    n_loops: usize,
}


#[allow(unused)]
impl<T: PID> Graph<T> {
    pub fn from( edges: Vec<(String, String)> ) -> Self {
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
            graph::Graph::<i64>::from(v)
        }
    };

    ( $([$v1: expr, $v2: expr; weight=$weight: expr]), *) => {
        {
            let mut v = Vec::new();
            $(
                v.push(($v1.to_string(), $v2.to_string(), $weight));
            )*
            crate::graph::WeightedSimpleGraph::from(v)
        }
    };

    ( embedding=$embedding: expr, weight_fn=$weight_fn: expr, $([$v1: expr, $v2: expr]), *) => {
        {
            let mut v = Vec::new();
            $(
                v.push(($v1.to_string(), $v2.to_string(), $weight_fn( *$embedding.get($v1).unwrap(), *$embedding.get($v2).unwrap())));
            )*
            crate::graph::WeightedSimpleGraph::from(v)
        }
    };
}


#[cfg(test)]
mod graph_test {
    use crate::graph;
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
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct VertexVector {
    pub start: String,
    pub end: String,
    pub direction: usize,
    pub pos: Rational,
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
#[derive(Clone, PartialEq, Eq, Hash)]
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


#[derive(Eq, PartialEq, Clone)]
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

    // function returning the raw index for self.data. The condition that i < j is required.
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

    pub fn degree_of(&self, idx: usize) -> usize {
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

    pub fn get_adjacent_vertices(&self, v: usize) -> Vec<usize> {
        self.adjacent_vertices_iter(v).collect()
    }

    fn get_edge_path_recc<const FORWARD: bool>(&self, end_cube: &mut Cube, path: &mut VecDeque<[usize; 2]>) {
        let smallest_unblocked = end_cube
            .iter()
            .enumerate()
            .filter(|&(_, &x)| !x.is_blocked(&end_cube, &self))
            .min_by(|(_, x), (_,y)| x.cmp(y) );

        let (idx, &smallest_unblocked) = if let Some(x) = smallest_unblocked {
            x
        } else {
            return;
        };

        let new_vertex = smallest_unblocked.flow(&self).terminal();
        end_cube[idx] = new_vertex;

        if FORWARD {
            let edge = [smallest_unblocked.vertex(), new_vertex.vertex()];
            path.push_back(edge);
        } else {
            let edge = [new_vertex.vertex(), smallest_unblocked.vertex()];
            path.push_front(edge);
        }

        // run the algorithm recursively
        self.get_edge_path_recc::<FORWARD>(end_cube, path);
    }

    pub fn get_edge_path(&self, c: Cube) -> EdgePath {
        // let edge = c.0.iter().find(|&&f| f.is_edge() ).unwrap().edge();
        // if self.factor_pool.maximal_tree_contains(edge) {
        //     self.get_edge_path_of_edge_in_maximal_tree(c) 
        // } else {
        //     self.get_edge_path_of_edge_not_in_maximal_tree(c)
        // }
        
        
        let mut start = Cube::from(c.iter().map(|f| match f {
            CubeFactor::Vertex(_) => *f,
            CubeFactor::Edge([_, x]) => CubeFactor::Vertex(*x),
        }).collect::<Vec<_>>());

        let mut end = Cube::from(c.iter().map(|f| match f {
            CubeFactor::Vertex(_) => *f,
            CubeFactor::Edge([x, _]) => CubeFactor::Vertex(*x),
        }).collect::<Vec<_>>());

        let critical_edge = c.iter().find(|f| f.is_edge()).unwrap().edge();
        let critical_edge = [critical_edge[1], critical_edge[0]];

        let mut path = VecDeque::from([critical_edge]);
        self.get_edge_path_recc::<false>(&mut start, &mut path);
        self.get_edge_path_recc::<true>(&mut end, &mut path);

        let path = EdgePath::from( 
            path.into_iter().map(|e| vec![e]).collect::<Vec<_>>() 
        ).reduce_to_geodesic();

        path
    }

    pub fn get_edge_path_to_from_base_pt(&self, mut c: Cube, to_not_from: bool) -> EdgePath {
        assert!(c.iter().all(|f| f.is_vertex() ));

        let mut path = VecDeque::new();
        self.get_edge_path_recc::<true>(&mut c, &mut path);
        
        let path = EdgePath::from(
            path.into_iter().map(|e| vec![e]).collect::<Vec<_>>() 
        ).reduce_to_geodesic();
        
        if to_not_from {
            path 
        } else {
            path.inverse()
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

// Debug
impl std::fmt::Debug for RawSimpleGraph {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.n_vertices() {
            for _ in 0..i {
                write!(f, "  ")?;
            }
            if i<10 { write!(f, " ")?; }
            write!(f, "{}", i)?;
            for j in i+1..self.n_vertices() {
                write!(f, " {}", if self.contains(i, j) {"+"} else {" "} )?;
            }
            write!(f, "\n")?;
        }
        write!(f, "")
    }
}


// symmetric action on 'RawSimpleGraph'
use alg::non_commutative::Gset;
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
    use alg::{non_commutative::{permutation::*, Group}, permutation};

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
        // -- if the vertex 0 is of degree one, we can just let it be the base point
        if self.degree_of(0) != 1 {
            // -- otherwise, we want an essential vertex v such that the graph with v deleted is connected.
            for i in 0..self.n_vertices() {
                // if i is not essential, then continue.
                if self.degree_of(i) ==1 {
                    self.swap(0, i);
                    // return the permutation that represents this swap
                    return permutation!{ (0, i) };
                } if self.degree_of(i) >= 3 {
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
                .filter(|&j| self.maximal_tree_contains([i,j]) )
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

    pub fn maximal_tree_contains(&self, [i, j]: [usize; 2]) -> bool {
        if !self.is_maximal_tree_chosen { panic!("Trying to access the maximal tree, but the tree is not chosen yet.") };

        if !self.contains(i, j) { return false; };

        if (i+1..j).any(|k| self.contains(k, j)) {
            return false;
        };

        return true;
    }
}


#[cfg(test)]
mod maximal_tree_test {
    use super::RawSimpleGraph;
    use alg::{non_commutative::{permutation::*, Group}, permutation};

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
            .take_while(|&j| graph.maximal_tree_contains([i, j]))
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
        assert_eq!(p, permutation!{ (0,4) });
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
        assert_eq!(p, permutation!{(0,4)(1,3)(5,7)} );

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

        assert_eq!(p, permutation!{(0,4)(1,3)(5,9,6,8)} );

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

    pub fn from(graph: &Graph<i64> ) -> Self {
        Self::subdivided_graph_from(graph, 2)
    }

    pub fn subdivided_graph_from(graph: &Graph<i64>, n_marked_points: usize) -> Self {
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
                    subdivided_graph.vertices_list.push(SimpleVertex::Redundant(VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[j].clone(), k, rational!( 1; n_marked_points-1))));

                    for l in 2..=n_mid_vertices {
                        subdivided_graph.data.add_edge(subdivided_graph.n_vertices()-1, subdivided_graph.n_vertices());
                        subdivided_graph.vertices_list.push(SimpleVertex::Redundant(VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[j].clone(), k, rational!( l; n_marked_points-1))));
                    }

                    subdivided_graph.data.add_edge(j, subdivided_graph.n_vertices()-1);
                }
            }

            // Finally, we process loops, i.e. edges of the form [i, i]
            let n_mid_vertices = n_marked_points;

            for (k, edge) in graph.data[i][i].iter().enumerate() {
                subdivided_graph.data.add_edge(i, subdivided_graph.n_vertices());
                subdivided_graph.vertices_list.push(SimpleVertex::Redundant(VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[i].clone(), k, rational!( 1; n_marked_points+1))));

                for l in 2..=n_mid_vertices {
                    subdivided_graph.data.add_edge(subdivided_graph.n_vertices()-1, subdivided_graph.n_vertices());
                    subdivided_graph.vertices_list.push(SimpleVertex::Redundant(VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[i].clone(), k, rational!( l; n_marked_points+1))));
                }

                subdivided_graph.data.add_edge(i, subdivided_graph.n_vertices()-1);
            }
        }

        let p = subdivided_graph.data.choose_maximal_tree();
        subdivided_graph.vertices_list <<= p;

        subdivided_graph
    }

    pub fn n_vertices(&self) -> usize {
        self.vertices_list.len()
    }

    pub fn n_essential_vertices(&self) -> usize {
        self.n_essential_vertices
    }

    pub fn modify_embedding<const N: usize>(&self, embedding: HashMap<String, ConstVector<f64, N>>) -> HashMap<SimpleVertex, ConstVector<f64,N>> {
        assert_eq!( self.n_essential_vertices(), embedding.len() );

        let mut out = HashMap::new();

        for v in self.vertices_list.iter() {
            match v {
                SimpleVertex::Essential(s) => out.insert(v.clone(), *embedding.get(s).unwrap()),
                SimpleVertex::Redundant(s) => {
                    let start = *embedding.get(&s.start).unwrap();
                    let end = *embedding.get(&s.end).unwrap();
                    let ratio: f64 = s.pos.into();
                    let pos = start + (end - start) * ratio;
                    out.insert(v.clone(), pos)
                },
            };
        }

        out
    }

    pub fn vertices_list(&self) -> &[SimpleVertex] {
        &self.vertices_list
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
    use crate::graph;
    #[test]
    fn init_test() {
        let graph = graph!{
            ["v0", "v0"],
            ["v0", "v1"],
            ["v0", "v1"],
            ["v1", "v1"],
            ["v1", "v2"]
        };

        let simple_graph = SimpleGraph::from(&graph);
        assert_eq!(simple_graph.n_vertices(), 10);
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
        for (i, [v, w]) in self.edge_iter().map(|e| e.edge()).enumerate() {
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
    Edge([usize; 2]),
    Vertex(usize)
}

impl CubeFactor {
    pub fn is_disjoint_from(self, other: Self) -> bool {
        match (self, other) {
            (Self::Edge([w,x]), Self::Edge([y,z])) => w != y && w != z && x != y && x != z,
            (Self::Vertex(x), Self::Edge([y,z])) =>  x != y && x != z,
            (Self::Edge([x, y]), Self::Vertex(z)) =>  x != z && y != z,
            (Self::Vertex(x), Self::Vertex(y)) =>  x != y
        }
    }

    pub fn is_in_maximal_tree(self, graph: &RawSimpleGraph) -> bool {
        match self {
            Self::Edge([x, y]) => graph.maximal_tree_contains([x, y]),
            Self::Vertex(_) => true,
        }
    }

    pub fn is_edge(self) -> bool {
        match self {
            Self::Edge([_, _]) => true,
            Self::Vertex(_) => false,
        }
    }

    pub fn is_vertex(self) -> bool {
        !self.is_edge()
    }

    pub fn edge(self) -> [usize; 2] {
        match self {
            Self::Edge([x, y]) => [x, y],
            Self::Vertex(x) => panic!("unwrap failed: ({}) is not an edge", x),
        }
    }

    pub fn vertex(self) -> usize {
        match self {
            Self::Edge([x, y])=> panic!("unwrap failed: ({}, {}) is not a vertex", x, y),
            Self::Vertex(x) => x,
        }
    }

    pub fn is_order_respecting(self) -> bool {
        // note that 'self' must be an edge
        let [i, j] = self.edge();
        return j - i == 1
    }

    pub fn is_disrespected_by(self, other: CubeFactor) -> bool {
        // 'self' must be a vertex, and 'other' has to be an edge
        let [i, j] = other.edge();
        let k = self.vertex();

        i<k && k<j
    }

    pub fn is_blocked(self, factors: &[CubeFactor], graph: &RawSimpleGraph) -> bool {
        if let CubeFactor::Vertex(v) = self {
            // if this is the basepoint of the tree, then this is trivially blocked.
            if v == 0 { return true; }

            // otherwise, let it flow and check
            let w = self.flow(graph).terminal();
            if factors.iter().any(|&f| (f.is_vertex() && f == w) || (f.is_edge() && ( f.terminal() == w || f.initial() == w )) ) {
                true
            } else {
                false
            }
        } else {
            panic!("failed to unwrap: Input must be a vertex, but it is an edge");
        }
    }

    pub fn terminal(self) -> CubeFactor {
        // note that 'self' must be an edge
        CubeFactor::Vertex( self.edge()[0] )
    }

    pub fn initial(self) -> CubeFactor {
        // note that 'self' must be an edge
        CubeFactor::Vertex( self.edge()[1] )
    }

    pub fn flow(self, graph: &RawSimpleGraph) -> CubeFactor {
        if let CubeFactor::Vertex(v) = self {
            if v==0 { panic!("Cannot flow from the basepoint of the maximal tree."); };
            let terminal = (0..v).rev().find(|&i| graph.contains(i, v) && graph.maximal_tree_contains([i, v]) ).unwrap();
            CubeFactor::Edge( [terminal, v] )
        } else {
            panic!("Cannot to flow. Input must be a vertex.");
        }
    }


    // This function returns the numbring of the factor in a graph.
    // For Vertex(x), it is x.
    // For Edge([x, _]), it is x,.
    #[inline]
    pub fn get_order(&self) -> usize {
        match self {
            CubeFactor::Vertex(x) => *x,
            CubeFactor::Edge([x,_]) => *x,
        }
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

        let out = CubeFactor::Edge([self.curr.0, self.curr.1]);
        
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
            (Self::Edge([x, _]), Self::Edge([y, _])) => Some(x.cmp(y)),
            (Self::Vertex(x), Self::Edge([y, _])) =>  Some(x.cmp(y)),
            (Self::Edge([x, _]), Self::Vertex(y)) => Some(x.cmp(y)),
            (Self::Vertex(x), Self::Vertex(y)) =>  Some(x.cmp(y)),
        }
    } 
}


impl fmt::Display for CubeFactor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CubeFactor::Edge([x, y]) => write!(f, "[{} {}]", x, y),
            CubeFactor::Vertex(x) => write!(f, "{}", x),
        }
    }
}


#[cfg(test)]
mod simple_graph_utility_test {
    use super::SimpleGraph;
    use super::CubeFactor;
    use crate::graph;

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
            CubeFactor::Edge([0,1]),
            CubeFactor::Edge([1,2]),
            CubeFactor::Edge([2,3]),
            CubeFactor::Edge([3,4]),
            CubeFactor::Edge([2,4])
        ];
        assert_eq!(edges, answer);

        let subdivided_graph = SimpleGraph::subdivided_graph_from(&graph, 3);
        let edges:Vec<_> = subdivided_graph.edge_iter().collect();
        let answer = vec![
            CubeFactor::Edge([0,1]),
            CubeFactor::Edge([1,2]),
            CubeFactor::Edge([2,3]),
            CubeFactor::Edge([3,4]),
            CubeFactor::Edge([4,5]),
            CubeFactor::Edge([2,5])
        ];
        assert_eq!(edges, answer);
    }
}



// The functions in this module provides the graphs with maximal trees for the tests in other files
pub mod graphs_for_tests {
    // we manually create the graph of K_{5,5} with a maximal tree.
    use crate::graph::*;

    #[allow(unused)]
    pub fn k_5_5_with_4_marked_points() -> SimpleGraph {
        let mut raw_graph = RawSimpleGraph::new(25);
        // add deleted edges
        raw_graph.add_edge(0, 8);
        raw_graph.add_edge(0, 15);
        raw_graph.add_edge(0, 24);
        raw_graph.add_edge(3, 13);
        raw_graph.add_edge(3, 22);
        raw_graph.add_edge(6, 20);
        // add non-order-respecting edges in the maximal tree
        raw_graph.add_edge(6, 9);
        raw_graph.add_edge(11, 14);
        raw_graph.add_edge(11, 16);
        raw_graph.add_edge(18, 21);
        raw_graph.add_edge(18, 23);
        (0..25-1).filter(|&i| i!=8 && i!=13 && i!=15 && i!=20 && i!=22 ).for_each(|i| raw_graph.add_edge(i, i+1));

        raw_graph.is_maximal_tree_chosen = true;

        SimpleGraph {
            data: raw_graph,
            vertices_list: vec![SimpleVertex::Essential("".to_string()); 25],
            n_essential_vertices: 0,
        }
    }

    #[allow(unused)]
    pub fn skeleton_of_cube_with_4_marked_points() -> SimpleGraph {
        let mut raw_graph = RawSimpleGraph::new(32);
        // add deleted edges
        raw_graph.add_edge(0, 11);
        raw_graph.add_edge(0, 19);
        raw_graph.add_edge(3, 24);
        raw_graph.add_edge(6, 29);
        raw_graph.add_edge(14, 31);
        // add non-order-respecting edges in the maximal tree
        raw_graph.add_edge(9, 12);
        raw_graph.add_edge(17, 20);
        raw_graph.add_edge(22, 25);
        raw_graph.add_edge(27, 30);
        // add order-respecting edges in the maximal tree
        (0..32-1).filter(|&i| i!=11 && i!=19 && i!=24 && i!=29 ).for_each(|i| raw_graph.add_edge(i, i+1));

        raw_graph.is_maximal_tree_chosen = true;

        SimpleGraph {
            data: raw_graph,
            vertices_list: vec![SimpleVertex::Essential("".to_string()); 32],
            n_essential_vertices: 0,
        }
    }
}

#[derive(Clone)]
pub struct WeightedSimpleGraph<T: PID> {
    data: Matrix<T>,
    vertices_labeling: Vec<String>,
}

impl<T: PID + Copy> WeightedSimpleGraph<T> {
    pub fn from(input: Vec<(String, String, T)>)-> Self {
        let mut vertices_labeling: Vec<String> = Vec::new();
        let mut weights = Vec::new();
        for (v1, v2, weight) in input {
            assert!( v1!=v2, "loop is not allowed in simple graphs" );
            let v1_idx =  vertices_labeling.iter().position( |v| v == &v1 );
            let v2_idx =  vertices_labeling.iter().position( |v| v == &v2 );
            // assert!( !(v1_idx!=None && v2_idx!=None), "Simple Graph cannot have more than one edge between the identical pair of eddges." );

            let v1_idx = match v1_idx {
                Some(idx) => idx,
                None => {
                    vertices_labeling.push( v1 );
                    vertices_labeling.len()-1
                },
            };

            let v2_idx = match v2_idx {
                Some(idx) => idx,
                None => {
                    vertices_labeling.push( v2 );
                    vertices_labeling.len()-1
                },
            };

            weights.push( (v1_idx, v2_idx, weight) )
        }

        // register the weights
        let mut data = Matrix::zero(vertices_labeling.len(), vertices_labeling.len());
        for (i, j, weight) in weights.into_iter() {
            data[(i,j)] = weight;
            data[(j,i)] = weight;
        }

        WeightedSimpleGraph {
            data: data,
            vertices_labeling: vertices_labeling,
        }
    }

    pub fn laplacian(self) -> Matrix<T> {
        let mut laplacian = self.data;
        let n = laplacian.size().0;
        for i in 0..n {
            laplacian[(i,i)] = (0..n).map( |j| laplacian[(i,j)] ).sum();
        }

        for i in 0..n {
            for j in 0..n {
                if i==j {
                    continue;
                }
                laplacian[(i,j)] *= -T::one();
            }
        }

        laplacian
    }

    pub fn n_vertices(&self) -> usize {
        self.vertices_labeling.len()
    }

    pub fn vertices(&self) -> &Vec<String> {
        &self.vertices_labeling
    }

    pub fn data(&self) -> &Matrix<T> {
        &self.data
    }
}

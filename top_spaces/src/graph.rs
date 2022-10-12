use util::BinaryPool;
use algebra::commutative::rational::*;
use std::fmt;

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

#[allow(unused)]
pub struct SimpleGraph {
    data: BinaryPool,
    vertices_list: Vec<SimpleVertex>,
    n_essential_vertices: usize,
    n_edges: usize,
}


// non public interfaces of SimpleGraph
#[allow(unused)]
impl SimpleGraph {
    fn new() -> SimpleGraph {
        SimpleGraph {
            data: BinaryPool::new(0),
            vertices_list: Vec::new(),
            n_essential_vertices: 0,
            n_edges: 0,
        }
    }

    fn raw_graph_from(graph: &Graph ) -> SimpleGraph {
        // allocation
        let n_vertices = 
            graph.n_vertices() 
            + graph.n_non_loop_edges() - graph.n_non_loop_edges_up_to_multiplicity()
            + graph.n_loops() * 2;
        let data = BinaryPool::new(n_vertices * (n_vertices-1)/2);

        let vertices_list:Vec<_> = (0..graph.n_vertices()).map(|i| SimpleVertex::Essential(graph.vertices_labeling[i].clone())).collect();

        let mut simple_graph = SimpleGraph {
            data: data,
            vertices_list: vertices_list,
            n_essential_vertices: graph.n_vertices(),
            n_edges: 0,
        };

        // Adding vertices
        for i in 0..simple_graph.n_essential_vertices() {
            // Iterating over non-loop edges
            for j in i+1..simple_graph.n_essential_vertices() {
                // if there is no edge on the two vertex, we simply continue to the next.
                if graph.data[i][j].len() == 0 {
                    continue;
                }

                simple_graph.add(i,j);

                // if the edge does not have any multiplicity, continue to the next.
                if graph.data[i][j].len() == 1 {
                    continue;
                }

                // else, there are multiple edges on the same pair of vertices.
                for (k, edge) in graph.data[i][j].iter().skip(1).enumerate() {
                    let mid_vertex = VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[j].clone(), k, 1.over(2));
                    let mid_vertex = SimpleVertex::Redundant(mid_vertex);
                    simple_graph.vertices_list.push(mid_vertex);

                    let mid_idx = simple_graph.vertices_list.len();
                    simple_graph.add(i, mid_idx);
                    simple_graph.add(j, mid_idx);
                }
            }

            // Finally, we process loops, i.e. edges of the form [i, i]
            for (k, edge) in graph.data[i][i].iter().enumerate() {
                let first_vertex = VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[i].clone(), k, 1.over(3));
                let second_vertex = VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[i].clone(), k, 2.over(3));
                let first_vertex = SimpleVertex::Redundant(first_vertex);
                let second_vertex = SimpleVertex::Redundant(second_vertex);
                simple_graph.vertices_list.push(first_vertex);
                simple_graph.vertices_list.push(second_vertex);

                let first_idx = simple_graph.vertices_list.len()-2;
                let second_idx = simple_graph.vertices_list.len()-1;
                simple_graph.add(i, first_idx);
                simple_graph.add(i, second_idx);
                simple_graph.add(first_idx, second_idx);
            }
        }

        simple_graph
    }

    fn raw_subdivided_graph_from(graph: &Graph, n_marked_points: usize) -> Self {
        // When n_marked_points <= 2, we can just return the simple graph
        if n_marked_points <= 2 { return SimpleGraph::from(graph) }

        let n_vertices = 
            graph.n_vertices() 
            + graph.n_non_loop_edges() * (n_marked_points-2)
            + graph.n_loops() * n_marked_points;
        let data = BinaryPool::new(n_vertices * (n_vertices-1) / 2 );

        // Adding vertices
        let vertices_list:Vec<_> = (0..graph.n_vertices()).map(|i| SimpleVertex::Essential(graph.vertices_labeling[i].clone())).collect();

        let mut subdivided_graph = SimpleGraph {
            data: data,
            vertices_list: vertices_list,
            n_essential_vertices: graph.n_vertices(),
            n_edges: 0,
        };

        // for each essential vertex...
        for i in 0..subdivided_graph.n_essential_vertices() {
            // Iterating over non-loop edges
            for j in i+1..subdivided_graph.n_essential_vertices() {
                if graph.data[i][j].len() == 0 {
                    continue;
                }

                for (k, edge) in graph.data[i][j].iter().enumerate() {
                    subdivided_graph.add(i, subdivided_graph.n_vertices());
                    subdivided_graph.vertices_list.push(SimpleVertex::Redundant(VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[j].clone(), k, 1.over(n_marked_points-1))));

                    for l in 2..n_marked_points-1 {
                        subdivided_graph.add(subdivided_graph.n_vertices()-1, subdivided_graph.n_vertices());
                        subdivided_graph.vertices_list.push(SimpleVertex::Redundant(VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[j].clone(), k, l.over(n_marked_points-1))));
                    }

                    subdivided_graph.add(j, subdivided_graph.n_vertices()-1);
                }
            }

            // Finally, we process loops, i.e. edges of the form [i, i]
            for (k, edge) in graph.data[i][i].iter().enumerate() {
                subdivided_graph.add(i, subdivided_graph.n_vertices());
                subdivided_graph.vertices_list.push(SimpleVertex::Redundant(VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[i].clone(), k, 1.over(n_marked_points+1))));

                for l in 2..=n_marked_points {
                    subdivided_graph.add(subdivided_graph.n_vertices()-1, subdivided_graph.n_vertices());
                    subdivided_graph.vertices_list.push(SimpleVertex::Redundant(VertexVector::new(graph.vertices_labeling[i].clone(), graph.vertices_labeling[i].clone(), k, l.over(n_marked_points+1))));
                }

                subdivided_graph.add(i, subdivided_graph.n_vertices()-1);
            }
        }

        subdivided_graph
    }

    // function returing the raw index for self.data. The condition that i < j is required.
    fn raw_idx(i: usize, j: usize) -> usize {
        debug_assert!(i < j, "The condition i < j is required, but not satisfied: i={} but j={}", i, j);
        j*(j-1)/2 + i
    }

    fn add(&mut self, i: usize, j: usize) {
        self.data.add(SimpleGraph::raw_idx(i, j));
        self.n_edges+=1;
    }

    fn subtract(&mut self, i: usize, j: usize) {
        self.data.subtract(SimpleGraph::raw_idx(i, j));
        self.n_edges-=1;
    }

    fn contains(&self, i: usize, j: usize) -> bool {
        self.data.contains(SimpleGraph::raw_idx(i, j))
    }

    fn erase_edges_adjacent_to(&mut self, i: usize) {
        let edges_to_delete: Vec<_> = self.adjacent_vertices_iter(i).map(|j| if i<j {(i,j)} else {(j,i)} ).collect();
        edges_to_delete.iter().for_each(|&(i, j)| self.subtract(i, j) );
    }

    fn swap(&mut self, i: usize, j: usize) {
        debug_assert!(i < j && j < self.n_vertices(), "Two condition i < j and j < self.n_vertices() are required, but not satisfied: i={}, j={}, self.n_vertices()={}.", i, j, self.n_vertices());

        let isotopy = self.contains(i,j);

        let tmp_i: Vec<_> = self.adjacent_vertices_iter(i).filter(|&v| v != j ).map(|v| if j < v {(j, v)} else {(v, j)} ).collect();
        let tmp_j: Vec<_> = self.adjacent_vertices_iter(j).filter(|&v| v != i ).map(|v| if i < v {(i, v)} else {(v, i)} ).collect();
        self.erase_edges_adjacent_to(i);
        self.erase_edges_adjacent_to(j);
        tmp_i.iter().for_each(|&(i,j)| self.add(i,j));
        tmp_j.iter().for_each(|&(i,j)| self.add(i,j));

        if isotopy==true { self.add(i,j); };

        self.vertices_list.swap(i, j);
    }

    fn degree_of(&self, idx: usize) -> usize {
        if idx >= self.vertices_list.len() { return 1; };

        (0..idx)
        .map(|j| {
            if self.contains(j, idx) {1} else {0}
        })
        .fold(0, |sum, n| sum + n)
    }

    fn choose_maximal_tree(&mut self) {
        for i in 0..self.n_vertices()-1 {
            for j in (0..=i).rev() {
                let result = self.fetch_right_and_swap(j, i+1);
                if result { break; }
            }
        }
    }


    // fetch the vertex in the right. This swaps it with i and return true if it finds it. If it did not find any vertex, then this returns false.
    fn fetch_right_and_swap( &mut self, i: usize, j: usize ) -> bool {
        let op_k = self.adjacent_vertices_iter(i).filter(|&k| k > j).next();
        match op_k {
            Some(k) => { self.swap(j, k); true },
            None => self.contains(i, j) ,
        }
    }

    fn are_adjacent(&self, v1: usize, v2: usize ) -> bool {
        if v1 < v2 {
            self.contains(v1, v2)
        } else if v1 > v2 {
            self.contains(v2, v1)
        } else {
            false
        }
    }
}


// public interfaces of SimpleGraph
impl SimpleGraph {
    pub fn n_vertices(&self) -> usize{
        self.vertices_list.len()
    }

    pub fn n_edges(&self) -> usize {
        self.n_edges
    }

    pub fn n_essential_vertices(&self) -> usize{
        self.n_essential_vertices
    }

    pub fn from(graph: &Graph) -> SimpleGraph {
        let mut simple_graph = Self::raw_graph_from(graph);
        simple_graph.choose_maximal_tree();
        simple_graph
    }

    pub fn subdivided_graph_from(graph: &Graph, n_marked_points: usize) -> SimpleGraph {
        let mut subdivided_graph = Self::raw_subdivided_graph_from(graph, n_marked_points);
        subdivided_graph.choose_maximal_tree();
        subdivided_graph
    }
}

pub struct AdjacentVerticesIter<'a> {
    graph: &'a SimpleGraph,
    curr: usize,
    vertex: usize,
}

impl<'a> AdjacentVerticesIter<'a> {
    fn new(graph: &'a SimpleGraph, vertex: usize) -> AdjacentVerticesIter<'a> {
        AdjacentVerticesIter { graph: graph, curr: 0, vertex: vertex }
    }
}

impl<'a> Iterator for AdjacentVerticesIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr >= self.graph.n_vertices() { return None; }

        let mut out = None;

        if self.vertex < self.curr {
            if self.graph.data.contains(SimpleGraph::raw_idx(self.vertex, self.curr)) {
                out = Some(self.curr);
            }
        } else if self.vertex > self.curr {
            if self.graph.data.contains(SimpleGraph::raw_idx(self.curr, self.vertex)) {
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

pub struct VertexIter<'a> {
    _graph: &'a SimpleGraph,
    iter: std::ops::Range<usize>,
}

impl<'a> VertexIter<'a> {
    fn new(graph: &'a SimpleGraph) -> Self {
        VertexIter {
            iter: (0..graph.n_vertices()),
            _graph: graph,
        }
    }
}

impl<'a> Iterator for VertexIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        return self.iter.next();
    }
}

pub struct EdgeIter<'a> {
    graph: &'a SimpleGraph,
    curr: (usize, usize),
    is_done: bool,
}

impl<'a> EdgeIter<'a> {
    fn new(graph: &'a SimpleGraph) -> Self {
        let edge = (0, 1);
        EdgeIter {
            graph: graph,
            curr: edge,
            is_done: graph.n_vertices() < 2,
        }
    }
}

impl<'a> Iterator for EdgeIter<'a> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_done {
            return None;
        };

        let out = self.curr;
        
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
impl<'a> SimpleGraph {
    pub fn adjacent_vertices_iter(&'a self, idx: usize) -> AdjacentVerticesIter<'a> {
        AdjacentVerticesIter::new(self, idx)
    }

    pub fn edge_iter(&'a self) -> EdgeIter<'a> {
        EdgeIter::new(self)
    }

    pub fn vertex_iter(&'a self) -> VertexIter<'a> {
        VertexIter::new(self)
    }
}


#[cfg(test)]
mod simple_graph_data_manip_test {
    use util::BinaryPool;
    use super::SimpleGraph;

    #[test]
    fn iter_and_erase_test() {
        let graph = graph! {
            ["v0", "v1"],
            ["v1", "v2"],
            ["v1", "v3"],
            ["v0", "v3"],
            ["v0", "v4"]
        };
        let mut simple_graph = SimpleGraph::raw_graph_from(&graph);
        assert_eq!(simple_graph.adjacent_vertices_iter(0).collect::<Vec<_>>(), vec![1,3,4]);
        simple_graph.erase_edges_adjacent_to(0);
        assert_eq!(simple_graph.adjacent_vertices_iter(0).collect::<Vec<_>>(), vec![]);
    }

    #[test]
    fn swap_test() {
        let graph = graph! {
            ["v0", "v1"],
            ["v1", "v2"],
            ["v1", "v3"],
            ["v0", "v3"],
            ["v0", "v4"]
        };
        let mut subdivided_graph = SimpleGraph::raw_subdivided_graph_from(&graph, 3);
        let mut answer_data = BinaryPool::new(45);
        {
            answer_data.add(SimpleGraph::raw_idx(0, 5));
            answer_data.add(SimpleGraph::raw_idx(1, 5));
            answer_data.add(SimpleGraph::raw_idx(0, 6));
            answer_data.add(SimpleGraph::raw_idx(3, 6));
            answer_data.add(SimpleGraph::raw_idx(0, 7));
            answer_data.add(SimpleGraph::raw_idx(4, 7));
            answer_data.add(SimpleGraph::raw_idx(1, 8));
            answer_data.add(SimpleGraph::raw_idx(2, 8));
            answer_data.add(SimpleGraph::raw_idx(1, 9));
            answer_data.add(SimpleGraph::raw_idx(3, 9));
        }
        assert_eq!(subdivided_graph.data, answer_data);

        subdivided_graph.swap(1, 5);
        let mut answer_data = BinaryPool::new(45);
        {
            answer_data.add(SimpleGraph::raw_idx(0, 1));
            answer_data.add(SimpleGraph::raw_idx(0, 6));
            answer_data.add(SimpleGraph::raw_idx(0, 7));
            answer_data.add(SimpleGraph::raw_idx(1, 5));
            answer_data.add(SimpleGraph::raw_idx(2, 8));
            answer_data.add(SimpleGraph::raw_idx(3, 6));
            answer_data.add(SimpleGraph::raw_idx(3, 9));
            answer_data.add(SimpleGraph::raw_idx(4, 7));
            answer_data.add(SimpleGraph::raw_idx(5, 8));
            answer_data.add(SimpleGraph::raw_idx(5, 9));
        }
        assert_eq!(subdivided_graph.data, answer_data);

        subdivided_graph.swap(2, 5);
        subdivided_graph.swap(3, 8);
        subdivided_graph.swap(4, 5);
        subdivided_graph.swap(5, 9);
        subdivided_graph.swap(6, 8);
        subdivided_graph.swap(7, 8);
        let mut answer_data = BinaryPool::new(45);
        {
            answer_data.add(SimpleGraph::raw_idx(0, 1));
            answer_data.add(SimpleGraph::raw_idx(0, 7));
            answer_data.add(SimpleGraph::raw_idx(0, 8));
            answer_data.add(SimpleGraph::raw_idx(1, 2));
            answer_data.add(SimpleGraph::raw_idx(2, 3));
            answer_data.add(SimpleGraph::raw_idx(2, 5));
            answer_data.add(SimpleGraph::raw_idx(3, 4));
            answer_data.add(SimpleGraph::raw_idx(5, 6));
            answer_data.add(SimpleGraph::raw_idx(6, 7));
            answer_data.add(SimpleGraph::raw_idx(8, 9));
        }
        assert_eq!(subdivided_graph.data, answer_data);
    }
}


#[cfg(test)]
mod simple_graph_test {
    use util::BinaryPool;
    use super::SimpleGraph;
    #[test]
    fn init_test() {
        let graph = graph!{
            ["v0", "v0"],
            ["v0", "v1"],
            ["v0", "v1"],
            ["v1", "v1"]
        };

        let simple_graph = SimpleGraph::raw_graph_from(&graph);
        assert_eq!(simple_graph.n_vertices(), 7);
    }

    #[test]
    fn subdivision_test() {
        // test 1
        let graph = graph!{
            ["v0", "v1"],
            ["v1", "v1"]
        };

        let mut subdivided_graph = SimpleGraph::raw_subdivided_graph_from(&graph, 3);
        assert_eq!(subdivided_graph.n_vertices(), 6);
        let mut answer_data = BinaryPool::new(15);
        {
            answer_data.add(SimpleGraph::raw_idx(0, 2));
            answer_data.add(SimpleGraph::raw_idx(1, 2));
            answer_data.add(SimpleGraph::raw_idx(1, 3));
            answer_data.add(SimpleGraph::raw_idx(1, 5));
            answer_data.add(SimpleGraph::raw_idx(3, 4));
            answer_data.add(SimpleGraph::raw_idx(4, 5));
        }
        assert_eq!(subdivided_graph.data, answer_data);

        subdivided_graph.choose_maximal_tree();
        let mut answer_data = BinaryPool::new(15);
        {
            answer_data.add(SimpleGraph::raw_idx(0, 1));
            answer_data.add(SimpleGraph::raw_idx(1, 2));
            answer_data.add(SimpleGraph::raw_idx(2, 3));
            answer_data.add(SimpleGraph::raw_idx(3, 4));
            answer_data.add(SimpleGraph::raw_idx(4, 5));
            answer_data.add(SimpleGraph::raw_idx(2, 5));
        }
        assert_eq!(subdivided_graph.data, answer_data);


        // test2
        let graph = graph!{
            ["v0", "v0"],
            ["v0", "v1"],
            ["v0", "v1"],
            ["v1", "v1"]
        };

        let mut subdivided_graph = SimpleGraph::raw_subdivided_graph_from(&graph, 3);
        assert_eq!(subdivided_graph.n_vertices(), 10);
        let mut answer_data = BinaryPool::new(45);
        {
            answer_data.add(SimpleGraph::raw_idx(0, 2));
            answer_data.add(SimpleGraph::raw_idx(0, 3));
            answer_data.add(SimpleGraph::raw_idx(0, 4));
            answer_data.add(SimpleGraph::raw_idx(0, 6));
            answer_data.add(SimpleGraph::raw_idx(1, 2));
            answer_data.add(SimpleGraph::raw_idx(1, 3));
            answer_data.add(SimpleGraph::raw_idx(1, 7));
            answer_data.add(SimpleGraph::raw_idx(1, 9));
            answer_data.add(SimpleGraph::raw_idx(4, 5));
            answer_data.add(SimpleGraph::raw_idx(5, 6));
            answer_data.add(SimpleGraph::raw_idx(7, 8));
            answer_data.add(SimpleGraph::raw_idx(8, 9));
        }
        assert_eq!(subdivided_graph.data, answer_data);

        subdivided_graph.choose_maximal_tree();
        let mut answer_data = BinaryPool::new(45);
        {
            answer_data.add(SimpleGraph::raw_idx(0, 1));
            answer_data.add(SimpleGraph::raw_idx(1, 2));
            answer_data.add(SimpleGraph::raw_idx(2, 3));
            answer_data.add(SimpleGraph::raw_idx(3, 4));
            answer_data.add(SimpleGraph::raw_idx(4, 5));
            answer_data.add(SimpleGraph::raw_idx(2, 5));
            answer_data.add(SimpleGraph::raw_idx(2, 6));
            answer_data.add(SimpleGraph::raw_idx(0, 6));
            answer_data.add(SimpleGraph::raw_idx(0, 7));
            answer_data.add(SimpleGraph::raw_idx(7, 8));
            answer_data.add(SimpleGraph::raw_idx(8, 9));
            answer_data.add(SimpleGraph::raw_idx(0, 9));
        }
        assert_eq!(subdivided_graph.data, answer_data);

        // test 3
        let graph = graph! {
            ["v0", "v1"],
            ["v1", "v2"],
            ["v1", "v3"],
            ["v0", "v3"],
            ["v0", "v4"]
        };
        let mut subdivided_graph = SimpleGraph::raw_subdivided_graph_from(&graph, 3);
        assert_eq!(subdivided_graph.n_vertices(), 10);
        let mut answer_data = BinaryPool::new(45);
        {
            answer_data.add(SimpleGraph::raw_idx(0, 5));
            answer_data.add(SimpleGraph::raw_idx(1, 5));
            answer_data.add(SimpleGraph::raw_idx(0, 6));
            answer_data.add(SimpleGraph::raw_idx(3, 6));
            answer_data.add(SimpleGraph::raw_idx(0, 7));
            answer_data.add(SimpleGraph::raw_idx(4, 7));
            answer_data.add(SimpleGraph::raw_idx(1, 8));
            answer_data.add(SimpleGraph::raw_idx(2, 8));
            answer_data.add(SimpleGraph::raw_idx(1, 9));
            answer_data.add(SimpleGraph::raw_idx(3, 9));
        }
        assert_eq!(subdivided_graph.data, answer_data);

        subdivided_graph.choose_maximal_tree();
        let mut answer_data = BinaryPool::new(45);
        {
            answer_data.add(SimpleGraph::raw_idx(0, 1));
            answer_data.add(SimpleGraph::raw_idx(0, 7));
            answer_data.add(SimpleGraph::raw_idx(0, 8));
            answer_data.add(SimpleGraph::raw_idx(1, 2));
            answer_data.add(SimpleGraph::raw_idx(2, 3));
            answer_data.add(SimpleGraph::raw_idx(2, 5));
            answer_data.add(SimpleGraph::raw_idx(3, 4));
            answer_data.add(SimpleGraph::raw_idx(5, 6));
            answer_data.add(SimpleGraph::raw_idx(6, 7));
            answer_data.add(SimpleGraph::raw_idx(8, 9));
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
        for (i, (v, w)) in self.edge_iter().enumerate() {
            // printing each cell
            write!(f, "[{}, {}]", self.vertices_list[v], self.vertices_list[w])?;
            if i != self.n_edges()-1 {
                write!(f, ", ")?;
            };
        };

        write!(f, "\n")
    }
}


#[cfg(test)]
mod simple_graph_utility_test {
    use super::SimpleGraph;

    #[test]
    fn display_test() {
        let graph = graph!{
            ["v0", "v1"],
            ["v1", "v1"]
        };

        let mut subdivided_graph = SimpleGraph::raw_subdivided_graph_from(&graph, 3);
        subdivided_graph.choose_maximal_tree();
        let description = subdivided_graph.to_string();
        let ans = 
            "vertices: v0, (1/2)(v0, v1), v1, (3/4)(v1, v1), (1/2)(v1, v1), (1/4)(v1, v1)\nedges: [v0, (1/2)(v0, v1)], [(1/2)(v0, v1), v1], [v1, (3/4)(v1, v1)], [(3/4)(v1, v1), (1/2)(v1, v1)], [(1/2)(v1, v1), (1/4)(v1, v1)], [v1, (1/4)(v1, v1)]\n"
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
            (0,1),
            (1,2),
            (2,3),
            (1,3)
        ];
        assert_eq!(edges, answer);

        let subdivided_graph = SimpleGraph::subdivided_graph_from(&graph, 3);
        let edges:Vec<_> = subdivided_graph.edge_iter().collect();
        let answer = vec![
            (0,1),
            (1,2),
            (2,3),
            (3,4),
            (4,5),
            (2,5)
        ];
        assert_eq!(edges, answer);
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

    pub fn is_in_maximal_tree(self) -> bool {
        match self {
            Self::Edge(x, y) => y - x == 1,
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
            Self::Vertex(_) => panic!("Not an edge"),
        }
    }

    pub fn vertex(self) -> usize {
        match self {
            Self::Edge(_, _) => panic!("Not a vertex"),
            Self::Vertex(x) => x,
        }
    }
}


pub struct FactorIter<'a> {
    v_iter: VertexIter<'a>,
    e_iter: EdgeIter<'a>,
}

impl<'a> FactorIter<'a> {
    fn new(graph: &'a SimpleGraph) -> FactorIter<'a>{
        FactorIter { 
            v_iter: graph.vertex_iter(), 
            e_iter: graph.edge_iter(),
        }
    }
}

impl<'a> Iterator for FactorIter<'a> {
    type Item = CubeFactor;

    fn next(&mut self) -> Option<Self::Item> {
        let out = self.v_iter.next();
        match out {
            Some(x) => Some(CubeFactor::Vertex(x)),
            None => match self.e_iter.next() {
                Some((x, y)) => Some(CubeFactor::Edge(x,y)),
                None => None,
            },
        }
    }
}

impl<'a> SimpleGraph {
    pub fn factor_iter(&'a self) -> FactorIter<'a> {
        FactorIter::new(self)
    }
}

use std::cmp::Ordering;
impl PartialOrd for CubeFactor {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // if the two cells are not dijoint, then they are incomparable.
        if !self.is_disjoint_from(*other) {
            return None;
        }
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
            CubeFactor::Edge(x, y) => write!(f, "({}, {})", x, y),
            CubeFactor::Vertex(x) => write!(f, "({})", x),
        }
    }
}
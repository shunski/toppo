use util::BinaryPool;
use algebra::commutative::rational::*;

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

    n_edges: usize,
    n_edges_up_to_multiplicity: usize,
    n_loops: usize,
}


#[allow(unused)]
impl Graph {
    fn from( edges: Vec<(String, String)> ) -> Graph {
        let n_edges = edges.len();
        let mut n_edges_up_to_multiplicity = edges.len();
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
            let (v_idx, w_idx) = if v_idx > w_idx {
                (v_idx, w_idx)
            } else {
                (w_idx, v_idx)
            };

            data[v_idx][w_idx].push(Edge::new());

            // update some members
            if data[v_idx][w_idx].len() > 1 {
                n_edges_up_to_multiplicity -= 1;
            } 
            if v_idx == w_idx {
                n_loops += 1;
            }
        });

        Graph {
            data: data,
            vertices_labeling: vertices_labeling,
            n_edges: n_edges,
            n_edges_up_to_multiplicity: n_edges_up_to_multiplicity,
            n_loops: n_loops,
        }
    }

    pub fn n_vertices(&self) -> usize { self.vertices_labeling.len() }

    pub fn n_edges(&self) -> usize { 
        self.n_edges
    }

    pub fn n_edges_up_to_multiplicity(&self) -> usize {
        self.n_edges_up_to_multiplicity
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
            use crate::graph::Graph;
            Graph::from(v)
        }
    };
}


mod graph_test {
    #[test]
    fn init_test() {
        let graph = graph!{
            ["v0", "v0"],
            ["v0", "v1"],
            ["v1", "v0"],
            ["v1", "v1"]
        };
        assert_eq!(graph.n_vertices(), 2);
        assert_eq!(graph.n_edges(), 4);
        assert_eq!(graph.n_edges_up_to_multiplicity(), 3);
        assert_eq!(graph.n_loops(), 2);
        assert_eq!(graph.data[0][0].len(), 1);
        assert_eq!(graph.data[1][0].len(), 2);
        assert_eq!(graph.data[1][1].len(), 1);
    }
}

#[allow(unused)]
#[derive(Copy, Clone)]
pub struct VertexVector {
    start: usize,
    end: usize,
    direction: usize,
    pos: Rational,
}

#[allow(unused)]
impl VertexVector {
    fn new(start: usize, end: usize, direction: usize, pos: Rational ) -> VertexVector {
        let (start, end) = if start >= end {
            (start, end)
        } else {
            (end, start)
        };
        VertexVector { start: start, end: end, direction: direction, pos: pos }
    }
}

#[allow(unused)]
#[derive(Copy, Clone)]
pub enum SimpleVertex {
    Essential(usize),
    Inessential(VertexVector),
}

#[allow(unused)]
pub struct SimpleGraph {
    data: BinaryPool,
    vertices_list: Vec<SimpleVertex>,
    n_essential_vertices: usize,
}

#[allow(unused)]
impl SimpleGraph {
    fn new() -> SimpleGraph {
        SimpleGraph {
            data: BinaryPool::new(0),
            vertices_list: Vec::new(),
            n_essential_vertices: 0,
        }
    }

    pub fn from(graph: &Graph ) -> SimpleGraph {
        // allocation
        let n_vertices = 
            graph.n_vertices() 
            + graph.n_edges() - graph.n_edges_up_to_multiplicity()
            + graph.n_loops()*2;
        let mut data = BinaryPool::new(n_vertices * (n_vertices-1));

        // Adding vertices
        let mut vertices_list:Vec<_> = (0..graph.n_vertices()).map(|i| SimpleVertex::Essential(i)).collect();
        for i in 0..graph.n_vertices() {
            // Iterating over non-loop edges
            for j in 0..i {
                // if there is no edge on the two vertex, we simply continue to the next.
                if graph.data[i][j].len() == 0 {
                    continue;
                }

                data.add(Self::raw_idx(i,j));

                // if the edge does not have any multiplicity, continue to the next.
                if graph.data[i][j].len() == 1 {
                    continue;
                }

                println!("multiple edges");

                // else, there are multiple edges on the same pair of vertices.
                for (k, edge) in graph.data[i][j].iter().skip(1).enumerate() {
                    let mid_vertex = VertexVector::new(i, j, k, 1.over(2));
                    let mid_vertex = SimpleVertex::Inessential(mid_vertex);
                    vertices_list.push(mid_vertex);

                    let mid_idx = vertices_list.len();
                    data.add(Self::raw_idx(mid_idx, i));
                    data.add(Self::raw_idx(mid_idx, j));
                }
            }

            // Finally, we process loops, i.e. edges of the form [i, i]
            for (k, edge) in graph.data[i][i].iter().enumerate() {
                println!("adding loops");
                let first_vertex = VertexVector::new(i, i, k, 1.over(3));
                let second_vertex = VertexVector::new(i, i, k, 2.over(3));
                let first_vertex = SimpleVertex::Inessential(first_vertex);
                let second_vertex = SimpleVertex::Inessential(second_vertex);
                vertices_list.push(first_vertex);
                vertices_list.push(second_vertex);

                let first_idx = vertices_list.len()-1;
                let second_idx = vertices_list.len();
                data.add(Self::raw_idx(first_idx, i));
                data.add(Self::raw_idx(second_idx, i));
                data.add(Self::raw_idx(second_idx, first_idx));
            }
        }

        SimpleGraph {
            data: data,
            vertices_list: vertices_list,
            n_essential_vertices: graph.n_vertices(),
        }
    }

    fn raw_idx(i: usize, j: usize) -> usize {
        if i <= j || i==0 { panic!("i must be greater than j, and i must be greater than 0: i={} but j={}", i, j); }
        i*(i-1)/2 + j
    }

    pub fn n_vertices(&self) -> usize{
        self.vertices_list.len()
    }

    fn degree_of(&self, idx: usize) -> usize {
        if idx >= self.vertices_list.len() { return 1; };

        (0..idx)
        .map(|j| {
            if self.data.contains(Self::raw_idx(idx,j)) {1} else {0}
        })
        .fold(0, |sum, n| sum + n)
    }

    pub fn subdivided_graph_from(graph: &Graph, n_marked_points: usize) -> Self {
        // When n_marked_points <= 2, we can just return the simple graph
        if n_marked_points <= 2 { return SimpleGraph::from(graph) }

        let n_vertices = 
            graph.n_vertices() 
            + graph.n_edges() - graph.n_edges_up_to_multiplicity() * (n_marked_points-2)
            + graph.n_loops() * n_marked_points;
        let mut data = BinaryPool::new(n_vertices * (n_vertices-1));

        // Adding vertices
        let mut vertices_list:Vec<_> = (0..graph.n_vertices()).map(|i| SimpleVertex::Essential(i)).collect();

        let n_essential_vertices = graph.n_vertices();

        // for each essential vertex...
        for i in 0..n_essential_vertices {
            // Iterating over non-loop edges
            for j in 0..i {
                if graph.data[i][j].len() == 0 {
                    continue;
                }

                for (k, edge) in graph.data[i][j].iter().enumerate() {
                    data.add(Self::raw_idx(vertices_list.len(), i));
                    vertices_list.push(SimpleVertex::Inessential(VertexVector::new(i, j, k, 1.over(n_marked_points-1))));

                    for l in 2..n_marked_points-1 {
                        data.add(Self::raw_idx(vertices_list.len(), vertices_list.len()-1));
                        vertices_list.push(SimpleVertex::Inessential(VertexVector::new(i, j, k, l.over(n_marked_points-1))));
                    }

                    data.add(Self::raw_idx(vertices_list.len(), j));
                }
            }

            // Finally, we process loops, i.e. edges of the form [i, i]
            for (k, edge) in graph.data[i][i].iter().enumerate() {
                data.add(Self::raw_idx(vertices_list.len(), i));
                vertices_list.push(SimpleVertex::Inessential(VertexVector::new(i, i, k, 1.over(n_marked_points+1))));

                for l in 2..n_marked_points+1 {
                    data.add(Self::raw_idx(vertices_list.len(), vertices_list.len()-1));
                    vertices_list.push(SimpleVertex::Inessential(VertexVector::new(i, i, k, l.over(n_marked_points+1))));
                }

                data.add(Self::raw_idx(vertices_list.len(), i));
            }
        }

        SimpleGraph {
            data: data,
            vertices_list: vertices_list,
            n_essential_vertices: n_essential_vertices,
        }
    }

    fn adjacent(&self, v1: usize, v2: usize ) -> bool {
        if v1 > v2 {
            self.data.contains(Self::raw_idx(v1, v2))
        } else if v1 < v2 {
            self.data.contains(Self::raw_idx(v2, v1))
        } else {
            false
        }
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

        if self.curr < self.vertex {
            if self.graph.data.contains(SimpleGraph::raw_idx(self.vertex, self.curr)) {
                out = Some(self.curr-1);
            }
        } else if self.curr > self.vertex {
            if self.graph.data.contains(SimpleGraph::raw_idx(self.curr, self.vertex)) {
                out = Some(self.curr-1);
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

#[allow(unused)]
impl<'a> SimpleGraph {
    fn adjacent_vertices_iter(&'a self, idx: usize) -> AdjacentVerticesIter<'a> {
        AdjacentVerticesIter::new(self, idx)
    }
}


#[cfg(test)]
mod simple_graph_test {
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
        assert_eq!(simple_graph.n_vertices(), 7);
    }

    #[test]
    fn subdivision_test() {
        let graph = graph!{
            ["v0", "v0"],
            ["v0", "v1"],
            ["v0", "v1"],
            ["v1", "v1"]
        };

        let subdivided_graph = SimpleGraph::subdivided_graph_from(&graph, 3);
        assert_eq!(subdivided_graph.n_vertices(), 10);
    }
}
use crate::Space;
use algebra::commutative::Field;
use crate::graph::*;

#[derive(Clone, PartialEq, Eq)]
struct CubicalCell {
    dim: usize,
    data: Vec<(usize, usize)>,
}

impl CubicalCell {
    fn new (mut data: Vec<(usize, usize)>) -> Self {
        data.iter_mut().for_each(|(mut x, mut y)| if x>y {let tmp = x; x = y; y = tmp; } );
        let dim = data.iter().map(|&(x, y)| if x==y {0} else {1} ).fold(0, |sum, n| sum+n );
        CubicalCell {
            data: data,
            dim: dim
        }
    }
}

struct CubicalCplx {
    edge_pool: SimpleGraph,
    maximal_cells: Vec<Vec<CubicalCell>>,
}

impl CubicalCplx {
    fn DC(graph: Graph, n_marked_points: usize) -> CubicalCplx {
        let subdivided_graph = SimpleGraph::subdivided_graph_from(&graph, n_marked_points);
        
        let mut maximal_cells = vec![Vec::new(); n_marked_points + 1];

        for dim in (0..=n_marked_points).rev() {
            
        }

        CubicalCplx {
            edge_pool: subdivided_graph,
            maximal_cells: maximal_cells,
        }
    }
}
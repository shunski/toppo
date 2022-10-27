use topo_spaces;
use topo_spaces::graph;
use topo_spaces::cubical::CubicalCplx;
use topo_spaces::*;

fn main() {
    // Non-hamiltonian graph #1
    let non_hamilonian = graph! {
        ["v1", "v2"],
        ["v1", "v3"],
        ["v1", "v4"],
        ["v2", "v3"],
        ["v2", "v4"],
        ["v3", "v4"],

        ["v1", "v5"],
        ["v1", "v6"],
        ["v5", "v6"]
    };
    let n_marked_points = 2;
    let cplx = CubicalCplx::udc(&non_hamilonian, n_marked_points);
    println!("Non_hamilonian graph #1 with {} marked points: \n{}", n_marked_points, H!(cplx));

    let n_marked_points = 3;
    let cplx = CubicalCplx::udc(&non_hamilonian, n_marked_points);
    println!("Non_hamilonian graph #1 with {} marked points: \n{}", n_marked_points, H!(cplx));
}
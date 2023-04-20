use topo_spaces;
use topo_spaces::graph;
use topo_spaces::cubical::UdnMorseCplx;
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
    for i in 2..8 {
        let cplx = UdnMorseCplx::new(&non_hamilonian, i, Some(4));
        println!("Non_hamilonian graph #1 with {} marked points: \n{}", i, H!(cplx));
    }
}
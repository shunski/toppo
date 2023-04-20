use topo_spaces;
use topo_spaces::graph;
use topo_spaces::cubical::UdnMorseCplx;
use topo_spaces::*;

fn main() {
    // K5
    let k5 = graph! {
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

    for i in 2..6 {
        let cplx = UdnMorseCplx::new(&k5, i, Some(4));
        println!("K_5 with {} marked points: \n{}", i, H!(cplx));
    }


    // K33
    let k33 = graph! {
        ["v1", "w1"],
        ["v1", "w2"],
        ["v1", "w3"],
        ["v2", "w1"],
        ["v2", "w2"],
        ["v2", "w3"],
        ["v3", "w1"],
        ["v3", "w2"],
        ["v3", "w3"]
    };

    for i in 2..6 {
        let cplx = UdnMorseCplx::new(&k33, i, Some(4));
        println!("K_33 with {} marked points: \n{}", i, H!(cplx));
    }
}
use topo_spaces;
use topo_spaces::graph;
use topo_spaces::cubical::CubicalCplx;
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
    let n_marked_points = 2;
    let cplx = CubicalCplx::udc(&k5, n_marked_points);
    println!("K_5 with {} marked points: \n{}", n_marked_points, H!(cplx));

    let n_marked_points = 3;
    let cplx = CubicalCplx::udc(&k5, n_marked_points);
    println!("K_5 with {} marked points: \n{}", n_marked_points, H!(cplx));

    // let n_marked_points = 4;
    // let cplx = CubicalCplx::udc(&k5, n_marked_points);
    // println!("K_5 with {} marked points: \n{}", n_marked_points, H!(cplx));


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
    let n_marked_points = 2;
    let cplx = CubicalCplx::udc(&k33, n_marked_points);
    println!("K_33 with {} marked points: \n{}", n_marked_points, H!(cplx));

    let n_marked_points = 3;
    let cplx = CubicalCplx::udc(&k33, n_marked_points);
    println!("K_33 with {} marked points: \n{}", n_marked_points, H!(cplx));

    let n_marked_points = 4;
    let cplx = CubicalCplx::udc(&k33, n_marked_points);
    println!("K_33 with {} marked points: \n{}", n_marked_points, H!(cplx));
}
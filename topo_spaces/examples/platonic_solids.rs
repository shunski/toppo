use topo_spaces;
use topo_spaces::graph;
use topo_spaces::cubical::CubicalCplx;
use topo_spaces::*;

fn main() {
    // tetrahedron
    let tetrahedron = graph! {
        ["v0", "v1"],
        ["v0", "v2"],
        ["v0", "v3"],
        ["v1", "v2"],
        ["v1", "v3"],
        ["v2", "v3"]
    };
    let n_marked_points = 2;
    let cplx = CubicalCplx::udc(&tetrahedron, n_marked_points);
    println!("Tetrahedron with {} marked points: \n{}", n_marked_points, H!(cplx));

    let n_marked_points = 3;
    let cplx = CubicalCplx::udc(&tetrahedron, n_marked_points);
    println!("Tetrahedron with {} marked points: \n{}", n_marked_points, H!(cplx));

    // let n_marked_points = 4;
    // let cplx = CubicalCplx::udc(&tetrahedron, n_marked_points);
    // println!("Tetrahedron with {} marked points: \n{}", n_marked_points, H!(cplx));


    // Cube 
    let cube = graph! {
        ["v0", "v1"],
        ["v1", "v2"],
        ["v2", "v3"],
        ["v3", "v0"],
        ["v4", "v5"],
        ["v5", "v6"],
        ["v6", "v7"],
        ["v7", "v4"],
        ["v0", "v4"],
        ["v1", "v5"],
        ["v2", "v6"],
        ["v3", "v7"]
    };
    let n_marked_points = 2;
    let cplx = CubicalCplx::udc(&cube, n_marked_points);
    println!("Cube with {} marked points: \n{}", n_marked_points, H!(cplx));

    let n_marked_points = 3;
    let cplx = CubicalCplx::udc(&cube, n_marked_points);
    println!("Cube with {} marked points: \n{}", n_marked_points, H!(cplx));

    // let n_marked_points = 4;
    // let cplx = CubicalCplx::udc(&cube, n_marked_points);
    // println!("Cube with {} marked points: \n{}", n_marked_points, H!(cplx));

    // Octahedron
    let octahedron = graph! {
        ["t", "m1"],
        ["t", "m2"],
        ["t", "m3"],
        ["t", "m4"],

        ["m1", "m2"],
        ["m2", "m3"],
        ["m3", "m4"],
        ["m4", "m1"],

        ["b", "m1"],
        ["b", "m2"],
        ["b", "m3"],
        ["b", "m4"]
    };
    let n_marked_points = 2;
    let cplx = CubicalCplx::udc(&octahedron, n_marked_points);
    println!("Octahedron with {} marked points: \n{}", n_marked_points, H!(cplx));

    let n_marked_points = 3;
    let cplx = CubicalCplx::udc(&octahedron, n_marked_points);
    println!("Octahedron with {} marked points: \n{}", n_marked_points, H!(cplx));

    // Dodecahedron
    let dodecahedron = graph! {
        ["b1", "b2"],
        ["b2", "b3"],
        ["b3", "b4"],
        ["b4", "b5"],
        ["b5", "b1"],

        ["t1", "t2"],
        ["t2", "t3"],
        ["t3", "t4"],
        ["t4", "t5"],
        ["t5", "t1"],

        ["m1", "m2"],
        ["m2", "m3"],
        ["m3", "m4"],
        ["m4", "m5"],
        ["m5", "m6"],
        ["m6", "m7"],
        ["m7", "m8"],
        ["m8", "m9"],
        ["m9", "m10"],
        ["m10", "m1"],

        ["b1", "m1"],
        ["b2", "m3"],
        ["b3", "m5"],
        ["b4", "m7"],
        ["b5", "m9"],

        ["t1", "m2"],
        ["t2", "m4"],
        ["t3", "m6"],
        ["t4", "m8"],
        ["t5", "m10"]
    };
    let n_marked_points = 2;
    let cplx = CubicalCplx::udc(&dodecahedron, n_marked_points);
    println!("Dodecahedron with {} marked points: \n{}", n_marked_points, H!(cplx));

    let n_marked_points = 3;
    let cplx = CubicalCplx::udc(&dodecahedron, n_marked_points);
    println!("Dodecahedron with {} marked points: \n{}", n_marked_points, H!(cplx));
}
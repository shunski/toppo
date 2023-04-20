use topo_spaces;
use topo_spaces::graph;
use topo_spaces::cubical::UdnMorseCplx;
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
    for i in 2..7 {
        let cplx = UdnMorseCplx::new(&tetrahedron, i, Some(4));
        println!("Tetrahedron with {} marked points: \n{}", i, H!(cplx));
    }


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

    for i in 2..5 {
        let cplx = UdnMorseCplx::new(&cube, i, Some(4));
        println!("Cube with {} marked points: \n{}", i, H!(cplx));
    }


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

    for i in 2..5 {
        let cplx = UdnMorseCplx::new(&octahedron, i, Some(3));
        println!("Octahedron with {} marked points: \n{}", i, H!(cplx));
    }


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

    for i in 2..5 {
        let cplx = UdnMorseCplx::new(&dodecahedron, i, Some(3));
        println!("Dodecahedron with {} marked points: \n{}", i, H!(cplx));
    }

    // Dodecahedron
    let icosahedron = graph! {
        ["v", "x1"],
        ["v", "x2"],
        ["v", "x3"],
        ["v", "x4"],
        ["v", "x5"],

        ["x1", "x2"],
        ["x2", "x3"],
        ["x3", "x4"],
        ["x4", "x5"],
        ["x5", "x1"],
        
        ["x1", "y1"],
        ["x2", "y1"],
        ["x2", "y2"],
        ["x3", "y2"],
        ["x3", "y3"],
        ["x4", "y3"],
        ["x4", "y4"],
        ["x5", "y4"],
        ["x5", "y5"],
        ["x1", "y5"],

        ["y1", "y2"],
        ["y2", "y3"],
        ["y3", "y4"],
        ["y4", "y5"],
        ["y5", "y1"],

        ["w", "y1"],
        ["w", "y2"],
        ["w", "y3"],
        ["w", "y4"],
        ["w", "y5"]
    };

    for i in 2..5 {
        let cplx = UdnMorseCplx::new(&icosahedron, i, Some(3));
        println!("Icosaheedron with {} marked points: \n{}", i, H!(cplx));
    }
}
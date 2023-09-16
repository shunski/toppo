#[cfg(test)]
mod simplicial_cell_test {
    use crate::simplicial::SimplicialCell;

    #[test]
    fn init_test() {
        let vertices = vec![10,6,2,3];
        let mut cell = SimplicialCell::new(30, vertices);

        let mut answer = 0;
        answer |= 1 << 2;
        answer |= 1 << 3;
        answer |= 1 << 6;
        answer |= 1 << 10;
        assert_eq!(cell.vertices[0], answer);
        assert_eq!(cell.size, 30);
        assert_eq!(cell.dim, Some(3));

        cell.cone(150);
        let mut answer2 = 0;
        answer2 |= 1 << 150 % 64;
        assert_eq!(cell.vertices[0], answer);
        assert_eq!(cell.vertices[1], 0);
        assert_eq!(cell.vertices[2], answer2);
        assert_eq!(cell.size, 151);
        assert_eq!(cell.dim, Some(4));
    }

    #[test]
    fn functionality_test(){
        let vertices = vec![1,2,3,10];
        let mut cell1 = SimplicialCell::new(100, vertices);

        let vertices = vec![9, 49, 99];
        let mut cell2 = SimplicialCell::new(100, vertices);

        // testing is_connected_with()
        assert_eq!(cell1.is_connected_with(&cell2), false);
        cell1.cone(99);
        assert_eq!(cell1.is_connected_with(&cell2), true);

        // testing contains()
        assert_eq!(cell1.contains(&cell2), false);
        cell1.cone(9);
        cell1.cone(49);
        assert_eq!(cell1.contains(&cell2), true);
        cell2.cone(33);
        assert_eq!(cell1.contains(&cell2), false);
    }

    #[test]
    fn sub_cell_test() {
        let vertices = vec![9, 49, 99, 200, 3000];
        let cell1 = SimplicialCell::new(4000, vertices);
        let sub_cells = cell1.get_face();
        let ans = vec![
            SimplicialCell::new(4000, vec![49, 99, 200, 3000]),
            SimplicialCell::new(4000, vec![9, 99, 200, 3000]),
            SimplicialCell::new(4000, vec![9, 49, 200, 3000]),
            SimplicialCell::new(4000, vec![9, 49, 99, 3000]),
            SimplicialCell::new(4000, vec![9, 49, 99, 200])
        ];
        assert_eq!(sub_cells, ans);
    }
}


#[cfg(test)]
mod simplicial_cplx_tests {
    use crate::{simplex, simplicial_cplx, simplicial::*};
    #[test]
    fn object_instantiation() {
        use crate::simplicial::*;

        let mut cplx = SimplicialCplx::new();
        cplx.attach( simplex!{"v0"} );
        cplx.attach( simplex!{"v1"} );
        assert_eq!(cplx.is_connected(), false);

        let simplex = simplex!{"v1", "v2", "v3"};
        assert_eq!(simplex.dim, 2);

        cplx.attach(simplex);
        let simplex = simplex!{"v0", "v1", "v2", "v3"};
        assert_eq!(
            cplx.get_labeled_zero_skeleton().len(), 
            4
        );

        cplx.attach(simplex);
        assert_eq!(
            cplx.get_labeled_zero_skeleton().len(), 
            4
        );

        let simplex = simplex!{"v0", "v1"};
        let mut cplx = SimplicialCplx::new_from(simplex);
        assert_eq!(cplx.dim(), 1);
        cplx.attach(simplex!{ "v0", "v2", "v3", "v4" });
        assert_eq!(cplx.dim(), 3);
        assert_eq!(cplx.is_connected(), true);
    }

    #[test]
    fn n_cells_of_dim() {
        let cplx = simplicial_cplx!(
            {"a"}
        );
        assert_eq!( cplx.n_cells_of_dim(0), 1 );


        let cplx = simplicial_cplx!(
            {"a", "b"},
            {"b", "c"},
            {"c", "a"}
        );
        assert_eq!( cplx.n_cells_of_dim(0), 3 );
        assert_eq!( cplx.n_cells_of_dim(1), 3 );
        assert_eq!( cplx.n_cells_of_dim(2), 0 );


        let cplx = simplicial_cplx! {
            {"v0", "v1", "v2"},
            {"v0", "v2", "v3"},
            {"v0", "v3", "v1"},
            {"v4", "v1", "v2"},
            {"v4", "v2", "v3"},
            {"v4", "v3", "v1"}
        };
        assert_eq!( cplx.n_cells_of_dim(0), 5 );
        assert_eq!( cplx.n_cells_of_dim(1), 9 );
        assert_eq!( cplx.n_cells_of_dim(2), 6 );
        assert_eq!( cplx.n_cells_of_dim(3), 0 );


        let cplx = simplicial_cplx!(
            {"x", "y", "z", "a", "b"},
            {"x", "y", "z", "c"},
            {"x", "y", "d"},
            {"x", "y", "e"}
        );
        assert_eq!( cplx.n_cells_of_dim(0), 8 );
        assert_eq!( cplx.n_cells_of_dim(1), 10 + 3 + 2 + 2 );
        assert_eq!( cplx.n_cells_of_dim(2), 10 + 4 -1 + 1 + 1 );
        assert_eq!( cplx.n_cells_of_dim(3), 5 + 1 );
    }
}


#[cfg(test)]
mod simpilicial_cplx_space_tests {
    use crate::simplicial::*;
    use crate::*;
    use alg::lin_alg::Matrix;
    use alg::matrix;

    #[test]
    fn boundary_map() {
        // test 1
        let circle = simplicial_cplx!{
            {"v1", "v2"},
            {"v2", "v3"},
            {"v1", "v3"}
        };

        let boundary_map = circle.boundary_map();
        let b0 = Matrix::zero(1,3);
        let b1 = matrix!(i128;
            [[ 1, 0, 1],
             [-1, 1, 0],
             [ 0,-1,-1]]
        );

        assert_eq!(boundary_map[0].0, b0);
        assert_eq!(boundary_map[1].0, b1);

        // test 2
        let square = simplicial_cplx! {
            {"v0", "v1", "v2"},
            {"v0", "v1", "v3"}
        };
        let boundary_map = square.boundary_map();

        let b1 = matrix!(i128;
            [[ 1, 1, 0, 1, 0],
             [-1, 0, 1, 0, 1],
             [ 0,-1,-1, 0, 0],
             [ 0, 0, 0,-1,-1]]
        );
        let b2 = matrix!(i128;
            [[ 1, 1],
             [-1, 0],
             [ 1, 0],
             [ 0,-1],
             [ 0, 1]]
        );
        assert_eq!(boundary_map[1].0, b1);
        assert_eq!(boundary_map[2].0, b2);
        
    }

    // #[test]
    fn homology_test() {
        // test1: with the circle (S^1)
        let circle = simplicial_cplx!{
            {"v1", "v2"},
            {"v2", "v3"},
            {"v1", "v3"}
        };
        let h = H!(circle);
        assert_eq!( h.to_string(), "dim 0: Z,\ndim 1: Z,\n" );
        assert_eq!( h.cycles[0], vec![ 1*simplex!("v1"), 1*simplex!("v2"), 1*simplex!("v3") ]);
        assert_eq!( h.cycles[1], vec![ simplex!("v1", "v2") + simplex!("v2", "v3") - simplex!("v1", "v3") ]);


        // test2: with the sphere (S^2)
        let sphere = simplicial_cplx! {
            {"v0", "v1", "v2"},
            {"v0", "v2", "v3"},
            {"v0", "v3", "v1"},
            {"v4", "v1", "v2"},
            {"v4", "v2", "v3"},
            {"v4", "v3", "v1"}
        };
        assert_eq!( H!(sphere).to_string(), "dim 0: Z,\ndim 1: 0,\ndim 2: Z,\n" );


        // test3: with the torus
        let torus = simplicial_cplx!{ 
            {"v0", "v1", "v3"},
            {"v1", "v3", "v5"},
            {"v1", "v2", "v5"},
            {"v0", "v2", "v7"},
            {"v0", "v3", "v7"},
            {"v2", "v5", "v7"},
            {"v3", "v4", "v5"},
            {"v4", "v5", "v6"},
            {"v5", "v6", "v7"},
            {"v6", "v7", "v8"},
            {"v3", "v7", "v8"},
            {"v3", "v4", "v8"},
            {"v2", "v4", "v8"},
            {"v0", "v2", "v4"},
            {"v1", "v2", "v8"},
            {"v1", "v6", "v8"},
            {"v0", "v1", "v6"},
            {"v0", "v4", "v6"}
        };
        assert_eq!(H!(torus).to_string(), "dim 0: Z,\ndim 1: Z^2,\ndim 2: Z,\n");

        // test4: with the RP2
        let rp2 = simplicial_cplx!{ 
            {"1", "3", "4"},
            {"2", "3", "4"},
            {"1", "4", "5"},
            {"2", "4", "6"},
            {"4", "5", "6"},
            {"1", "2", "6"},
            {"1", "3", "6"},
            {"3", "5", "6"},
            {"2", "3", "5"},
            {"1", "2", "5"}
        };
        assert_eq!(H!(rp2).to_string(), "dim 0: Z,\ndim 1: Z_2,\n");

        // test5: with the klein bottle
        let klein_bottle = simplicial_cplx!{ 
            {"0", "1", "4"},
            {"0", "1", "6"},
            {"0", "2", "6"},
            {"0", "2", "8"},
            {"0", "3", "4"},
            {"0", "3", "8"},
            {"1", "2", "5"},
            {"1", "2", "7"},
            {"1", "4", "5"},
            {"1", "6", "7"},
            {"2", "5", "6"},
            {"2", "7", "8"},
            {"3", "4", "7"},
            {"3", "5", "6"},
            {"3", "5", "8"},
            {"3", "6", "7"},
            {"4", "5", "8"},
            {"4", "7", "8"}
        };
        assert_eq!(H!(klein_bottle).to_string(), "dim 0: Z,\ndim 1: Z + Z_2,\n");
    }
}
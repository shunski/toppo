use alg::{lin_alg::Matrix, lin_alg::SubMatrix, matrix};
use util::TupleIterable;


#[derive(Clone, Copy, PartialEq)]
struct Sphere<'a> {
    radius: f64,
    center: &'a SubMatrix<f64>,
}

impl Sphere<'_> {
    fn contains(&self, p: &SubMatrix<f64>) -> bool {
        (p - &*self.center).two_norm() < self.radius
    }
}


#[derive(Clone, Copy, PartialEq)]
pub struct Triangle<'a>([&'a SubMatrix<f64>; 3]);

impl Triangle<'_> {
    fn circum_circle_contains(&self, p: &SubMatrix<f64>) -> bool {
        let orthogonal = |mut m: Matrix<f64>| {
            let tmp = m[(0,0)];
            m[(0,0)] = -m[(1,0)];
            m[(1,0)] = tmp;
            m
        };
        let inverse = |m: &SubMatrix<f64>| {
            let det = m[(0,0)] * m[(1,1)] - m[(1,0)] * m[(0,1)];
            let inverse = matrix!(f64;
                [[ m[(1,1)], -m[(0,1)]],
                 [-m[(1,0)],  m[(0,0)]]]
            ) / det;
            debug_assert!( (m* &*inverse-Matrix::identity(2)).frobenius_norm()/2_f64.sqrt()*100.0< 0.000000001 );
            inverse
        };
        let a = self.0[0] + (self.0[1]-self.0[0]) / 2.0;
        let b = orthogonal( self.0[1]-self.0[0] );
        let c = self.0[0] + (self.0[2]-self.0[0]) / 2.0;
        let d = orthogonal( self.0[2]-self.0[0] );

        let a = &*a; let b = &*b; let c = &*c; let d = &*d; 
        
        let mut m = Matrix::zero(2, 2);
        m[(..,0)].write( b );
        m[(..,1)].write( &(d*(-1.0)) );
        let t = inverse(&m)[(0,..)].dot( &(c-a) );

        let center = a + b*t;
        let center = &*center;
        let radius = (self.0[0] - center).two_norm();

        debug_assert!({
            let s = inverse(&m)[(1,..)].dot( &(c-a) );
            (center - (c+d*s)).two_norm()/(center).two_norm()*100.0 < 0.0000000001
        });

        Sphere { radius, center }.contains(p)
    }
}

pub fn delaunay_triangulation(points: &[Matrix<f64>]) -> (Vec<Triangle<'_>>, Vec<[usize; 3]>) {
    let original_points = {
        let center = points.iter()
            .fold(Matrix::zero(2,1), |accum, x| accum+&**x )
            / points.len() as f64;
        let r = points.iter()
            .map(|p| (&**p-&*center).two_norm())
            .max_by(|p, q| p.partial_cmp(q).unwrap())
            .unwrap()
             + 1_f64;
        let p1 = &*center + matrix!(f64; [[0.0], [2.0*r]]);
        let p2 = &*center + matrix!(f64; [[ 3_f64.sqrt()*r], [-r]]);
        let p3 = &*center + matrix!(f64; [[-3_f64.sqrt()*r], [-r]]);
        [p1,p2,p3]
    };
    
    let index = |i: usize| {
        if i<points.len(){
            &points[i]
        } else {
            &original_points[i-points.len()]
        }
    };

    let original_triangle = [points.len(), points.len()+1, points.len()+2];
    let mut triangulation: Vec<[usize; 3]> = vec![original_triangle];

    for i in 0..points.len() {
        let mut bad_triangles = Vec::new();
        for &t in &triangulation {
            let triangle = Triangle([&index(t[0]), &index(t[1]), &index(t[2])]);
            if triangle.circum_circle_contains(&points[i]) {
                bad_triangles.push( t );
            }
        }
        
        let mut polygon = Vec::new();
        for &triangle in &bad_triangles {
            for edge in triangle.tuple_iter(2) {
                if bad_triangles
                    .iter()
                    .filter(|&&t| t!=triangle)
                    .map(|t| t.tuple_iter(2))
                    .flatten()
                    .all(|e| e!=edge) 
                {
                    polygon.push([*edge[0], *edge[1]]);
                }
            }
        }
        
        triangulation.retain(|x| !bad_triangles.contains(x));

        for edge in polygon {
            let mut new_t = [edge[0], edge[1], i];
            new_t.sort();
            triangulation.push( new_t );
        }
    }

    triangulation.retain(|t| t.iter().all(|&i| i<points.len()) );

    let mut out = Vec::new();
    for t in triangulation.iter() {
        assert!(t[0]!=t[1] && t[0]!=t[2] && t[1]!=t[2] );
        out.push( Triangle([&points[t[0]], &points[t[1]], &points[t[2]]]) );
    }

    (out, triangulation)
}

#[cfg(test)]
mod geom_test {
    use alg::{lin_alg::Matrix, matrix};
    use plotters::prelude::*;
    use rand::prelude::*;

    use crate::Triangle;

    #[test]
    fn circum_circle_contains() {
        let triangle = [
            matrix!(f64; [[0.0], [0.0]]),
            matrix!(f64; [[2.0], [0.0]]),
            matrix!(f64; [[0.1], [1.0]])
        ];
        let triangle = Triangle([&triangle[0], &triangle[1], &triangle[2]]);
        assert!( triangle.circum_circle_contains( &matrix!(f64; [[1.0], [0.0]]) ));
        assert!( triangle.circum_circle_contains( &matrix!(f64; [[0.5], [0.5]]) ));
        assert!(!triangle.circum_circle_contains( &matrix!(f64; [[2.2], [0.0]]) ));
        assert!(!triangle.circum_circle_contains( &matrix!(f64; [[1.0], [-1.2]]) ));
    }

    #[test]
    fn delaunay_triangulation() -> Result<(), Box<dyn std::error::Error>> {
        let n = 100;

        let mut rng = rand::thread_rng();
        let points: Vec<_> = (0..n)
            .map(|_| matrix!(f64; [[rng.gen_range(0.0..1.0)], [rng.gen_range(0.0..1.0)]]) )
            .collect();

        let (triangulation, triangle_indeces) = super::delaunay_triangulation(&points);

        for triangle in &triangulation {
            assert!( points.iter().filter(|p| triangle.circum_circle_contains(p)).count() <= 3 );
        }

        assert!(!triangulation.is_empty());
        assert!(triangle_indeces.iter().all(|t| t[0]<t[1] && t[1]<t[2] ));

        // ------------------ Drawing ------------------------- 
        let root = BitMapBackend::new("test_output/delaunay_triangulation.jpg", (600, 600)).into_drawing_area();
        root.fill(&WHITE)?;
        let mut graphic = ChartBuilder::on(&root)
            .caption("delaunay triangulation", ("san-serif", 40))
            .build_cartesian_2d(0.0..1.0, 0.0..1.0)?;

        // draw vertices
        graphic.draw_series( points.iter()
            .map(|v| Circle::new((v[(0,0)], v[(1,0)]), 4, BLACK.filled()))
        )?;

        // draw edges
        for triangle in triangulation {
            let mut triangle: Vec<_> = triangle.0.iter().map(|v| (v[(0,0)], v[(1,0)])).collect();
            triangle.push(triangle[0]);
            graphic.draw_series(LineSeries::new(
                triangle.into_iter(),
                BLACK.filled().stroke_width(2)
            ))?;
        }

        root.present()?;

        Ok(())
    }
}
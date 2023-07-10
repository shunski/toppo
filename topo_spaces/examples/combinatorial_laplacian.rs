use topo_spaces::{simplicial_cplx, simplicial::*, Complex};
use alg::lin_alg::Matrix;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
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

    let boundaries = torus.boundary_map();
    let laplacian = {
        let mut m = boundaries
            .iter()
            .map(|(x, _)| x)
            .zip( boundaries.iter().map(|(x, _)| x).skip(1) )
            .map(|(x, y)| {
                let x = &**x;
                let y = &**y;
                // compute the laplacian
                x.transpose() * x + y * y.transpose()
            } )
            .collect::<Vec<_>>();

        let highest_map = &*boundaries.last().unwrap().0;
        m.push( highest_map.transpose() * highest_map );

        let mut n = Vec::new();
        for map in m {
            let mut real_map = Matrix::zero( map.size().0, map.size().1 );
            for (i, j) in (0..map.size().0).map(|i| (0..map.size().1).map(move |j| (i,j))).flatten() {
                real_map[(i,j)] = map[(i,j)] as f64;
            }
            n.push( real_map );
        }
        n
    };


    let eigenvals = laplacian
        .into_iter()
        .map( |map| map.spectrum() )
        .collect::<Vec<_>>();

    draw( &eigenvals, "combinatorial_laplacian.png" )?;

    Ok(())
}


// visualization using plotter-rs
fn draw( eigenvals: &Vec<Matrix<f64>>, file_name: &str ) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(file_name, (1500, 600)).into_drawing_area();

    root.fill(&WHITE)?;
    let root = root.margin(5, 5, 5, 5);
    let areas = root.split_evenly((1, eigenvals.len()));

    //  draw the histogram on each of the panels
    for (i, (area, data)) in areas.iter().zip( eigenvals.iter() ).enumerate() {
        let mut histogram = ChartBuilder::on(&area)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .margin(2)
            .caption(format!("Spectrum of Dim {}", i), ("sans-serif", 30.0))
            .build_cartesian_2d((0..data.size().0).into_segmented(), 0f64..12f64)?;

        histogram
            .configure_mesh()
            .disable_x_mesh()
            .disable_y_mesh()
            .y_desc("")
            .x_desc("")
            .draw()?;

        histogram.draw_series(
            Histogram::vertical(&histogram)
                .style(RED.filled())
                .data((0..data.size().0).map( |i| (i, data[(i,0)]) )),
        )?;
    }

    root.present().expect("");

    Ok(())
}


use topo_spaces::{simplicial_cplx, simplicial::*, Complex};
use alg::lin_alg::{Matrix, SubMatrix};
use plotters::prelude::*;

static FILE_NAME: &str = "output.jpg";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // create the base space
    let cplx = simplicial_cplx!{ 
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

    // compute boundary maps
    let boundaries = cplx.boundary_map();

    // a little closure to convert integral matrices to real matrices
    let to_real = { |m: &SubMatrix<i64>| {
            let mut real_map = Matrix::zero( m.size().0, m.size().1 );
            for (i, j) in (0..m.size().0).map(|i| (0..m.size().1).map(move |j| (i,j))).flatten() {
                real_map[(i,j)] = m[(i,j)] as f64;
            }
            real_map
        }
    };

    // n is the number of iteration over weights
    let n = 20;

    // spectrum that is plotted on the graph
    // first axis: dimension
    // second dim: distinct eigenvalues (smallest to largest)
    // third dim: eigenvalue with respect to weights
    let mut output: Vec<Vec<Vec<f64>>> = (0..boundaries.len()).map(|i|
        vec![Vec::new(); cplx.n_cells_of_dim(i)]
    ).collect();

    let weights: Vec<_> = (0..=cplx.dim()).map( |i| {
        (1..=n).map(|t| {
                let mut m = Matrix::identity( cplx.n_cells_of_dim(i) );
                if i==1 {
                    m[(0,0)] /= t as f64;
                    m[(1,1)] /= t as f64;
                    m[(2,2)] /= t as f64;
                } else if i==2 {
                    m[(0,0)] /= t as f64;
                }
                m
            })
            .collect::<Vec<_>>()
    }).collect();

    // iterate over dimensions of the chain complex
    for i in 0..boundaries.len() {
        // boundary maps
        let b_low = to_real(&boundaries[i].0);
        let b_low = &*b_low;
        let b_high = if i<boundaries.len()-1 {
            to_real(&boundaries[i+1].0)
        } else {
            Matrix::zero(cplx.n_cells_of_dim(cplx.dim()), 1)
        };
        let b_high = &*b_high;

        // iterate over the scences and collect the output
        for t in 0..n {
            // weights
            let w_low = if i>0 { 
                weights[i-1][t].clone()
            } else {
                Matrix::identity(1)
            };
            let w_low = &*w_low;

            let w = weights[i][t].clone();
            let w_inv = {
                let mut m = w.clone();
                (0..m.size().0).for_each(|i| m[(i,i)]=1.0/m[(i,i)] );
                m
            };
            let w = &*w;
            let w_inv = &*w_inv;

            let w_high_inv = if i < boundaries.len()-1 {
                let mut m = weights[i+1][t].clone();
                (0..m.size().0).for_each(|i| m[(i,i)]=1.0/m[(i,i)] );
                m
            } else {
                Matrix::<f64>::identity(1)
            };
            let w_high_inv = &*w_high_inv;

            // computation
            let l_down = w_inv * b_low.transpose() * w_low * b_low;
            let l_up = b_high * w_high_inv * b_high.transpose() * w;
            let weighted_laplacian = l_down + l_up;
            let eigenvals = weighted_laplacian.spectrum();

            for j in 0..eigenvals.size().0 {
                output[i][j].push(eigenvals[(j,0)]);
            }
        };
    }


    // ------------ drawing ---------------
    let root = BitMapBackend::new(FILE_NAME, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.margin(5, 5, 5, 5);
    let n_cols = ((cplx.dim()+1) as f64 - 0.1).sqrt() as usize + 1;
    let n_rows = cplx.dim()/n_cols + 1;
    let panels = root.split_evenly((n_rows, n_cols));

    // iterate over the dimensions (or the panels)
    for (i, panel) in panels.iter().enumerate() {
        if i >= output.len() {
            break;
        }
        // create the chart
        let mut chart = ChartBuilder::on(&panel)
            .set_label_area_size(LabelAreaPosition::Left, 60)
            .set_label_area_size(LabelAreaPosition::Bottom, 60)
            .caption(format!( "Spectrum of L_{}", i ), ("sans-serif", 20))
            .build_cartesian_2d(0..n, -1.0..30.0)?;

        // configure mesh
        chart
            .configure_mesh()
            .disable_x_mesh()
            .disable_y_mesh()
            .draw()?;

        // plot the data
        for eigenvalue in output[i].iter() {
            chart.draw_series(
                LineSeries::new(
                    (0..).zip(eigenvalue.iter()).map(|(x, y)| (x, *y)),
                    &RED.mix(0.5),
                )
            )?;
        }

    }

    root.present()?;

    Ok(())
}
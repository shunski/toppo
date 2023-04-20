use std::collections::HashMap;
use topo_spaces::graph::WeightedSimpleGraph;
use topo_spaces::graph;
use alg::lin_alg::Matrix;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let embedding: HashMap<&str, (f64, f64)> = HashMap::from([
        ("v1", (2.4, 1.3)),
        ("v2", (2.2, 0.8)),
        ("v3", (1.9, 0.7)),
        ("v4", (1.8, 1.1)),

        ("v5", (-0.8, 2.0)),
        ("v6", (-1.0, 1.6)),
        ("v7", (-1.2, 2.4)),
        ("v8", (-1.4, 1.7)),
        ("v9", (-0.7, 1.9)),

        ("v10", (-0.9, -2.0)),
        ("v11", (-1.1, -1.7)),
        ("v12", (-0.6, -2.4)),
        ("v13", (-0.8, -1.5)),
        ("v14", (-1.3, -2.3)),

        ("v15", ( 0.1, -0.0)),
        ("v16", (-0.1, -0.3)),
        ("v17", ( 0.3,  0.3)),
        ("v18", (-0.4, -0.4)),
        ("v19", ( 0.3,  0.0)),
        ("v20", ( 0.0,  0.2)),
        ("v21", ( 0.3, -0.2)),
    ]);

    let distance = |x:(f64,f64), y:(f64,f64)| {
        ((y.0-x.0)*(y.0-x.0) + (y.1-x.1)*(y.1-x.1)).sqrt()
    };

    let graph = graph! { 
        embedding=embedding, 
        weight_fn=distance,
        [ "v1",  "v2"],
        [ "v1",  "v3"],
        [ "v1",  "v4"],
        [ "v2",  "v3"],
        [ "v2",  "v4"],
        [ "v3",  "v4"],

        [ "v4",  "v9"],
        [ "v5",  "v6"],
        [ "v5",  "v7"],
        [ "v5",  "v8"],
        [ "v6",  "v7"],
        [ "v6",  "v8"],
        [ "v6",  "v9"],
        [ "v7",  "v8"],
        [ "v7",  "v9"],

        [ "v3", "v12"],
        ["v10", "v11"],
        ["v10", "v12"],
        ["v10", "v13"],
        ["v10", "v14"],
        ["v11", "v12"],
        ["v11", "v13"],
        ["v11", "v14"],
        ["v12", "v14"],
        ["v12", "v13"],
        ["v13", "v14"],

        ["v15", "v16"],
        ["v15", "v17"],
        ["v15", "v18"],
        ["v15", "v20"],
        ["v16", "v17"],
        ["v16", "v19"],
        ["v16", "v21"],
        ["v17", "v18"],
        ["v17", "v19"],
        ["v18", "v19"],
        ["v18", "v21"],
        ["v20", "v21"],
        [ "v6", "v19"],
        ["v13", "v18"]
    };

    // get the laplacian and compute the eigenvalue
    let l = graph.clone().laplacian();
    let (eigenvals, eigenvecs) = l.spectrum_with_invariant_space();

    // visualization
    draw(&graph, embedding, &eigenvals, &eigenvecs, "graph_laplacian.png")?;

    Ok(())
}


// visualization using plotter-rs
fn draw( graph: &WeightedSimpleGraph<f64>, embedding: HashMap<&str, (f64, f64)>, eigenvals: &Matrix<f64>, eigenvec: &Matrix<f64>, file_name: &str ) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(file_name, (1500, 600)).into_drawing_area();

    root.fill(&WHITE)?;
    let root = root.margin(5, 5, 5, 5);
    let (left, right) = root.split_horizontally(30.percent_width());

    //  draw the histogram on the left panel
    let mut histogram = ChartBuilder::on(&left)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .margin(2)
        .caption("Spectrum", ("sans-serif", 30.0))
        .build_cartesian_2d((0..eigenvals.size.0).into_segmented(), 0f64..12f64)?;

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
            .data((0..eigenvals.size.0).map( |i| (i, eigenvals[(i,0)]) )),
    )?;

    
    // draw the colored graphs on the right
    // 'offset' is the dimention of the null space.
    let offset = (0..).find(|&i| eigenvals[(i,0)].abs() > 0.00001).unwrap();
    let areas_for_graphs = right.split_evenly((2, 4));
    // iterating over the 8 panels in the right...
    for (k, area) in areas_for_graphs.iter().enumerate() {
        let mut chart = ChartBuilder::on( area )
            .caption(format!("u{}", 2*k+offset), ("sans-serif", 20))
            .margin(20)
            .build_cartesian_2d(-2.5f64..2.5f64, -2.5f64..2.5f64)?;

        chart.draw_series(
            graph.vertices().iter()
                .enumerate()
                .map(|(i, v)| Circle::new(
                    *embedding.get(&v[..] ).unwrap(), 
                    6,
                    HSLColor( 1.8/3.6 + (eigenvec[(i,2*k+offset)] + 1.0) / 3.0, 0.8, 0.5).filled()
                ) ),
        )?;

        let n = graph.vertices().len();
        let edge_iter = (0..n)
            .map(|i| (0..i).map(move |j| (i,j)))
            .flatten()
            .filter(|&(i,j)| graph.data()[(i,j)].abs()>0.000001 )
            .map(|(i,j)| ( &graph.vertices()[i][..], &graph.vertices()[j][..]) );

        for (v1, v2) in edge_iter {
            chart.draw_series(LineSeries::new(
                (0..2).map(|i| *embedding.get( if i==0 { v1 } else { v2 }).unwrap()),
                &BLACK,
            ))?;
        }
    }

    root.present().expect("");

    Ok(())
}

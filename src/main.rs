mod clustering;
mod detector_module;
mod event_generator;
mod file_reader;
mod file_writer;
mod matrix;
use plotters::prelude::*;

const DETECTOR_FILE_NAME: &str = "data/boxDetector.txt";
// const HITS_FILE_NAME: &str = "data/someHits.txt";
const OUT_FILE_NAME: &str = "data/boxDetectorRandomRays.png";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let layers = vec![5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.5, 15.0, 20.0, 25.0, 35.0];
    let plot_dim = 40.0;
    file_writer::create_box_detector(&layers, DETECTOR_FILE_NAME)?;

    let modules = file_reader::read_modules(DETECTOR_FILE_NAME)?;
    // let hits = file_reader::read_hits(HITS_FILE_NAME)?;

    let root = BitMapBackend::new(OUT_FILE_NAME, (1920, 1080)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("boxDetector", ("sans-serif", 40))
        .build_cartesian_3d(
            -plot_dim / 2.0..plot_dim / 2.0,
            -plot_dim / 2.0..plot_dim / 2.0,
            -plot_dim / 2.0..plot_dim / 2.0,
        )?;

    // chart.draw_series(std::iter::once(PathElement::new(
    //     vec![(0., 0., 0.), (2., 0., 0.)],
    //     RED,
    // )))?; // X axis
    // chart.draw_series(std::iter::once(PathElement::new(
    //     vec![(0., 0., 0.), (0., 2., 0.)],
    //     BLUE,
    // )))?; // Y axis
    // chart.draw_series(std::iter::once(PathElement::new(
    //     vec![(0., 0., 0.), (0., 0., 2.)],
    //     YELLOW,
    // )))?; // Z axis

    for (_, module) in &modules {
        #[allow(unused_mut)]
        let mut verts: Vec<_> = module
            .vertices()
            .iter()
            .map(|&it| (it.data[0], it.data[1], it.data[2]))
            .collect();
        // verts.push(verts[0]);

        chart.draw_series(std::iter::once(PathElement::new(verts, RED)))?;
    }

    // for hit in &hits {
    //     let mut verts: Vec<_> = (&modules)[&hit.module_id]
    //         .pixel_vertices(hit.pos)?
    //         .iter()
    //         .map(|&it| (it.data[0], it.data[1], it.data[2]))
    //         .collect();
    //     verts.push(verts[0]);
    //     chart.draw_series(std::iter::once(PathElement::new(verts, BLUE)))?;
    // }

    // let points = clustering::clustering(&modules, hits)?;

    // for point in points {
    //     chart.plotting_area().draw(&plotters::element::Circle::new(
    //         (point.data[0], point.data[1], point.data[2]),
    //         5,
    //         BLACK.filled(),
    //     ))?;
    // }

    let dist = event_generator::Distributions::Gauss(event_generator::Gauss::new(50., 40.));
    let rays = event_generator::generate_random_rays(500., 5., dist);
    for ray in &rays {
        chart.draw_series(std::iter::once(PathElement::new(
            [
                (0., 0., 0.),
                ray.at_radius(
                    (plot_dim / 2f64)
                        .min(ray.energy / event_generator::ENERGY_LOSS_PER_UNIT_DISTANCE),
                )?,
            ],
            BLACK,
        )))?;
    }

    for ray in &rays {
        for (_, module) in &modules {
            if let Some((pos, _)) = ray.intersect(module) {
                chart.plotting_area().draw(&plotters::element::Circle::new(
                    pos,
                    3,
                    GREEN.filled(),
                ))?;
            }
        }
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Saved output to file: {}", OUT_FILE_NAME);

    Ok(())
}

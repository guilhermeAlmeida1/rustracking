mod detector_module;
mod file_reader;
mod file_writer;
mod matrix;
use plotters::prelude::*;

const DETECTOR_FILE_NAME: &str = "data/boxDetector.txt";
const OUT_FILE_NAME: &str = "data/boxDetector.png";
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let layers = vec![5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.5, 15.0, 20.0, 25.0, 35.0];
    let plot_dim = 40.0;
    file_writer::create_box_detector(&layers, DETECTOR_FILE_NAME)?;

    let modules = file_reader::read_modules(DETECTOR_FILE_NAME)?;

    let root = BitMapBackend::new(OUT_FILE_NAME, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("boxDetector", ("sans-serif", 40))
        .build_cartesian_3d(
            -plot_dim / 2.0..plot_dim / 2.0,
            -plot_dim / 2.0..plot_dim / 2.0,
            -plot_dim / 2.0..plot_dim / 2.0,
        )?;

    for module in modules {
        let projection: Vec<_> = module
            .vertices()
            .iter()
            .map(|&it| (it.data[0], it.data[1], it.data[2]))
            .collect();

        chart.draw_series(std::iter::once(PathElement::new(projection, RED)))?;
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Saved output to file: {}", OUT_FILE_NAME);

    Ok(())
}

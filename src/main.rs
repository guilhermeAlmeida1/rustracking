mod detector_module;
mod file_reader;
mod matrix;
mod plotter;
use detector_module::*;
use matrix::*;
use plotters::prelude::*;

const OUT_FILE_NAME: &str = "plot.png";
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let modules = file_reader::read_modules("data/boxDetector.txt")?;

    let root = BitMapBackend::new(OUT_FILE_NAME, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("caption", ("sans-serif", 50))
        .build_cartesian_2d(-1.0..6.0, -1.0..6.0)?;

    for module in modules {
        // let projection: Vec<_> = module
        //     .vertices()
        //     .iter()
        //     .map(|&it| (it.data[0], it.data[1]))
        //     .collect();
        let projection: Vec<_> = module.vertices()
            .iter()
            .map(|&it| (it.data[1], it.data[2]))
            .collect();
        println!("{:?}", module.vertices());

        chart.draw_series(std::iter::once(PathElement::new(projection, RED)))?;
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");

    Ok(())
}

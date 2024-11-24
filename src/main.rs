#![allow(dead_code)]

mod clustering;
mod detector_module;
mod event_generator;
mod field;
mod file_reader;
mod file_writer;
mod matrix;
mod particle;
mod plotting;
mod utils;

const DETECTOR_FILE_NAME: &str = "data/boxDetector.txt";
const PARTICLES_FILE_NAME: &str = "data/magfield/particles.txt";
const HITS_FILE_NAME: &str = "data/magfield/hits.txt";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // let layers = vec![2.5, 3., 3.5, 4., 4.5, 5., 6.25, 7.5, 10., 12.5, 17.5];
    let plot_dim = 40.0;
    let plot_dims: plotting::Boundary = (
        -plot_dim / 2.0..plot_dim / 2.0,
        -plot_dim / 2.0..plot_dim / 2.0,
        -plot_dim / 2.0..plot_dim / 2.0,
    );

    // file_writer::create_box_detector(&layers, DETECTOR_FILE_NAME)?;

    let modules = file_reader::read_modules(DETECTOR_FILE_NAME)?;

    // let dist = event_generator::Distributions::Gauss(event_generator::Gauss::new(5000., 2500.));
    // let particles = event_generator::generate_random_particles(30000., crate::particle::PROTON_REST_MASS, dist)?;
    // let hits = event_generator::create_hits(particles.clone(), &modules, &plot_dims);

    // file_writer::write_particles(PARTICLES_FILE_NAME, &particles)?;
    // file_writer::write_hits(HITS_FILE_NAME, &hits)?;

    let particles = file_reader::read_particles(PARTICLES_FILE_NAME)?;
    let hits = file_reader::read_hits(HITS_FILE_NAME)?;

    let points = clustering::clustering(&modules, &hits)?;

    plotting::plot_all(
        "data/magfield",
        &plot_dims,
        &modules,
        None,
        Some(&points),
        Some(particles.clone()),
    )?;

    Ok(())
}

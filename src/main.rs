mod clustering;
mod detector_module;
mod event_generator;
mod file_reader;
mod file_writer;
mod matrix;
mod plotting;

const DETECTOR_FILE_NAME: &str = "data/boxDetector.txt";
const RAYS_FILE_NAME: &str = "data/randomHits/rays.txt";
const HITS_FILE_NAME: &str = "data/randomHits/hits.txt";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let layers = vec![2.5, 3., 3.5, 4., 4.5, 5., 6.25, 7.5, 10., 12.5, 17.5];
    let plot_dim = 40.0;
    let plot_dims = (
        -plot_dim / 2.0..plot_dim / 2.0,
        -plot_dim / 2.0..plot_dim / 2.0,
        -plot_dim / 2.0..plot_dim / 2.0,
    );

    file_writer::create_box_detector(&layers, DETECTOR_FILE_NAME)?;

    let modules = file_reader::read_modules(DETECTOR_FILE_NAME)?;

    let dist = event_generator::Distributions::Gauss(event_generator::Gauss::new(10., 10.));
    let (rays, hits) = event_generator::generate_random_event(70., 2.5, dist, &modules);

    file_writer::write_rays(RAYS_FILE_NAME, &rays)?;
    file_writer::write_hits(HITS_FILE_NAME, &hits)?;

    let rays = file_reader::read_rays(RAYS_FILE_NAME)?;
    let hits = file_reader::read_hits(HITS_FILE_NAME)?;

    // let rays = vec![event_generator::StraightRay {
    //     theta: 0.,
    //     phi: 0.,
    //     origin: (0., 0., 0.),
    //     energy: 15.,
    // }];
    // let hits = event_generator::create_hits(&rays, &modules);

    let points = clustering::clustering(&modules, &hits)?;

    plotting::plot_3d(
        "data/boxDetector.svg",
        &plot_dims,
        &modules,
        None,
        false,
        None,
        None,
    )?;

    plotting::plot_all(
        "data/randomHits",
        &plot_dims,
        &modules,
        Some(&rays),
        true,
        Some(&hits),
        // None,
        Some(&points),
    )?;

    Ok(())
}

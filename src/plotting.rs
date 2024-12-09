use crate::clustering::{Hit, SpacePoint};
use crate::detector_module::DetectorModule;
use plotters::prelude::*;

use crate::particle::*;

use std::collections::HashMap;

pub type Boundary = (
    std::ops::Range<f64>,
    std::ops::Range<f64>,
    std::ops::Range<f64>,
);

pub fn plot_3d(
    filename: &str,
    plot_dims: &Boundary,
    modules: &HashMap<u64, DetectorModule>,
    hits: Option<&Vec<Hit>>,
    clustered_points: Option<&Vec<SpacePoint>>,
    particles: Option<Vec<Particle>>,
    particles_as_points: Option<&Vec<Vec<(f64, f64, f64)>>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new(filename, (1920, 1080)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("boxDetector", ("sans-serif", 40))
        .build_cartesian_3d(
            plot_dims.0.clone(),
            plot_dims.1.clone(),
            plot_dims.2.clone(),
        )?;

    for (_, module) in modules {
        #[allow(unused_mut)]
        let mut verts: Vec<_> = module
            .vertices()
            .iter()
            .map(|&it| (it.data[0], it.data[1], it.data[2]))
            .collect();
        // verts.push(verts[0]);

        chart.draw_series(std::iter::once(PathElement::new(verts, RED)))?;
    }

    if let Some(hits) = hits {
        for hit in hits {
            let mut verts: Vec<_> = (&modules)[&hit.module_id]
                .pixel_vertices(hit.pos)?
                .iter()
                .map(|&it| (it.data[0], it.data[1], it.data[2]))
                .collect();
            verts.push(verts[0]);
            chart.draw_series(std::iter::once(PathElement::new(verts, BLUE)))?;
        }
    }

    if let Some(points) = clustered_points {
        for point in points {
            chart.plotting_area().draw(&plotters::element::Circle::new(
                (point.0[0], point.0[1], point.0[2]),
                3,
                BLACK.filled(),
            ))?;
        }
    }

    if let Some(mut particles) = particles {
        for particle in &mut particles {
            let mut points = vec![(
                particle.ray.origin[0],
                particle.ray.origin[1],
                particle.ray.origin[2],
            )];
            for _ in 0.. {
                if particle.energy - particle.rest_mass < 0. {
                    break;
                }
                particle.do_step();
                points.push((
                    particle.ray.origin[0],
                    particle.ray.origin[1],
                    particle.ray.origin[2],
                ));
            }
            chart
                .draw_series(std::iter::once(PathElement::new(points, BLUE)))
                .unwrap();
        }
    }

    if let Some(particles_as_points) = particles_as_points {
        for points_vec in particles_as_points {
            chart
                .draw_series(std::iter::once(PathElement::new(points_vec.clone(), BLUE)))
                .unwrap();
        }
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Saved output to file: {}", filename);

    Ok(())
}

pub fn plot_2d_xy(
    filename: &str,
    plot_dims: &Boundary,
    modules: &HashMap<u64, DetectorModule>,
    hits: Option<&Vec<Hit>>,
    clustered_points: Option<&Vec<SpacePoint>>,
    particles: Option<Vec<Particle>>,
    particles_as_points: Option<&Vec<Vec<(f64, f64, f64)>>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new(filename, (1920, 1080)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("boxDetector", ("sans-serif", 40))
        .build_cartesian_2d(plot_dims.0.clone(), plot_dims.1.clone())?;

    for (_, module) in modules {
        #[allow(unused_mut)]
        let mut verts: Vec<_> = module
            .vertices()
            .iter()
            .map(|&it| (it.data[0], it.data[1]))
            .collect();
        // verts.push(verts[0]);

        chart.draw_series(std::iter::once(PathElement::new(verts, RED)))?;
    }

    if let Some(hits) = hits {
        for hit in hits {
            let mut verts: Vec<_> = (&modules)[&hit.module_id]
                .pixel_vertices(hit.pos)?
                .iter()
                .map(|&it| (it.data[0], it.data[1]))
                .collect();
            verts.push(verts[0]);
            chart.draw_series(std::iter::once(PathElement::new(verts, BLUE)))?;
        }
    }

    if let Some(points) = clustered_points {
        for point in points {
            chart.plotting_area().draw(&plotters::element::Circle::new(
                (point.0[0], point.0[1]),
                3,
                BLACK.filled(),
            ))?;
        }
    }

    if let Some(mut particles) = particles {
        for particle in &mut particles {
            let mut points = vec![(particle.ray.origin[0], particle.ray.origin[1])];
            for _ in 0.. {
                if particle.energy - particle.rest_mass < 0. {
                    break;
                }
                particle.do_step();
                points.push((particle.ray.origin[0], particle.ray.origin[1]));
            }
            chart
                .draw_series(std::iter::once(PathElement::new(points, BLUE)))
                .unwrap();
        }
    }

    if let Some(particles_as_points) = particles_as_points {
        for points_vec in particles_as_points {
            chart
                .draw_series(std::iter::once(PathElement::new(
                    points_vec
                        .into_iter()
                        .map(|itr| (itr.0, itr.1))
                        .collect::<Vec<_>>(),
                    BLUE,
                )))
                .unwrap();
        }
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Saved output to file: {}", filename);

    Ok(())
}

pub fn plot_2d_xz(
    filename: &str,
    plot_dims: &Boundary,
    modules: &HashMap<u64, DetectorModule>,
    hits: Option<&Vec<Hit>>,
    clustered_points: Option<&Vec<SpacePoint>>,
    particles: Option<Vec<Particle>>,
    particles_as_points: Option<&Vec<Vec<(f64, f64, f64)>>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new(filename, (1920, 1080)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("boxDetector", ("sans-serif", 40))
        .build_cartesian_2d(plot_dims.0.clone(), plot_dims.2.clone())?;

    for (_, module) in modules {
        #[allow(unused_mut)]
        let mut verts: Vec<_> = module
            .vertices()
            .iter()
            .map(|&it| (it.data[0], it.data[2]))
            .collect();
        // verts.push(verts[0]);

        chart.draw_series(std::iter::once(PathElement::new(verts, RED)))?;
    }

    if let Some(hits) = hits {
        for hit in hits {
            let mut verts: Vec<_> = (&modules)[&hit.module_id]
                .pixel_vertices(hit.pos)?
                .iter()
                .map(|&it| (it.data[0], it.data[2]))
                .collect();
            verts.push(verts[0]);
            chart.draw_series(std::iter::once(PathElement::new(verts, BLUE)))?;
        }
    }

    if let Some(points) = clustered_points {
        for point in points {
            chart.plotting_area().draw(&plotters::element::Circle::new(
                (point.0[0], point.0[2]),
                3,
                BLACK.filled(),
            ))?;
        }
    }

    if let Some(mut particles) = particles {
        for particle in &mut particles {
            let mut points = vec![(particle.ray.origin[0], particle.ray.origin[2])];
            for _ in 0.. {
                if particle.energy - particle.rest_mass < 0. {
                    break;
                }
                particle.do_step();
                points.push((particle.ray.origin[0], particle.ray.origin[2]));
            }
            chart
                .draw_series(std::iter::once(PathElement::new(points, BLUE)))
                .unwrap();
        }
    }

    if let Some(particles_as_points) = particles_as_points {
        for points_vec in particles_as_points {
            chart
                .draw_series(std::iter::once(PathElement::new(
                    points_vec
                        .into_iter()
                        .map(|itr| (itr.0, itr.2))
                        .collect::<Vec<_>>(),
                    BLUE,
                )))
                .unwrap();
        }
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Saved output to file: {}", filename);

    Ok(())
}

pub fn plot_2d_yz(
    filename: &str,
    plot_dims: &Boundary,
    modules: &HashMap<u64, DetectorModule>,
    hits: Option<&Vec<Hit>>,
    clustered_points: Option<&Vec<SpacePoint>>,
    particles: Option<Vec<Particle>>,
    particles_as_points: Option<&Vec<Vec<(f64, f64, f64)>>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new(filename, (1920, 1080)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("boxDetector", ("sans-serif", 40))
        .build_cartesian_2d(plot_dims.1.clone(), plot_dims.2.clone())?;

    for (_, module) in modules {
        #[allow(unused_mut)]
        let mut verts: Vec<_> = module
            .vertices()
            .iter()
            .map(|&it| (it.data[1], it.data[2]))
            .collect();
        // verts.push(verts[0]);

        chart.draw_series(std::iter::once(PathElement::new(verts, RED)))?;
    }

    if let Some(hits) = hits {
        for hit in hits {
            let mut verts: Vec<_> = (&modules)[&hit.module_id]
                .pixel_vertices(hit.pos)?
                .iter()
                .map(|&it| (it.data[1], it.data[2]))
                .collect();
            verts.push(verts[0]);
            chart.draw_series(std::iter::once(PathElement::new(verts, BLUE)))?;
        }
    }

    if let Some(points) = clustered_points {
        for point in points {
            chart.plotting_area().draw(&plotters::element::Circle::new(
                (point.0[1], point.0[2]),
                3,
                BLACK.filled(),
            ))?;
        }
    }

    if let Some(mut particles) = particles {
        for particle in &mut particles {
            let mut points = vec![(particle.ray.origin[1], particle.ray.origin[2])];
            for _ in 0.. {
                if particle.energy - particle.rest_mass < 0. {
                    break;
                }
                particle.do_step();
                points.push((particle.ray.origin[1], particle.ray.origin[2]));
            }
            chart
                .draw_series(std::iter::once(PathElement::new(points, BLUE)))
                .unwrap();
        }
    }

    if let Some(particles_as_points) = particles_as_points {
        for points_vec in particles_as_points {
            chart
                .draw_series(std::iter::once(PathElement::new(
                    points_vec
                        .into_iter()
                        .map(|itr| (itr.1, itr.2))
                        .collect::<Vec<_>>(),
                    BLUE,
                )))
                .unwrap();
        }
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Saved output to file: {}", filename);

    Ok(())
}

pub fn plot_all(
    file_dir: &str,
    plot_dims: &Boundary,
    modules: &HashMap<u64, DetectorModule>,
    hits: Option<&Vec<Hit>>,
    clustered_points: Option<&Vec<SpacePoint>>,
    particles: Option<Vec<Particle>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let filename_3d = file_dir.to_string() + "/3d.svg";
    let filename_2d_xy = file_dir.to_string() + "/2d_xy.svg";
    let filename_2d_xz = file_dir.to_string() + "/2d_xz.svg";
    let filename_2d_yz = file_dir.to_string() + "/2d_yz.svg";

    let mut particles_as_points = Vec::new();

    if let Some(mut particles) = particles {
        for particle in &mut particles {
            let mut points = vec![(
                particle.ray.origin[0],
                particle.ray.origin[1],
                particle.ray.origin[2],
            )];
            for _ in 0.. {
                if particle.energy - particle.rest_mass < 0. {
                    break;
                }
                particle.do_step();
                points.push((
                    particle.ray.origin[0],
                    particle.ray.origin[1],
                    particle.ray.origin[2],
                ));
            }
            particles_as_points.push(points);
        }
    }

    plot_3d(
        &filename_3d,
        plot_dims,
        modules,
        hits,
        clustered_points,
        None,
        Some(&particles_as_points),
    )?;
    plot_2d_xy(
        &filename_2d_xy,
        plot_dims,
        modules,
        hits,
        clustered_points,
        None,
        Some(&particles_as_points),
    )?;
    plot_2d_xz(
        &filename_2d_xz,
        plot_dims,
        modules,
        hits,
        clustered_points,
        None,
        Some(&particles_as_points),
    )?;
    plot_2d_yz(
        &filename_2d_yz,
        plot_dims,
        modules,
        hits,
        clustered_points,
        None,
        Some(&particles_as_points),
    )?;

    Ok(())
}

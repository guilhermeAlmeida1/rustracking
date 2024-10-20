use crate::clustering::Hit;
use crate::detector_module::DetectorModule;
use crate::event_generator::StraightRay;
use crate::matrix::Vector3;
use plotters::prelude::*;

use std::collections::HashMap;

type PlotDims = (
    std::ops::Range<f64>,
    std::ops::Range<f64>,
    std::ops::Range<f64>,
);

fn at_max_radius(
    ray: &StraightRay,
    plot_dims: &PlotDims,
) -> Result<(f64, f64, f64), Box<dyn std::error::Error>> {
    let mut result = ray.end();
    if result.0 < plot_dims.0.start {
        result =
            ray.at_radius((plot_dims.0.start / ray.theta.sin() / ray.phi.cos()).abs() - 0.01)?;
    }
    if result.0 > plot_dims.0.end {
        result = ray.at_radius((plot_dims.0.end / ray.theta.sin() / ray.phi.cos()).abs() - 0.01)?;
    }
    if result.1 < plot_dims.1.start {
        result =
            ray.at_radius((plot_dims.1.start / ray.theta.sin() / ray.phi.sin()).abs() - 0.01)?;
    }
    if result.1 > plot_dims.1.end {
        result = ray.at_radius((plot_dims.1.end / ray.theta.sin() / ray.phi.sin()).abs() - 0.01)?;
    }
    if result.2 < plot_dims.2.start {
        result = ray.at_radius(plot_dims.2.start / ray.theta.cos().abs() - 0.01)?;
    }
    if result.2 > plot_dims.2.end {
        result = ray.at_radius(plot_dims.2.end / ray.theta.cos().abs() - 0.01)?;
    }

    Ok(result)
}

pub fn plot_3d(
    filename: &str,
    plot_dims: &PlotDims,
    modules: &HashMap<u64, DetectorModule>,
    rays: Option<&Vec<StraightRay>>,
    plot_real_intersections: bool,
    hits: Option<&Vec<Hit>>,
    clustered_points: Option<&Vec<Vector3<f64>>>,
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

    if let Some(rays) = rays {
        for ray in rays {
            chart.draw_series(std::iter::once(PathElement::new(
                [(0., 0., 0.), at_max_radius(ray, &plot_dims)?],
                BLACK,
            )))?;
        }

        if plot_real_intersections {
            for ray in rays {
                for (_, module) in modules {
                    if let Some(((pos, _), _)) = ray.intersect(module) {
                        chart.plotting_area().draw(&plotters::element::Circle::new(
                            pos,
                            3,
                            GREEN.filled(),
                        ))?;
                    }
                }
            }
        }
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
                (point.data[0], point.data[1], point.data[2]),
                3,
                BLACK.filled(),
            ))?;
        }
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Saved output to file: {}", filename);

    Ok(())
}

pub fn plot_2d_xy(
    filename: &str,
    plot_dims: &PlotDims,
    modules: &HashMap<u64, DetectorModule>,
    rays: Option<&Vec<StraightRay>>,
    plot_real_intersections: bool,
    hits: Option<&Vec<Hit>>,
    clustered_points: Option<&Vec<Vector3<f64>>>,
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

    if let Some(rays) = rays {
        for ray in rays {
            chart.draw_series(std::iter::once(PathElement::new(
                [
                    (0., 0.),
                    (
                        at_max_radius(ray, &plot_dims)?.0,
                        at_max_radius(ray, &plot_dims)?.1,
                    ),
                ],
                BLACK,
            )))?;
        }

        if plot_real_intersections {
            for ray in rays {
                for (_, module) in modules {
                    if let Some(((pos, _), _)) = ray.intersect(module) {
                        chart.plotting_area().draw(&plotters::element::Circle::new(
                            (pos.0, pos.1),
                            3,
                            GREEN.filled(),
                        ))?;
                    }
                }
            }
        }
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
                (point.data[0], point.data[1]),
                3,
                BLACK.filled(),
            ))?;
        }
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Saved output to file: {}", filename);

    Ok(())
}

pub fn plot_2d_xz(
    filename: &str,
    plot_dims: &PlotDims,
    modules: &HashMap<u64, DetectorModule>,
    rays: Option<&Vec<StraightRay>>,
    plot_real_intersections: bool,
    hits: Option<&Vec<Hit>>,
    clustered_points: Option<&Vec<Vector3<f64>>>,
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

    if let Some(rays) = rays {
        for ray in rays {
            chart.draw_series(std::iter::once(PathElement::new(
                [
                    (0., 0.),
                    (
                        at_max_radius(ray, &plot_dims)?.0,
                        at_max_radius(ray, &plot_dims)?.2,
                    ),
                ],
                BLACK,
            )))?;
        }

        if plot_real_intersections {
            for ray in rays {
                for (_, module) in modules {
                    if let Some(((pos, _), _)) = ray.intersect(module) {
                        chart.plotting_area().draw(&plotters::element::Circle::new(
                            (pos.0, pos.2),
                            3,
                            GREEN.filled(),
                        ))?;
                    }
                }
            }
        }
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
                (point.data[0], point.data[2]),
                3,
                BLACK.filled(),
            ))?;
        }
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Saved output to file: {}", filename);

    Ok(())
}

pub fn plot_2d_yz(
    filename: &str,
    plot_dims: &PlotDims,
    modules: &HashMap<u64, DetectorModule>,
    rays: Option<&Vec<StraightRay>>,
    plot_real_intersections: bool,
    hits: Option<&Vec<Hit>>,
    clustered_points: Option<&Vec<Vector3<f64>>>,
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

    if let Some(rays) = rays {
        for ray in rays {
            chart.draw_series(std::iter::once(PathElement::new(
                [
                    (0., 0.),
                    (
                        at_max_radius(ray, &plot_dims)?.1,
                        at_max_radius(ray, &plot_dims)?.2,
                    ),
                ],
                BLACK,
            )))?;
        }

        if plot_real_intersections {
            for ray in rays {
                for (_, module) in modules {
                    if let Some(((pos, _), _)) = ray.intersect(module) {
                        chart.plotting_area().draw(&plotters::element::Circle::new(
                            (pos.1, pos.2),
                            3,
                            GREEN.filled(),
                        ))?;
                    }
                }
            }
        }
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
                (point.data[1], point.data[2]),
                3,
                BLACK.filled(),
            ))?;
        }
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Saved output to file: {}", filename);

    Ok(())
}

pub fn plot_all(
    file_dir: &str,
    plot_dims: &PlotDims,
    modules: &HashMap<u64, DetectorModule>,
    rays: Option<&Vec<StraightRay>>,
    plot_real_intersections: bool,
    hits: Option<&Vec<Hit>>,
    clustered_points: Option<&Vec<Vector3<f64>>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let filename_3d = file_dir.to_string() + "/3d.svg";
    let filename_2d_xy = file_dir.to_string() + "/2d_xy.svg";
    let filename_2d_xz = file_dir.to_string() + "/2d_xz.svg";
    let filename_2d_yz = file_dir.to_string() + "/2d_yz.svg";

    plot_3d(
        &filename_3d,
        plot_dims,
        modules,
        rays,
        plot_real_intersections,
        hits,
        clustered_points,
    )?;
    plot_2d_xy(
        &filename_2d_xy,
        plot_dims,
        modules,
        rays,
        plot_real_intersections,
        hits,
        clustered_points,
    )?;
    plot_2d_xz(
        &filename_2d_xz,
        plot_dims,
        modules,
        rays,
        plot_real_intersections,
        hits,
        clustered_points,
    )?;
    plot_2d_yz(
        &filename_2d_yz,
        plot_dims,
        modules,
        rays,
        plot_real_intersections,
        hits,
        clustered_points,
    )?;

    Ok(())
}

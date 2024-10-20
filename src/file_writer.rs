use std::fs;
use std::io::Write;

use crate::clustering::Hit;
use crate::event_generator::StraightRay;

const XCONST: (f64, f64, f64) = (0.0, -90.0, 0.0);
const YCONST: (f64, f64, f64) = (90.0, 0.0, 0.0);
const ZCONST: (f64, f64, f64) = (0.0, 0.0, 0.0);

pub fn create_box_detector(
    layers: &Vec<f64>,
    filename: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut id: u64 = 0;
    let mut file = fs::File::create(filename)?;

    let rotations = [XCONST, XCONST, YCONST, YCONST, ZCONST, ZCONST];
    for &radius in layers {
        let translations = [
            (-radius, -radius, -radius),
            (radius, -radius, -radius),
            (-radius, -radius, -radius),
            (-radius, radius, -radius),
            (-radius, -radius, -radius),
            (-radius, -radius, radius),
        ];
        let dim = radius * 2.;
        for (transl, rot) in translations.iter().zip(rotations) {
            file.write_fmt(format_args!(
                "{} {} {} {} {} {} {} {} {} {} {}\n",
                id,
                dim,
                dim,
                transl.0,
                transl.1,
                transl.2,
                dim as u64 * 5,
                dim as u64 * 5,
                rot.0,
                rot.1,
                rot.2
            ))?;
            id += 1;
        }
    }
    Ok(())
}

pub fn write_rays(filename: &str, rays: &Vec<StraightRay>) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = fs::File::create(filename)?;
    for ray in rays {
        file.write_fmt(format_args!("{} {} {}\n", ray.energy, ray.theta, ray.phi,))?;
    }
    Ok(())
}

pub fn write_hits(filename: &str, hits: &Vec<Hit>) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = fs::File::create(filename)?;
    for hit in hits {
        file.write_fmt(format_args!(
            "{} {} {}\n",
            hit.module_id, hit.pos.0, hit.pos.1,
        ))?;
    }
    Ok(())
}

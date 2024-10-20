use crate::clustering::Hit;
use crate::detector_module::DetectorModule;
use crate::event_generator::StraightRay;
use crate::matrix::{Matrix3, Vector3};
use std::collections::HashMap;
use std::fs;

pub fn read_modules(
    filename: &str,
) -> Result<HashMap<u64, DetectorModule>, Box<dyn std::error::Error>> {
    let mut result = HashMap::new();
    for line in fs::read_to_string(filename)?.lines() {
        let line_vec: Vec<_> = line.split_whitespace().collect();
        if line_vec.len() < 11 {
            let err: Box<dyn std::error::Error> = format!(
                "Failure to read file: {}. Line {} did not have enough parameters.",
                filename, line
            )
            .into();
            return Err(err);
        }
        let id: u64 = line_vec[0].parse()?;
        let dims = (line_vec[1].parse()?, line_vec[2].parse()?);
        let transl = Vector3::new(
            line_vec[3].parse()?,
            line_vec[4].parse()?,
            line_vec[5].parse()?,
        );
        let pixel_dims = (line_vec[6].parse()?, line_vec[7].parse()?);
        let angles: (f64, f64, f64) = (
            line_vec[8].parse()?,
            line_vec[9].parse()?,
            line_vec[10].parse()?,
        );
        let rot = Matrix3::from_angles(
            angles.0.to_radians(),
            angles.1.to_radians(),
            angles.2.to_radians(),
        )?;
        let module = DetectorModule::new(dims, pixel_dims, transl, rot)?;
        result.insert(id, module);
    }
    Ok(result)
}

pub fn read_hits(filename: &str) -> Result<Vec<Hit>, Box<dyn std::error::Error>> {
    let mut result = Vec::new();
    for line in fs::read_to_string(filename)?.lines() {
        let line_vec: Vec<_> = line.split_whitespace().collect();
        if line_vec.len() < 3 {
            let err: Box<dyn std::error::Error> = format!(
                "Failure to read file: {}. Line {} did not have enough parameters.",
                filename, line
            )
            .into();
            return Err(err);
        }
        let module_id = line_vec[0].parse()?;
        let pos = (line_vec[1].parse()?, line_vec[2].parse()?);
        result.push(Hit { module_id, pos });
    }
    Ok(result)
}

pub fn read_rays(filename: &str) -> Result<Vec<StraightRay>, Box<dyn std::error::Error>> {
    let mut result = Vec::new();
    for line in fs::read_to_string(filename)?.lines() {
        let line_vec: Vec<_> = line.split_whitespace().collect();
        if line_vec.len() < 3 {
            let err: Box<dyn std::error::Error> = format!(
                "Failure to read file: {}. Line {} did not have enough parameters.",
                filename, line
            )
            .into();
            return Err(err);
        }
        result.push(StraightRay {
            energy: line_vec[0].parse()?,
            theta: line_vec[1].parse()?,
            phi: line_vec[2].parse()?,
        });
    }
    Ok(result)
}

use crate::detector_module::DetectorModule;
use crate::matrix::*;
use std::fs;

pub fn read_modules(filename: &str) -> Result<Vec<DetectorModule>, Box<dyn std::error::Error>> {
    let mut vec = Vec::new();
    for line in fs::read_to_string(filename)?.lines() {
        let line_vec: Vec<_> = line.split_whitespace().collect();
        if line_vec.len() < 11 {
            // return Err(format!("Failure to read file: line {} did not have enough parameters.", line));
        }
        let id: u64 = line_vec[0].parse()?;
        let dims = (line_vec[1].parse()?, line_vec[2].parse()?);
        let transl = [
            line_vec[3].parse()?,
            line_vec[4].parse()?,
            line_vec[5].parse()?,
        ]
        .into_vector3()?;
        let pixels_dims = (line_vec[6].parse()?, line_vec[7].parse()?);
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
        let module = DetectorModule::new(id, dims, pixels_dims, vec![], transl, rot)?;
        vec.push(module);
    }
    Ok(vec)
}

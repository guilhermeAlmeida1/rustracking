// generate rays (from 1 point (assume starting from origin) or from angles)
// find where they intersect all the detector planes (from translation + rotation)
// check whether intersections are within detector boundaries
// create hits with a random spread

use crate::clustering::Hit;
use crate::config::ENERGY_HIT_FACTOR;
use crate::detector_module::{DetectorModule, PixelPosition};
use crate::matrix::Vector3;
use crate::particle::Particle;
use crate::plotting::Boundary;
use rand::distributions::Distribution;
use rand::Rng;
use std::collections::HashMap;
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug)]
pub struct StraightRay {
    pub theta: f64, // in radians
    pub phi: f64,   // in radians
    pub origin: Vector3<f64>,
}

impl StraightRay {
    pub fn new(theta: f64, phi: f64, origin: Vector3<f64>) -> StraightRay {
        StraightRay { theta, phi, origin }
    }

    // Returns a 3D point of intersection between the ray and the module
    pub fn intersect(
        &self,
        module: &DetectorModule,
        max_radius: f64,
    ) -> Option<((f64, f64, f64), PixelPosition)> {
        /*
         * Any point along the ray can be defined as r * e_r + origin,
         * with: e_r = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
         *
         * Any point inside a detector module can be defined as:
         * translation + rotation * (x, y, 0)
         * with x, y ranging from 0 to the module dimensions.
         *
         * This function equalizes the two, analytically solving the equation:
         * r * e_r + origin = translation + rotation * (x, y, 0) <=>
         *
         *     |r sin(theta)cos(phi) + K|     |A B 0|   |x|     |G|
         * <=> |r sin(theta)sin(phi) + L|  =  |C D 0| * |y|  +  |H|
         *     |r cos(theta)         + M|     |E F 0|   |0|     |I|
         *
         * see the demonstration in: misc/interset_proof.pdf
         */
        let theta = self.theta;
        let phi = self.phi;
        let k = self.origin[0];
        let l = self.origin[1];
        let m = self.origin[2];
        let dims = module.dims();
        let pixel_dims = module.pixel_dims();
        let rotation = module.rotation();
        let translation = module.translation();
        let a = rotation[0];
        let b = rotation[1];
        let c = rotation[3];
        let d = rotation[4];
        let e = rotation[6];
        let f = rotation[7];
        let g = translation[0];
        let h = translation[1];
        let i = translation[2];

        let r = ((b * e - a * f) * (h - l) + (c * f - d * e) * (g - k) + (a * d - b * c) * (i - m))
            / (theta.cos() * (a * d - b * c)
                + theta.sin() * ((c * f - d * e) * phi.cos() + (b * e - a * f) * phi.sin()));

        if r > 0. && r < max_radius {
            let (x, y, z) = self.at_radius(r);
            // Three different factorizations of the solution need to be calculated
            // because depending on the rotation matrix up to 2 of these are
            // invalid because they result in a division by 0
            if (a * d - b * c - 0.).abs() > 10. * std::f64::EPSILON {
                let module_x = (d * (x - g) + b * (h - y)) / (a * d - b * c);
                let module_y = (c * (x - g) + a * (h - y)) / (b * c - a * d);

                if module_x >= 0. && module_x <= dims.0 && module_y >= 0. && module_y <= dims.1 {
                    let pixel_pos = (
                        (module_x * pixel_dims.0 as f64 / dims.0) as u64,
                        (module_y * pixel_dims.1 as f64 / dims.1) as u64,
                    );
                    return Some(((x, y, z), pixel_pos));
                }
            } else if (a * f - b * e - 0.).abs() > 10. * std::f64::EPSILON {
                let module_x = (f * (x - g) + b * (i - z)) / (a * f - b * e);
                let module_y = (e * (x - g) + a * (i - z)) / (b * e - a * f);

                if module_x >= 0. && module_x <= dims.0 && module_y >= 0. && module_y <= dims.1 {
                    let pixel_pos = (
                        (module_x * pixel_dims.0 as f64 / dims.0) as u64,
                        (module_y * pixel_dims.1 as f64 / dims.1) as u64,
                    );
                    return Some(((x, y, z), pixel_pos));
                }
            } else if (c * f - d * e - 0.).abs() > 10. * std::f64::EPSILON {
                let module_x = (f * (y - h) + d * (i - z)) / (c * f - d * e);
                let module_y = (e * (y - h) + c * (i - z)) / (d * e - c * f);

                if module_x >= 0. && module_x <= dims.0 && module_y >= 0. && module_y <= dims.1 {
                    let pixel_pos = (
                        (module_x * pixel_dims.0 as f64 / dims.0) as u64,
                        (module_y * pixel_dims.1 as f64 / dims.1) as u64,
                    );
                    return Some(((x, y, z), pixel_pos));
                }
            }
        }
        None
    }

    pub fn at_radius(&self, r: f64) -> (f64, f64, f64) {
        (
            r * self.theta.sin() * self.phi.cos() + self.origin[0],
            r * self.theta.sin() * self.phi.sin() + self.origin[1],
            r * self.theta.cos() + self.origin[2],
        )
    }
}

pub enum Distributions {
    Uniform(rand::distributions::Uniform<f64>),
    Gauss(Gauss),
    // Poisson(f64),
    Exponential(Exponential),
}

// I would expect there to be a better way to do this
impl Distribution<f64> for Distributions {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        match self {
            Self::Uniform(dist) => dist.sample(rng),
            Self::Gauss(dist) => dist.sample(rng),
            Self::Exponential(dist) => dist.sample(rng),
        }
    }
}

pub struct Gauss {
    mean: f64,
    sigma: f64,
}
impl Distribution<f64> for Gauss {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        let r: [f64; 12] = rng.gen();
        let result = r.iter().sum::<f64>() - 6.;
        let result = result * self.sigma + self.mean;
        result
    }
}
impl Gauss {
    pub fn new(mean: f64, sigma: f64) -> Self {
        Self { mean, sigma }
    }
}
pub struct Exponential {
    lambda: f64,
}
impl Exponential {
    pub fn new(lambda: f64) -> Self {
        Self { lambda }
    }
}
impl Distribution<f64> for Exponential {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        // Apply the inverse commulative distribution function
        // to a random number between 0 and 1
        let r: f64 = rng.gen();
        -f64::ln(1. - r) / self.lambda
    }
}

pub fn generate_random_particles(
    total_energy: f64,
    min_energy: f64,
    dist: Distributions,
) -> Result<Vec<Particle>, &'static str> {
    let mut result = Vec::new();
    let mut rng = rand::thread_rng();

    let angle_dist = rand::distributions::Uniform::new(-PI, PI);
    let origin_dist = Gauss::new(0., 0.1);

    let mut current_energy = 0.;
    let mut switch = false;
    while current_energy < total_energy {
        let energy = dist.sample(&mut rng);
        let origin = Vector3::new(
            origin_dist.sample(&mut rng),
            origin_dist.sample(&mut rng),
            origin_dist.sample(&mut rng),
        );
        current_energy += energy;
        if energy < min_energy {
            continue;
        }
        let theta = angle_dist.sample(&mut rng);
        let phi = angle_dist.sample(&mut rng);
        let ray = StraightRay::new(theta, phi, origin);
        if switch {
            result.push(Particle::new_electron(energy, ray)?);
        } else {
            result.push(Particle::new_proton(energy, ray)?);
        }
        switch = !switch;
    }
    Ok(result)
}

pub fn create_hits(
    mut particles: Vec<Particle>,
    modules: &HashMap<u64, DetectorModule>,
    bounds: &Boundary,
) -> Vec<Hit> {
    static DIST: Gauss = Gauss {
        mean: 0.,
        sigma: 0.15,
    };
    let mut rng = rand::thread_rng();

    let mut hits = Vec::new();
    for particle in &mut particles {
        while particle.validate_energy() {
            for (id, module) in modules {
                if let Some((_, centre_pixel)) = particle.intersect(module) {
                    let mut partial = Vec::new();
                    let pixel_dims = module.pixel_dims();
                    let dims = module.dims();
                    let num_hits =
                        1 + (5. * (1. - f64::exp(-particle.energy / ENERGY_HIT_FACTOR))) as u64;
                    for _ in 0..num_hits {
                        let x = (centre_pixel.0 as f64
                            + (DIST.sample(&mut rng) * pixel_dims.0 as f64 / dims.0))
                            .round() as u64;
                        let y = (centre_pixel.1 as f64
                            + (DIST.sample(&mut rng) * pixel_dims.1 as f64 / dims.1))
                            .round() as u64;
                        let hit = Hit::new(*id, (x, y));
                        if x >= pixel_dims.0 || y >= pixel_dims.1 || partial.contains(&hit) {
                            continue;
                        }
                        partial.push(hit);
                    }
                    hits.append(&mut partial);
                }
            }
            particle.do_step();
        }
        let origin = particle.ray.origin;
        if origin[0] < bounds.0.start || origin[0] > bounds.0.end {
            continue;
        }
        if origin[1] < bounds.1.start || origin[1] > bounds.1.end {
            continue;
        }
        if origin[2] < bounds.2.start || origin[2] > bounds.2.end {
            continue;
        }
    }
    hits
}

pub fn generate_random_event(
    total_energy: f64,
    min_energy: f64,
    dist: Distributions,
    modules: &HashMap<u64, DetectorModule>,
    bounds: &Boundary,
) -> Result<(Vec<Particle>, Vec<Hit>), &'static str> {
    let particles = generate_random_particles(total_energy, min_energy, dist)?;
    let hits = create_hits(particles.clone(), modules, bounds);
    Ok((particles, hits))
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::matrix::{Matrix3, Vector3};

    #[test]
    #[ignore]
    fn gauss_distribution() {
        const TOLERANCE: f64 = 0.01;
        let mut rng = rand::thread_rng();
        let mean = 3.7;
        let sigma = 2.5;
        let dist = Gauss { mean, sigma };

        // Workaround for a flaky test
        let mut assertion = false;
        for _ in 0..20 {
            let mut vec = vec![0.; 100000];
            vec.iter_mut()
                .for_each(|val: &mut f64| *val = dist.sample(&mut rng));

            let mean_val = vec.iter().sum::<f64>() / vec.len() as f64;
            let sigma_val = (vec
                .iter()
                .map(|val| (*val - mean_val) * (*val - mean_val))
                .sum::<f64>()
                / vec.len() as f64)
                .sqrt();

            if (mean - mean_val).abs() < TOLERANCE && (sigma - sigma_val).abs() < TOLERANCE {
                assertion = true;
                break;
            }
        }
        assert!(assertion);
    }

    #[test]
    #[ignore]
    fn exponential_distribution() {
        const TOLERANCE: f64 = 0.01;
        let mut rng = rand::thread_rng();
        let lambda = 3.7;
        let dist = Exponential { lambda };
        let mean_expected = 1. / lambda;
        let sigma_expected = 1. / lambda;

        // Workaround for a flaky test
        let mut assertion = false;
        for _ in 0..20 {
            let mut vec = vec![0.; 100000];
            vec.iter_mut()
                .for_each(|val: &mut f64| *val = dist.sample(&mut rng));

            let mean_val = vec.iter().sum::<f64>() / vec.len() as f64;
            let sigma_val = (vec
                .iter()
                .map(|val| (*val - mean_val) * (*val - mean_val))
                .sum::<f64>()
                / vec.len() as f64)
                .sqrt();

            if (mean_expected - mean_val).abs() < TOLERANCE
                && (sigma_expected - sigma_val).abs() < TOLERANCE
            {
                assertion = true;
                break;
            }
        }
        assert!(assertion);
    }

    #[test]
    fn at_radius() {
        let ray = StraightRay::new(PI / 2., PI / 2., Vector3::new(0., 0., 0.));
        let result = ray.at_radius(5.);
        let expected = (0., 5., 0.);

        assert!((result.0 - expected.0).abs() < 10. * std::f64::EPSILON);
        assert!((result.1 - expected.1).abs() < 10. * std::f64::EPSILON);
        assert!((result.2 - expected.2).abs() < 10. * std::f64::EPSILON);
    }

    #[test]
    fn intersect() {
        let rotation = Matrix3::from_angles(0., -PI / 2., 0.).unwrap(); // xconst
        let translation = Vector3::new(10., -5., -5.);
        let module = DetectorModule::new((10., 10.), (10, 10), translation, rotation).unwrap();

        let ray = StraightRay::new(PI / 2., 0., Vector3::new(0., 0., 0.));
        let (result, pixel) = ray.intersect(&module, std::f64::MAX).unwrap();

        let expected = (10., 0., 0.);
        println!("result: {:?}", result);
        assert!((result.0 - expected.0).abs() < 10. * std::f64::EPSILON);
        assert!((result.1 - expected.1).abs() < 10. * std::f64::EPSILON);
        assert!((result.2 - expected.2).abs() < 10. * std::f64::EPSILON);
        assert_eq!(pixel, (5, 5));

        let ray = StraightRay::new(PI / 2., 0., Vector3::new(0., 0., 0.));
        assert!(ray.intersect(&module, 9.).is_none());

        let ray = StraightRay::new(0., 0., Vector3::new(0., 0., 0.));
        assert!(ray.intersect(&module, std::f64::MAX).is_none());

        let translation = Vector3::new(10., 1., 1.);
        let module = DetectorModule::new((10., 10.), (10, 10), translation, rotation).unwrap();

        let ray = StraightRay::new(PI / 2., 0., Vector3::new(0., 0., 0.));
        assert!(ray.intersect(&module, std::f64::MAX).is_none());

        let rotation = Matrix3::from_angles(PI / 2., 0., 0.).unwrap(); // yconst
        let translation = Vector3::new(-5., 10., -5.);
        let module = DetectorModule::new((10., 10.), (10, 10), translation, rotation).unwrap();

        let ray = StraightRay::new(PI / 2., PI / 2., Vector3::new(0., 0., 0.));
        let (result, pixel) = ray.intersect(&module, std::f64::MAX).unwrap();

        let expected = (0., 10., 0.);
        println!("result: {:?}", result);
        assert!((result.0 - expected.0).abs() < 10. * std::f64::EPSILON);
        assert!((result.1 - expected.1).abs() < 10. * std::f64::EPSILON);
        assert!((result.2 - expected.2).abs() < 10. * std::f64::EPSILON);
        assert_eq!(pixel, (5, 5));

        let rotation = Matrix3::identity(); // zconst
        let translation = Vector3::new(-5., -5., 10.);
        let module = DetectorModule::new((10., 10.), (10, 10), translation, rotation).unwrap();

        let ray = StraightRay::new(0., 0., Vector3::new(0., 0., 0.));
        let (result, pixel) = ray.intersect(&module, std::f64::MAX).unwrap();

        let expected = (0., 0., 10.);
        println!("result: {:?}", result);
        assert!((result.0 - expected.0).abs() < 10. * std::f64::EPSILON);
        assert!((result.1 - expected.1).abs() < 10. * std::f64::EPSILON);
        assert!((result.2 - expected.2).abs() < 10. * std::f64::EPSILON);
        assert_eq!(pixel, (5, 5));

        let ray = StraightRay::new(0., 0., Vector3::new(3., 3., 2.));
        let (result, _) = ray.intersect(&module, std::f64::MAX).unwrap();

        let expected = (3., 3., 10.);
        println!("result: {:?}", result);
        assert!((result.0 - expected.0).abs() < 10. * std::f64::EPSILON);
        assert!((result.1 - expected.1).abs() < 10. * std::f64::EPSILON);
        assert!((result.2 - expected.2).abs() < 10. * std::f64::EPSILON);
    }
}

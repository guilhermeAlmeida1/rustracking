// generate rays (from 1 point (assume starting from origin) or from angles)
// find where they intersect all the detector planes (from translation + rotation)
// check whether intersections are within detector boundaries
// create hits with a random spread

use rand::distributions::Distribution;
use rand::Rng;

pub struct Ray {
    pub theta: f64,
    pub phi: f64,
    // pub charge: f64,
    pub energy: f64,
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

///
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
        Self {mean, sigma}
    }
}
pub struct Exponential {
    lambda: f64,
}
impl Exponential {
    pub fn new(lambda: f64) -> Self {
        Self {lambda}
    }
}
impl Distribution<f64> for Exponential {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        let r: f64 = rng.gen();
        -f64::ln(1. - r) / self.lambda
    }
}

pub fn generate_random_rays(total_energy: f64, min_energy: f64, dist: Distributions) -> Vec<Ray> {
    use std::f64::consts::PI;
    let mut result = Vec::new();
    let mut rng = rand::thread_rng();

    let angle_dist = rand::distributions::Uniform::new(-PI, PI);

    let mut current_energy = 0.;
    while current_energy < total_energy {
        let energy = dist.sample(&mut rng);
        current_energy += energy;
        if energy < min_energy { // do not store these values
            continue;
        }
        let theta = angle_dist.sample(&mut rng);
        let phi = angle_dist.sample(&mut rng);
        result.push(Ray {theta, phi, energy});
    }
    result
}

#[cfg(test)]
mod test {
    use super::*;

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
}

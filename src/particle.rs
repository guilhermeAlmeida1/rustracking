use crate::event_generator::StraightRay;
use crate::field::{ConstantField, Field};

pub const STEP_TIME: f64 = 1e-10; // in s

const ELECTRON_CHARGE: f64 = -1.; // in q_e
const ELECTRON_REST_MASS: f64 = 0.5110; // in MeV/c^2

const PROTON_CHARGE: f64 = -ELECTRON_CHARGE; // in q_e
const PROTON_REST_MASS: f64 = 938.8; // in MeV/c^2

pub const SPEED_OF_LIGHT: f64 = 299792458.; // in unit distance/s
pub const ENERGY_LOSS_PER_SECOND: f64 = 1e11; // in MeV/s

pub const CONST_MAGNETIC_FIELD: ConstantField = ConstantField::new(1e8, 0., 0.); // in TBD units

// galmeida: change this whole file to use crate::matrix::Vector3

fn to_cartesian(polar: [f64; 3]) -> [f64; 3] {
    [
        polar[0] * polar[1].sin() * polar[2].cos(),
        polar[0] * polar[1].sin() * polar[2].sin(),
        polar[0] * polar[1].cos(),
    ]
}

fn to_polar(cartesian: [f64; 3]) -> [f64; 3] {
    let r = f64::sqrt(cartesian.iter().map(|v| v.powi(2)).sum());
    if (r - 0.).abs() < f64::EPSILON {
        return [0., 0., 0.];
    }
    let x = cartesian[0];
    let y = cartesian[1];
    let z = cartesian[2];

    let theta = if z == 0. {
        std::f64::consts::FRAC_PI_2
    } else if z < 0. {
        (x.powi(2) + y.powi(2)) / z
    } else {
        (x.powi(2) + y.powi(2)) / z + std::f64::consts::PI
    };

    let phi = if x == 0. {
        if y >= 0. {
            std::f64::consts::FRAC_PI_2
        } else {
            -std::f64::consts::FRAC_PI_2
        }
    } else if x > 0. {
        (y / x).atan()
    } else {
        if y >= 0. {
            (y / x).atan() + std::f64::consts::PI
        } else {
            (y / x).atan() - std::f64::consts::PI
        }
    };
    [r, theta, phi]
}

// cross product between two vectors in cartesian coordinates
fn cross_product(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let [a1, a2, a3] = a;
    let [b1, b2, b3] = b;
    [a2 * b3 - b2 * a3, a3 * b1 - b3 * a1, a1 * b2 - b1 * a2]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn mult_scalar(a: [f64; 3], b: f64) -> [f64; 3] {
    a.map(|v| v * b)
}

#[derive(Debug)]
pub struct Particle {
    pub charge: f64,
    pub rest_mass: f64,
    pub ray: StraightRay,
    velocity_cartesian: [f64; 3],
}

impl Particle {
    pub fn new_electron(ray: StraightRay) -> Result<Particle, &'static str> {
        let mut result = Self {
            charge: ELECTRON_CHARGE,
            rest_mass: ELECTRON_REST_MASS,
            ray,
            velocity_cartesian: [0., 0., 0.],
        };
        result.velocity_cartesian = to_cartesian(result.velocity_polar()?);
        Ok(result)
    }

    pub fn new_proton(ray: StraightRay) -> Result<Particle, &'static str> {
        let mut result = Self {
            charge: PROTON_CHARGE,
            rest_mass: PROTON_REST_MASS,
            ray,
            velocity_cartesian: [0., 0., 0.],
        };
        result.velocity_cartesian = to_cartesian(result.velocity_polar()?);
        Ok(result)
    }

    pub fn validate_energy(&self) -> bool {
        self.ray.energy >= self.rest_mass
    }

    // in MeV/c
    pub fn abs_momentum(&self) -> Result<f64, &'static str> {
        if !self.validate_energy() {
            return Err("Particle's energy may not be smaller than it's rest mass.");
        }
        Ok(f64::sqrt(
            self.ray.energy.powi(2) - (self.rest_mass).powi(2),
        ))
    }

    // in c
    pub fn abs_velocity(&self) -> Result<f64, &'static str> {
        let momentum = self.abs_momentum()?;
        Ok(momentum / f64::sqrt(self.rest_mass.powi(2) + momentum.powi(2)))
    }

    fn velocity_polar(&self) -> Result<[f64; 3], &'static str> {
        Ok([self.abs_velocity()?, self.ray.theta, self.ray.phi])
    }

    fn acceleration_cartesian(&self, magnetic_field: impl Field) -> Result<[f64; 3], &'static str> {
        let b_field = magnetic_field.at(self.ray.origin.into());
        let cross = cross_product(self.velocity_cartesian, b_field);

        Ok(cross.map(|val| val * self.charge))
    }

    pub fn do_step(&mut self) -> Result<(), &'static str> {
        let accel = self.acceleration_cartesian(CONST_MAGNETIC_FIELD)?; // acceleration vector in natural units
                                                                        // println!("accel: {:?}", accel);

        let linear_term = mult_scalar(self.velocity_cartesian, STEP_TIME * SPEED_OF_LIGHT);
        self.ray.origin = add(self.ray.origin.into(), linear_term).into(); // O(t^2) term (accel) is not considered

        let delta_v = mult_scalar(accel, STEP_TIME);
        self.velocity_cartesian = add(self.velocity_cartesian, delta_v);
        // this will always be gaining energy
        let new_vel_polar = to_polar(self.velocity_cartesian);

        self.ray.theta = new_vel_polar[1];
        self.ray.phi = new_vel_polar[2];
        self.ray.energy -= ENERGY_LOSS_PER_SECOND * STEP_TIME;

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::utils::assert_near;

    fn _new_basic_proton(energy: f64) -> Result<Particle, &'static str> {
        Particle::new_proton(StraightRay::new(0., 0., (0., 0., 0.), energy))
    }

    #[test]
    fn new_basic_proton() {
        let particle = _new_basic_proton(5000.).unwrap();
        assert_near(particle.charge, PROTON_CHARGE, 0.1);
        assert_near(particle.rest_mass, PROTON_REST_MASS, 0.1);
        assert_near(particle.ray.energy, 5000., 0.1);
    }

    #[test]
    fn to_polar_() {
        let cartesian = [f64::sqrt(2.) / 2., f64::sqrt(2.) / 2., 0.];
        let result = to_polar(cartesian);
        let expected = [1., std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_4];
        for i in 0..3 {
            assert_near(result[i], expected[i], 10. * f64::EPSILON);
        }
    }

    #[test]
    fn cross_product_() {
        let result = cross_product([1., 0., 0.], [1., 2., 4.]);
        let expected = [0., -4., 2.];
        for i in 0..3 {
            assert_near(result[i], expected[i], 10. * f64::EPSILON);
        }
    }

    #[test]
    fn given_energy_less_than_rest_mass_then_ctor_returns_error() {
        assert!(_new_basic_proton(PROTON_REST_MASS - 1.).is_err());
    }

    #[test]
    fn given_energy_less_than_rest_mass_then_abs_momentum_returns_error() {
        let mut particle = _new_basic_proton(1000.).unwrap();
        particle.ray.energy = 0.;
        assert!(particle.abs_momentum().is_err());
    }

    #[test]
    fn abs_momentum_zero() {
        let particle = _new_basic_proton(PROTON_REST_MASS).unwrap();
        assert_near(particle.abs_momentum().unwrap(), 0., 1.);
    }

    #[test]
    fn abs_momentum() {
        let particle = _new_basic_proton(5000.).unwrap();
        assert_near(particle.abs_momentum().unwrap(), 4911., 1.);
    }

    #[test]
    fn abs_velocity_zero() {
        let particle = _new_basic_proton(PROTON_REST_MASS).unwrap();
        assert_near(particle.abs_velocity().unwrap(), 0., 0.001);
    }

    #[test]
    fn abs_velocity_max() {
        let particle = _new_basic_proton(f64::MAX).unwrap();
        assert_near(particle.abs_velocity().unwrap(), 1., 0.001);
    }
}

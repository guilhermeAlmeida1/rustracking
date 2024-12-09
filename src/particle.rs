use crate::config::{
    CONST_MAGNETIC_FIELD, ELECTRON_CHARGE, ELECTRON_REST_MASS, ENERGY_LOSS_PER_SECOND,
    PROTON_CHARGE, PROTON_REST_MASS, SPEED_OF_LIGHT, STEP_TIME,
};
use crate::detector_module::{DetectorModule, PixelPosition};
use crate::event_generator::StraightRay;
use crate::field::Field;
use crate::matrix::Vector3;

#[derive(Clone, Copy, Debug)]
pub struct Particle {
    pub charge: f64,
    pub rest_mass: f64,
    pub energy: f64,
    pub ray: StraightRay,
    velocity_cartesian: Vector3<f64>,
}

impl Particle {
    pub fn new(
        charge: f64,
        rest_mass: f64,
        energy: f64,
        ray: StraightRay,
    ) -> Result<Particle, &'static str> {
        let mut result = Self {
            charge,
            rest_mass,
            energy,
            ray,
            velocity_cartesian: Vector3::default(),
        };
        if !result.validate_energy() {
            return Err("Energy of particle must be larger than its rest mass.");
        }
        result.velocity_cartesian = result.velocity_polar().to_cartesian();
        Ok(result)
    }

    pub fn new_electron(energy: f64, ray: StraightRay) -> Result<Particle, &'static str> {
        Self::new(ELECTRON_CHARGE, ELECTRON_REST_MASS, energy, ray)
    }

    pub fn new_proton(energy: f64, ray: StraightRay) -> Result<Particle, &'static str> {
        Self::new(PROTON_CHARGE, PROTON_REST_MASS, energy, ray)
    }

    pub fn validate_energy(&self) -> bool {
        self.energy >= self.rest_mass
    }

    // in MeV/c
    pub fn abs_momentum(&self) -> f64 {
        f64::sqrt(self.energy.powi(2) - (self.rest_mass).powi(2))
    }

    // in c
    pub fn abs_velocity(&self) -> f64 {
        let momentum = self.abs_momentum();
        momentum / f64::sqrt(self.rest_mass.powi(2) + momentum.powi(2))
    }

    fn velocity_polar(&self) -> Vector3<f64> {
        Vector3::new(self.abs_velocity(), self.ray.theta, self.ray.phi)
    }

    fn acceleration_cartesian(&self, magnetic_field: impl Field) -> Vector3<f64> {
        let b_field = magnetic_field.at(self.ray.origin.into());
        let cross = self.velocity_cartesian.cross_cartesian(b_field);

        cross.map(|val| val * self.charge)
    }

    pub fn do_step(&mut self) -> () {
        let accel = self.acceleration_cartesian(CONST_MAGNETIC_FIELD); // acceleration vector in natural units
                                                                       // println!("accel: {:?}", accel);

        let linear_term = self.velocity_cartesian * (STEP_TIME * SPEED_OF_LIGHT);
        self.ray.origin = self.ray.origin + linear_term; // O(t^2) term (accel) is not considered

        let delta_v = accel * STEP_TIME;
        self.velocity_cartesian = self.velocity_cartesian + delta_v;
        let new_vel_spherical = self.velocity_cartesian.to_spherical();

        self.ray.theta = new_vel_spherical[1];
        self.ray.phi = new_vel_spherical[2];
        self.energy -= ENERGY_LOSS_PER_SECOND * STEP_TIME;
    }

    fn step_distance(&self) -> f64 {
        self.abs_velocity() * STEP_TIME * SPEED_OF_LIGHT
    }

    pub fn intersect(&self, module: &DetectorModule) -> Option<((f64, f64, f64), PixelPosition)> {
        self.ray.intersect(module, self.step_distance())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::utils::assert_near;

    fn _new_basic_proton(energy: f64) -> Result<Particle, &'static str> {
        Particle::new_proton(energy, StraightRay::new(0., 0., Vector3::new(0., 0., 0.)))
    }

    #[test]
    fn new_basic_proton() {
        let particle = _new_basic_proton(5000.).unwrap();
        assert_near(particle.charge, PROTON_CHARGE, 0.1);
        assert_near(particle.rest_mass, PROTON_REST_MASS, 0.1);
        assert_near(particle.energy, 5000., 0.1);
    }

    #[test]
    fn given_energy_less_than_rest_mass_then_ctor_returns_error() {
        assert!(_new_basic_proton(PROTON_REST_MASS - 1.).is_err());
    }

    #[test]
    fn abs_momentum_zero() {
        let particle = _new_basic_proton(PROTON_REST_MASS).unwrap();
        assert_near(particle.abs_momentum(), 0., 1.);
    }

    #[test]
    fn abs_momentum() {
        let particle = _new_basic_proton(5000.).unwrap();
        assert_near(particle.abs_momentum(), 4911., 1.);
    }

    #[test]
    fn abs_velocity_zero() {
        let particle = _new_basic_proton(PROTON_REST_MASS).unwrap();
        assert_near(particle.abs_velocity(), 0., 0.001);
    }

    #[test]
    fn abs_velocity_max() {
        let particle = _new_basic_proton(f64::MAX).unwrap();
        assert_near(particle.abs_velocity(), 1., 0.001);
    }
}

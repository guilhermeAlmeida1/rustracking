use crate::field::ConstantField;

pub const STEP_TIME: f64 = 1e-10; // in s

pub const ELECTRON_CHARGE: f64 = -1.; // in q_e
pub const ELECTRON_REST_MASS: f64 = 0.5110; // in MeV/c^2

pub const PROTON_CHARGE: f64 = -ELECTRON_CHARGE; // in q_e
pub const PROTON_REST_MASS: f64 = 938.8; // in MeV/c^2

pub const SPEED_OF_LIGHT: f64 = 299792458.; // in unit distance/s
pub const ENERGY_LOSS_PER_SECOND: f64 = 1e11; // in MeV/s

pub const CONST_MAGNETIC_FIELD: ConstantField = ConstantField::new(5e7, 0., 0.); // in TBD units

pub const ENERGY_HIT_FACTOR: f64 = 700.; // in MeV

pub const SEEDING_DELTA_R: f64 = 3.; // in unit distance
pub const SEEDING_DELTA_THETA: f64 = 0.3; // in radian
pub const SEEDING_DELTA_PHI: f64 = SEEDING_DELTA_THETA; // in unit distance

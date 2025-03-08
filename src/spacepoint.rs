use crate::config::{SEEDING_DELTA_PHI, SEEDING_DELTA_R, SEEDING_DELTA_THETA};
use crate::matrix::Vector3;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpacePoint(pub Vector3<f64>);
impl SpacePoint {
    #[inline]
    pub fn new(data: Vector3<f64>) -> Self {
        Self { 0: data }
    }

    pub fn is_compatible(&self, other: &SpacePoint) -> bool {
        let self_spher = self.0.to_spherical();
        let other_spher = other.0.to_spherical();
        let delta_spher = (other_spher - self_spher).map(|v| v.abs());

        return delta_spher[0] < SEEDING_DELTA_R
            && (delta_spher[1] < SEEDING_DELTA_THETA
                || delta_spher[1] > std::f64::consts::PI - SEEDING_DELTA_THETA)
            && (delta_spher[2] < SEEDING_DELTA_PHI
                || delta_spher[2] > std::f64::consts::PI - SEEDING_DELTA_PHI);
    }
}

impl Into<Vector3<f64>> for SpacePoint {
    #[inline]
    fn into(self) -> Vector3<f64> {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn new_sp(x: f64, y: f64, z: f64) -> SpacePoint {
        Vector3::new(x, y, z).into()
    }

    #[test]
    fn is_compatible() {
        let sp0 = new_sp(0., 1., 1.);
        let sp1 = new_sp(1e-6, 1., 1.);
        let sp2 = new_sp(1e-6, 1.00001, 1.00001);

        let sp3 = sp0.clone();
        let sp4 = new_sp(1000., 1., 1.);
        let sp5 = new_sp(0., 5., 1.);
        let sp6 = new_sp(0., 1., 5.);

        assert!(sp0.is_compatible(&sp1));
        assert!(sp0.is_compatible(&sp2));
        assert!(sp1.is_compatible(&sp2));
        assert!(sp0.is_compatible(&sp3));

        assert!(!sp0.is_compatible(&sp4));
        assert!(!sp0.is_compatible(&sp5));
        assert!(!sp0.is_compatible(&sp6));
    }

    #[test]
    fn is_compatible_x_axis_edge() {
        let sp0 = new_sp(1e-6, 0., 0.);
        let sp1 = new_sp(-1e-6, 0., 0.);

        assert!(sp0.is_compatible(&sp1));
    }

    #[test]
    fn is_compatible_y_axis_edge() {
        let sp0 = new_sp(0., 1e-6, 0.);
        let sp1 = new_sp(0., -1e-6, 0.);

        assert!(sp0.is_compatible(&sp1));
    }

    #[test]
    fn is_compatible_z_axis_edge() {
        let sp0 = new_sp(0., 0., 1e-6);
        let sp1 = new_sp(0., 0., -1e-6);

        assert!(sp0.is_compatible(&sp1));
    }
}

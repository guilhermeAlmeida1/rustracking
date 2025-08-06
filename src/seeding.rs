use crate::config::{SEEDING_DELTA_PHI, SEEDING_DELTA_R, SEEDING_DELTA_THETA};
use crate::spacepoint::SpacePoint;

#[derive(Debug, PartialEq)]
pub struct Doublet<'a>(pub &'a SpacePoint, pub &'a SpacePoint);

impl<'a> Doublet<'a> {
    pub fn is_compatible(&self, other: &Doublet) -> bool {
        // assume self.1 == other.0
        let origin_spher = self.0 .0.to_spherical();
        let end_spher = other.1 .0.to_spherical();

        let delta_spher = (end_spher - origin_spher).map(|v| v.abs());

        assert!(origin_spher[0] < end_spher[0]);

        return delta_spher[0] > SEEDING_DELTA_R
            && delta_spher[0] < 2. * SEEDING_DELTA_R
            && (delta_spher[1] < 2. * SEEDING_DELTA_THETA
                || delta_spher[1] > std::f64::consts::PI - 2. * SEEDING_DELTA_THETA)
            && (delta_spher[2] < 2. * SEEDING_DELTA_PHI
                || delta_spher[2] > std::f64::consts::PI - 2. * SEEDING_DELTA_PHI);
    }
}

pub fn find_doublets<'a>(spacepoints: &'a Vec<SpacePoint>) -> Vec<Doublet<'a>> {
    if let Some(result) = spacepoints
        .iter()
        .map(|sp| {
            spacepoints
                .iter()
                .filter(move |sp2| sp != *sp2 && sp.is_compatible(sp2))
                .map(move |sp2| Doublet(sp, sp2))
                .collect::<Vec<Doublet>>()
        })
        .reduce(|mut v1, mut v2| {
            v1.append(&mut v2);
            v1
        })
    {
        return result;
    }
    return Vec::new();
}

#[derive(Debug, PartialEq)]
struct Seed<'a>(pub &'a SpacePoint, pub &'a SpacePoint, pub &'a SpacePoint);

impl<'a> Seed<'a> {}

pub fn find_seeds<'a>(doublets: &'a Vec<Doublet>) -> Vec<Seed<'a>> {
    if let Some(result) = doublets
        .iter()
        .map(|doub| {
            doublets
                .iter()
                .filter(move |doub2| doub != *doub2 && doub.is_compatible(doub2))
                .map(move |doub2| Seed(doub.0, doub.1, doub2.1))
                .collect::<Vec<Seed>>()
        })
        .reduce(|mut v1, mut v2| {
            v1.append(&mut v2);
            v1
        })
    {
        return result;
    }
    return Vec::new();
}

// Find the centre and radius of a sphere that passes through the origin and the 3 other spacepoints
// see the demonstration in: misc/find_sphere_proof.pdf
pub fn find_sphere(sp1: &SpacePoint, sp2: &SpacePoint, sp3: &SpacePoint) -> (SpacePoint, f64) {
    let a = sp1.0;
    let b = sp2.0;
    let c = sp3.0;

    let mut k: f64 = a[0] * (b[1] * c[2] - b[2] * c[1])
        + a[1] * (b[2] * c[0] - b[0] * c[2])
        + a[2] * (b[0] * c[1] - b[1] * c[0]);

    let mut l: f64 = a.norm_sq() * (b[2] * c[1] - b[1] * c[2])
        + b.norm_sq() * (a[1] * c[2] - a[2] * c[1])
        + c.norm_sq() * (a[2] * b[1] - a[1] * b[2]);

    let mut m: f64 = a.norm_sq() * (b[0] * c[2] - b[2] * c[0])
        + b.norm_sq() * (a[2] * c[0] - a[0] * c[2])
        + c.norm_sq() * (a[0] * b[2] - a[2] * b[0]);

    let mut n: f64 = a.norm_sq() * (b[1] * c[0] - b[0] * c[1])
        + b.norm_sq() * (a[0] * c[1] - a[1] * c[0])
        + c.norm_sq() * (a[1] * b[0] - a[0] * b[1]);

    if k < 0. {
        k = -k;
        l = -l;
        m = -m;
        n = -n;
    }

    let radius = (l * l + m * m + n * n).sqrt() / (2. * k);
    let centre = SpacePoint::new(-l / (2. * k), -m / (2. * k), -n / (2. * k));

    (centre, radius)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::matrix::Vector3;
    use crate::utils::assert_near;

    fn new_sp(x: f64, y: f64, z: f64) -> SpacePoint {
        Vector3::new(x, y, z).into()
    }

    #[test]
    fn find_doublets_() {
        let spacepoints = vec![
            new_sp(0., 0., 0.),
            new_sp(1., 2., 3.),
            new_sp(1., 2., 3.01),
            new_sp(1., 2.05, 3.),
            new_sp(0.01, 0., 0.01),
            new_sp(0., 0.01, 0.),
        ];

        let doublets = find_doublets(&spacepoints);

        let mut expected_doublets = Vec::new();
        for sp in &spacepoints {
            for sp2 in &spacepoints {
                if sp2 != sp && sp.is_compatible(sp2) {
                    expected_doublets.push(Doublet(&sp, &sp2));
                }
            }
        }

        assert!(
            expected_doublets.len() >= 5,
            "Use more significant test data."
        );

        assert_eq!(doublets.len(), expected_doublets.len());
        assert_eq!(doublets, expected_doublets);
    }

    #[test]
    fn test_find_sphere() {
        let sp1 = SpacePoint::new(1., 2., 1. + 10f64.sqrt());
        let sp2 = SpacePoint::new(2., 1., 1. - 12f64.sqrt());
        let sp3 = SpacePoint::new(3. + 13f64.sqrt(), 3., 1.);

        let (centre, radius) = find_sphere(&sp1, &sp2, &sp3);

        assert_near(centre.0[0], 3.0, 100. * f64::EPSILON);
        assert_near(centre.0[1], 2.0, 100. * f64::EPSILON);
        assert_near(centre.0[2], 1.0, 100. * f64::EPSILON);
        assert_near(radius.powi(2), 14., 100. * f64::EPSILON);
    }
}

use std::collections::HashMap;

use crate::config::{SEEDING_DELTA_PHI, SEEDING_DELTA_R, SEEDING_DELTA_THETA};
use crate::{detector_module::*, matrix::Vector3};

#[derive(Debug, PartialEq)]
pub struct Hit {
    pub module_id: u64,
    pub pos: PixelPosition,
}

impl Hit {
    pub fn new(module_id: u64, pos: PixelPosition) -> Self {
        Hit { module_id, pos }
    }

    fn adjacent(&self, other: &Hit) -> bool {
        if self.module_id != other.module_id {
            return false;
        }

        let mut a = [self.pos.0, other.pos.0];
        let mut b = [self.pos.1, other.pos.1];
        a.sort();
        b.sort();

        if a[1] - a[0] > 1 || b[1] - b[0] > 1 {
            return false;
        }
        true
    }
}

fn find_root(indices: &mut Vec<usize>, e: usize) -> usize {
    let mut r = e;
    while indices[r] != r {
        r = indices[r];
    }
    r
}

fn make_union(indices: &mut Vec<usize>, e1: usize, e2: usize) -> usize {
    match e1 < e2 {
        true => {
            indices[e2] = e1;
            return e1;
        }
        false => {
            indices[e1] = e2;
            return e2;
        }
    }
}

fn ccl(hits: &Vec<Hit>) -> Vec<usize> {
    let mut indices = vec![0; hits.len()];

    // first scan: pixel association
    for i in 0..indices.len() {
        indices[i] = i;
        let mut ai = i;
        for j in 0..i {
            if hits[i].adjacent(&hits[j]) {
                let root = find_root(&mut indices, j);
                ai = make_union(&mut indices, ai, root);
            }
        }
    }

    // second scan: transitive closure
    for i in 0..indices.len() {
        if indices[i] != i {
            indices[i] = indices[indices[i]];
        }
    }

    indices
}

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

pub fn clustering(
    modules: &HashMap<u64, DetectorModule>,
    hits: &Vec<Hit>,
) -> Result<Vec<SpacePoint>, &'static str> {
    let ccl = ccl(hits);

    // possibly the most hideously beautiful code I've created
    let result: Vec<SpacePoint> = ccl
        .iter()
        .enumerate()
        .filter(|(i, &val)| *i == val)
        .map(|(_, &it)| {
            ccl.iter()
                .enumerate()
                .filter(|(_, &it2)| it2 == it)
                .map(|(idx, _)| (hits[idx].module_id, hits[idx].pos, 1.))
                .reduce(|a, b| (a.0, (a.1 .0 + b.1 .0, a.1 .1 + b.1 .1), a.2 + b.2))
                .unwrap()
        })
        .map(|agg| {
            (&modules[&agg.0])
                .pixel_position_to_vector3((
                    agg.1 .0 as f64 / agg.2 + 0.5,
                    agg.1 .1 as f64 / agg.2 + 0.5,
                ))
                .unwrap()
                .into()
        })
        .collect();

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::matrix::*;

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

    #[test]
    fn adjacent() {
        let centre = Hit::new(0, (2, 2));

        assert!(Hit::new(0, (1, 1)).adjacent(&centre));
        assert!(Hit::new(0, (1, 2)).adjacent(&centre));
        assert!(Hit::new(0, (1, 3)).adjacent(&centre));
        assert!(Hit::new(0, (2, 1)).adjacent(&centre));
        assert!(Hit::new(0, (2, 3)).adjacent(&centre));
        assert!(Hit::new(0, (3, 1)).adjacent(&centre));
        assert!(Hit::new(0, (3, 2)).adjacent(&centre));
        assert!(Hit::new(0, (3, 3)).adjacent(&centre));

        assert!(Hit::new(0, (2, 2)).adjacent(&centre));

        assert!(!Hit::new(1, (0, 2)).adjacent(&centre));
        assert!(!Hit::new(0, (2, 0)).adjacent(&centre));
        assert!(!Hit::new(0, (2, 4)).adjacent(&centre));
        assert!(!Hit::new(0, (4, 2)).adjacent(&centre));
    }

    #[test]
    fn ccl_() {
        let hits = vec![
            Hit::new(0, (2, 2)),
            Hit::new(0, (3, 3)),
            Hit::new(0, (9, 9)),
            Hit::new(1, (3, 3)),
            Hit::new(1, (2, 2)),
            Hit::new(0, (5, 5)),
            Hit::new(0, (4, 4)),
        ];

        assert_eq!(ccl(&hits), vec![0, 0, 2, 3, 3, 0, 0]);
    }

    #[test]
    fn clustering_() {
        let module0 = DetectorModule::new(
            (10., 10.),
            (10, 10),
            Vector3::default(),
            Matrix3::identity(),
        )
        .unwrap();
        let module1 = DetectorModule::new(
            (20., 20.),
            (10, 10),
            Vector3::new(1., 2., 3.),
            Matrix3::from_angles(std::f64::consts::PI / 2., -std::f64::consts::PI / 2., 0.)
                .unwrap(),
        )
        .unwrap();
        let mut modules = HashMap::new();
        modules.insert(0, module0);
        modules.insert(1, module1);
        let hits: Vec<Hit> = vec![
            Hit::new(0, (4, 4)),
            Hit::new(0, (0, 0)),
            Hit::new(0, (5, 5)),
            Hit::new(1, (5, 5)),
        ];
        let result: Vec<SpacePoint> = clustering(&modules, &hits).unwrap();

        let expected = [[5., 5., 0.], [0.5, 0.5, 0.], [-10., 2., 14.]];

        println!("{:?}", result);
        for i in 0..result.len() {
            for j in 0..3 {
                assert!((result[i].0[j] - expected[i][j]).abs() <= 10. * std::f64::EPSILON);
            }
        }
    }
}

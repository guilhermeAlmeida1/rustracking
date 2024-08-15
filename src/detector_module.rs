use crate::matrix::*;

type Pos2 = (u64, u64);

#[derive(Debug)]
pub struct DetectorModule {
    id: u64,
    dims: (f64, f64),
    pixels_dims: Pos2,
    hits: Vec<Pos2>,
    translation: Vector3<f64>,
    rotation: Matrix3<f64>,
}

impl DetectorModule {
    pub fn new(
        id: u64,
        dims: (f64, f64),
        pixels_dims: Pos2,
        translation: Vector3<f64>,
        rotation: Matrix3<f64>,
    ) -> Result<Self, &'static str> {
        let inv = rotation.inverse()?;
        let transp = rotation.transpose();
        for i in 0..9 {
            if (inv[i] - transp[i]).abs() >= std::f64::EPSILON {
                return Err("Rotation matrix must be orthogonal.");
            }
        }

        Ok(DetectorModule {
            id,
            dims,
            pixels_dims,
            hits : vec![],
            translation,
            rotation,
        })
    }

    pub fn vertices(&self) -> [Vector3<f64>; 4] {
        let v1: Vector3<f64> = self.translation;
        let vec_x = Vector3::<f64> {
            data: [self.dims.0, 0., 0.],
        };
        let v2: Vector3<f64> = self.translation + self.rotation * vec_x;
        let vec_y = Vector3::<f64> {
            data: [0., self.dims.1, 0.],
        };
        let v4: Vector3<f64> = self.translation + self.rotation * vec_y;
        let v3: Vector3<f64> = v4 + v2 - v1;
        [v1, v2, v3, v4]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ctor_invalid_det_0() {
        let transl = [1., 2., 3.].into_vector3().unwrap();
        let rot = [1., 0., 0., 0., 0.5, 0.5, 0., 0.5, 0.5]
            .into_matrix3()
            .unwrap();
        let a = DetectorModule::new(0, (1., 2.), (10, 20), transl, rot);
        assert_eq!(
            a.unwrap_err(),
            "Matrix is not invertible. Determinant is 0."
        );
    }

    #[test]
    fn ctor_invalid_not_orthogonal() {
        let transl = [1., 2., 3.].into_vector3().unwrap();
        let rot = [1., 0., 0., 0., 0.5, 0.5, 1., 2., 0.5]
            .into_matrix3()
            .unwrap();
        let a = DetectorModule::new(0, (1., 2.), (10, 20), transl, rot);
        assert_eq!(a.unwrap_err(), "Rotation matrix must be orthogonal.");
    }

    #[test]
    fn ctor_valid() {
        let transl = [1., 2., 3.].into_vector3().unwrap();
        let rot = [1., 0., 0., 0., 0., -1., 0., 1., 0.]
            .into_matrix3()
            .unwrap();
        let _a = DetectorModule::new(0, (1., 2.), (10, 20), transl, rot).unwrap();

        let rot2 = Matrix3::from_angles(std::f64::consts::PI / 2., 0., 0.).unwrap();
        let _b = DetectorModule::new(0, (0., 0.), (0, 0), transl, rot2).unwrap();
    }

    #[test]
    fn vertices() {
        let dims = (10., 20.);
        let transl = [1., 2., 3.].into_vector3().unwrap();
        let rot = Matrix3::from_angles(std::f64::consts::PI / 2., -std::f64::consts::PI / 2., 0.)
            .unwrap();
        let verts = DetectorModule::new(0, dims, (1, 0), transl, rot)
            .unwrap()
            .vertices();

        println!("vertices: {:?}", verts);
        let expected = [[1., 2., 3.], [1., 2., 13.], [-19., 2., 13.], [-19., 2., 3.]];
        for i in 0..4 {
            for j in 0..3 {
                assert!((verts[i][j] - expected[i][j]).abs() <= 10. * std::f64::EPSILON);
            }
        }
    }
}

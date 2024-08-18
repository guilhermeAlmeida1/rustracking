use crate::matrix::*;

pub type PixelPosition = (u64, u64);

#[derive(Debug)]
pub struct DetectorModule {
    id: u64,
    dims: (f64, f64),
    pixels_dims: PixelPosition,
    translation: Vector3<f64>,
    rotation: Matrix3<f64>,
}

impl DetectorModule {
    pub fn new(
        id: u64,
        dims: (f64, f64),
        pixels_dims: PixelPosition,
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

    pub fn pixel_position_to_vector3(&self, pos: (f64, f64)) -> Result<Vector3<f64>, &'static str> {
        let pixels_dims = (self.pixels_dims.0 as f64, self.pixels_dims.1 as f64);
        if pos.0 > pixels_dims.0 || pos.1 > pixels_dims.1 {
            return Err("Pixel position larger than dimensions of the module.");
        }

        let pos = (
            pos.0 as f64 / self.pixels_dims.0 as f64 * self.dims.0,
            pos.1 as f64 / self.pixels_dims.1 as f64 * self.dims.1,
            0.0,
        )
            .into_vector3()?;
        let pos = self.translation + self.rotation * pos;
        Ok(pos)
    }

    pub fn pixel_vertices(&self, pos: PixelPosition) -> Result<[Vector3<f64>; 4], &'static str> {
        if pos.0 > self.pixels_dims.0 || pos.1 > self.pixels_dims.1 {
            return Err("Pixel position larger than dimensions of the module.");
        }
        let step_x = self.rotation
            * [self.dims.0 / self.pixels_dims.0 as f64, 0., 0.]
                .into_vector3()
                .unwrap();
        let step_y = self.rotation
            * [0., self.dims.1 / self.pixels_dims.1 as f64, 0.]
                .into_vector3()
                .unwrap();

        let v1: Vector3<f64> = self.translation
            + self.rotation
                * [
                    pos.0 as f64 * self.dims.0 / self.pixels_dims.0 as f64,
                    pos.1 as f64 * self.dims.1 / self.pixels_dims.1 as f64,
                    0.,
                ]
                .into_vector3()
                .unwrap();
        let v2: Vector3<f64> = v1 + step_x;
        let v3: Vector3<f64> = v2 + step_y;
        let v4: Vector3<f64> = v3 - step_x;
        Ok([v1, v2, v3, v4])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ctor_invalid_det_0() {
        let rot = [1., 0., 0., 0., 0.5, 0.5, 0., 0.5, 0.5]
            .into_matrix3()
            .unwrap();
        let a = DetectorModule::new(0, (0., 0.), (0, 0), Vector3::identity(), rot);
        assert_eq!(
            a.unwrap_err(),
            "Matrix is not invertible. Determinant is 0."
        );
    }

    #[test]
    fn ctor_invalid_not_orthogonal() {
        let rot = [1., 0., 0., 0., 0.5, 0.5, 1., 2., 0.5]
            .into_matrix3()
            .unwrap();
        let a = DetectorModule::new(0, (0., 0.), (0, 0), Vector3::identity(), rot);
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

    #[test]
    fn pixel_pos_to_vector3() {
        let dims = (10., 20.);
        let transl = [1., 2., 3.].into_vector3().unwrap();
        let rot = Matrix3::from_angles(std::f64::consts::PI / 2., -std::f64::consts::PI / 2., 0.)
            .unwrap();
        let module = DetectorModule::new(0, dims, (10, 10), transl, rot).unwrap();

        let points: Vec<_> = [(0., 0.), (10., 0.), (0., 10.), (10., 10.), (5., 5.)]
            .iter()
            .map(|it| module.pixel_position_to_vector3(*it).unwrap())
            .collect();

        let expected = [
            [1., 2., 3.],
            [1., 2., 13.],
            [-19., 2., 3.],
            [-19., 2., 13.],
            [-9., 2., 8.],
        ];

        println!("points: {:?}", points);
        for i in 0..points.len() {
            for j in 0..3 {
                assert!((points[i][j] - expected[i][j]).abs() <= 10. * std::f64::EPSILON);
            }
        }
    }

    #[test]
    fn pixel_pos_to_vector3_invalid() {
        let module = DetectorModule::new(
            0,
            (0., 0.),
            (10, 10),
            Vector3::identity(),
            Matrix3::identity(),
        )
        .unwrap();

        assert_eq!(
            module.pixel_position_to_vector3((0., 11.)).unwrap_err(),
            "Pixel position larger than dimensions of the module."
        );
    }

    #[test]
    fn pixel_vertices() {
        let dims = (10., 20.);
        let transl = [1., 2., 3.].into_vector3().unwrap();
        let rot = Matrix3::from_angles(std::f64::consts::PI / 2., -std::f64::consts::PI / 2., 0.)
            .unwrap();
        let module = DetectorModule::new(0, dims, (10, 10), transl, rot).unwrap();

        let expected = [[1., 2., 3.], [1., 2., 4.], [-1., 2., 4.], [-1., 2., 3.]];
        let vertices = module.pixel_vertices((0, 0)).unwrap();
        println!("vertices: {:?}", vertices);

        for i in 0..vertices.len() {
            for j in 0..3 {
                assert!((vertices[i][j] - expected[i][j]).abs() <= 10. * std::f64::EPSILON);
            }
        }

        let expected1 = [
            [-17., 2., 12.],
            [-17., 2., 13.],
            [-19., 2., 13.],
            [-19., 2., 12.],
        ];
        let vertices1 = module.pixel_vertices((9, 9)).unwrap();
        println!("vertices1: {:?}", vertices1);

        for i in 0..vertices1.len() {
            for j in 0..3 {
                assert!((vertices1[i][j] - expected1[i][j]).abs() <= 10. * std::f64::EPSILON);
            }
        }
    }

    #[test]
    fn pixel_vertices_invalid() {
        let module = DetectorModule::new(
            0,
            (0., 0.),
            (10, 10),
            Vector3::identity(),
            Matrix3::identity(),
        )
        .unwrap();

        assert_eq!(
            module.pixel_vertices((0, 11)).unwrap_err(),
            "Pixel position larger than dimensions of the module."
        );
    }
}

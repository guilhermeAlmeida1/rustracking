use crate::matrix::{Matrix3, IntoMatrix3, Vector3};

type Pos2 = (i32, i32);

#[derive(Debug)]
pub struct DetectorModule {
    dims: (f64, f64),
    pixels_dims: Pos2,
    cells: Vec<Pos2>,
    translation: Vector3<f64>,
    rotation: Matrix3<f64>,
}

impl DetectorModule {
    pub fn new(
        dims: (f64, f64),
        pixels_dims: Pos2,
        cells: Vec<Pos2>,
        translation: Vector3<f64>,
        rotation: Matrix3<f64>,
    ) -> Result<Self, &'static str> {

        match rotation.inverse() {
            Ok(inv) => if inv != rotation.transpose() {
                return Err("Rotation matrix must be orthogonal.");
            }
            Err(err) => {
                return Err(err);
            }
        }

        Ok(DetectorModule {
            dims, pixels_dims, cells, translation, rotation,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ctor_invalid_det_0() {
        let transl = [1.,2.,3.];
        let rot = [1.,0.,0.,0., 0.5, 0.5, 0., 0.5, 0.5].into_matrix3().unwrap();
        let a = DetectorModule::new((1., 2.), (10,20), vec![(0,0)],
        transl, rot);
        assert_eq!(a.unwrap_err(), "Matrix is not invertible. Determinant is 0.");
    }

    #[test]
    fn ctor_invalid_not_orthogonal() {
        let transl = [1.,2.,3.];
        let rot = [1.,0.,0.,0., 0.5, 0.5, 1., 2., 0.5].into_matrix3().unwrap();
        let a = DetectorModule::new((1., 2.), (10,20), vec![(0,0)],
        transl, rot);
        assert_eq!(a.unwrap_err(), "Rotation matrix must be orthogonal.");
    }

    #[test]
    fn ctor_valid() {
        let transl = [1.,2.,3.];
        let rot = [1.,0.,0.,0., 0., -1., 0., 1., 0.].into_matrix3().unwrap();
        let _a = DetectorModule::new((1., 2.), (10,20), vec![(0,0)],
        transl, rot).unwrap();
    }

}

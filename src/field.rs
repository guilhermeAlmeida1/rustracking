use crate::matrix::Vector3;

pub trait Field {
    fn at(&self, pos: Vector3<f64>) -> Vector3<f64>;
}

pub struct ConstantField {
    pub field: Vector3<f64>,
}
impl Field for ConstantField {
    fn at(&self, _: Vector3<f64>) -> Vector3<f64> {
        self.field
    }
}
impl ConstantField {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self {
            field: Vector3::new(x, y, z),
        }
    }
}

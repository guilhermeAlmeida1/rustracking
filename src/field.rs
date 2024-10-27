pub trait Field {
    fn at(&self, pos: [f64; 3]) -> [f64; 3];
}

pub struct ConstantField {
    pub field: [f64; 3],
}
impl Field for ConstantField {
    fn at(&self, _: [f64; 3]) -> [f64; 3] {
        self.field
    }
}
impl ConstantField {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { field: [x, y, z] }
    }
}

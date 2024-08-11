use std::f64::consts::PI;
use std::slice::SliceIndex;
use std::{default::Default, fmt::Debug, ops};

pub trait Good<T>:
    Default
    + Copy
    + ops::Add<Output = T>
    + ops::AddAssign
    + ops::Sub<Output = T>
    + ops::SubAssign
    + ops::Neg<Output = T>
    + ops::Mul<Output = T>
    + ops::MulAssign
    + ops::Div<Output = T>
    + ops::DivAssign
    + From<i8>
    + PartialEq
    + std::fmt::Debug
{
}
impl<T> Good<T> for T where
    T: Default
        + Copy
        + ops::Add<T, Output = T>
        + ops::AddAssign<T>
        + ops::Sub<T, Output = T>
        + ops::SubAssign<T>
        + ops::Neg<Output = T>
        + ops::Mul<T, Output = T>
        + ops::MulAssign<T>
        + ops::Div<T, Output = T>
        + ops::DivAssign<T>
        + From<i8>
        + PartialEq
        + std::fmt::Debug
{
}

#[derive(Debug, PartialEq, Clone, Copy, Default)]
pub struct Vector3<T: Good<T>> {
    pub data: [T; 3],
}

impl Vector3<f64> {
    pub fn distance(&self, other: &Vector3<f64>) -> f64 {
        f64::sqrt(
            self.data
                .iter()
                .zip(other.data.iter())
                .map(|(a, b)| (a - b) * (a - b))
                .sum(),
        )
    }
}

impl<T: Good<T>> ops::Add<Vector3<T>> for Vector3<T> {
    type Output = Vector3<T>;
    fn add(self, rhs: Vector3<T>) -> Self::Output {
        let mut out: Self::Output = Default::default();
        for i in 0..3 {
            out[i] = self[i] + rhs[i];
        }
        out
    }
}

impl<T: Good<T>> ops::Sub<Vector3<T>> for Vector3<T> {
    type Output = Vector3<T>;
    fn sub(self, rhs: Vector3<T>) -> Self::Output {
        let mut out: Self::Output = Default::default();
        for i in 0..3 {
            out[i] = self[i] - rhs[i];
        }
        out
    }
}

impl<T: Good<T>, I: SliceIndex<[T]>> ops::Index<I> for Vector3<T> {
    type Output = I::Output;
    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        std::ops::Index::index(&self.data, index)
    }
}

impl<T: Good<T>, I: SliceIndex<[T]>> ops::IndexMut<I> for Vector3<T> {
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        std::ops::IndexMut::index_mut(&mut self.data, index)
    }
}

pub trait IntoVector3<T: Good<T>> {
    fn into_vector3(&self) -> Result<Vector3<T>, &'static str>;
}

impl<T: Good<T>> IntoVector3<T> for [T; 3] {
    fn into_vector3(&self) -> Result<Vector3<T>, &'static str> {
        Ok(Vector3 { data: self.clone() })
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Matrix3<T: Good<T>> {
    pub data: [T; 9],
}

impl<T: Good<T>> Matrix3<T> {
    pub fn new() -> Self {
        let data: [T; 9] = Default::default();
        Self { data }
    }
    pub fn data(&self) -> [T; 9] {
        self.data
    }
    pub fn set(&mut self, value: T, i: usize, j: usize) {
        assert!(i < 3);
        assert!(j < 3);
        self.data[i * 3 + j] = value;
    }
    pub fn set_global(&mut self, value: T, i: usize) {
        assert!(i < 9);
        self.data[i] = value;
    }
    pub fn get(&self, i: usize, j: usize) -> T {
        assert!(i < 3);
        assert!(j < 3);
        self.data[i * 3 + j]
    }
    pub fn get_global(&self, i: usize) -> T {
        assert!(i < 9);
        self.data[i]
    }
    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut T {
        assert!(i < 3);
        assert!(j < 3);
        &mut self.data[i * 3 + j]
    }

    pub fn row(&self, i: usize) -> [T; 3] {
        assert!(i < 3);
        let g = i * 3;
        [self.data[g], self.data[g + 1], self.data[g + 2]]
    }
    pub fn col(&self, j: usize) -> [T; 3] {
        assert!(j < 3);
        [self.data[j], self.data[j + 3], self.data[j + 6]]
    }
}

impl<T: Good<T>> Matrix3<T> {
    pub fn determinant(&self) -> T {
        self.get(0, 0) * self.get(1, 1) * self.get(2, 2)
            + self.get(0, 1) * self.get(1, 2) * self.get(2, 0)
            + self.get(0, 2) * self.get(1, 0) * self.get(2, 1)
            - self.get(0, 2) * self.get(1, 1) * self.get(2, 0)
            - self.get(0, 1) * self.get(1, 0) * self.get(2, 2)
            - self.get(0, 0) * self.get(1, 2) * self.get(2, 1)
    }

    pub fn inverse(&self) -> Result<Matrix3<T>, &'static str> {
        let det = self.determinant();
        if det == 0i8.into() {
            return Err("Matrix is not invertible. Determinant is 0.");
        }
        let mut r = Matrix3::<T>::new();
        for i in 0..3 {
            for j in 0..3 {
                let m2by2: Vec<T> = self
                    .into_iter()
                    .enumerate()
                    .filter(|&(idx, _)| idx / 3 != i && idx % 3 != j)
                    .map(|(_, item)| *item)
                    .collect();
                let unit: T = 1i8.into();
                let mut result: T = unit / det * ((m2by2[0] * m2by2[3]) - (m2by2[1] * m2by2[2]));
                if (i + j) % 2 != 0 {
                    result = -result;
                }
                // Note: j, i swapped equivalent to transpose
                *r.get_mut(j, i) = result;
            }
        }
        Ok(r)
    }

    pub fn transpose(&self) -> Matrix3<T> {
        let mut r = Matrix3::<T>::new();
        for i in 0..3 {
            for j in 0..3 {
                *r.get_mut(i, j) = self.get(j, i);
            }
        }
        r
    }
}

impl<T: Good<T>, I: SliceIndex<[T]>> ops::Index<I> for Matrix3<T> {
    type Output = I::Output;
    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        std::ops::Index::index(&self.data, index)
    }
}

impl<T: Good<T>, I: SliceIndex<[T]>> ops::IndexMut<I> for Matrix3<T> {
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        std::ops::IndexMut::index_mut(&mut self.data, index)
    }
}

impl<T: Good<T>> IntoIterator for Matrix3<T> {
    type Item = T;
    type IntoIter = std::array::IntoIter<T, 9>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}
impl<'a, T: Good<T>> IntoIterator for &'a Matrix3<T> {
    type Item = &'a T;
    type IntoIter = std::array::IntoIter<&'a T, 9>;

    fn into_iter(self) -> Self::IntoIter {
        let arr: [_; 9] = [
            &self.data[0],
            &self.data[1],
            &self.data[2],
            &self.data[3],
            &self.data[4],
            &self.data[5],
            &self.data[6],
            &self.data[7],
            &self.data[8],
        ];
        arr.into_iter()
    }
}
// galmeida: fix
// impl<'a, T : Good<T>> IntoIterator for &'a mut Matrix3<T> {
//     type Item = &'a mut T;
//     type IntoIter = std::array::IntoIter<&'a mut T, 9>;

//     fn into_iter(self) -> Self::IntoIter {
//         let mut arr : [&mut T; 9];
//         unsafe {
//             arr = [
//                 &mut self.data[0],
//                 &mut self.data[1],
//                 &mut self.data[2],
//                 &mut self.data[3],
//                 &mut self.data[4],
//                 &mut self.data[5],
//                 &mut self.data[6],
//                 &mut self.data[7],
//                 &mut self.data[8],
//             ];
//         };
//         arr.into_iter()
//     }
// }

impl<T: Good<T>> ops::Mul<Vector3<T>> for Matrix3<T> {
    type Output = Vector3<T>;

    fn mul(self, rhs: Vector3<T>) -> Self::Output {
        let mut r = Vector3::<T>::default();
        for i in 0..3 {
            for j in 0..3 {
                r[i] += self.get(i, j) * rhs[j];
            }
        }
        r
    }
}

impl<T: Good<T>> ops::Mul<Matrix3<T>> for Matrix3<T> {
    type Output = Matrix3<T>;

    fn mul(self, rhs: Matrix3<T>) -> Self::Output {
        let mut r = Matrix3::<T>::new();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    *r.get_mut(i, j) += self.get(i, k) * rhs.get(k, j);
                }
            }
        }
        r
    }
}

impl Matrix3<f64> {
    // Euler angles used.
    // a: alpha, b: beta, g:gamma
    pub fn from_angles(a: f64, b: f64, g: f64) -> Result<Self, &'static str> {
        if !(-2. * PI <= a
            && 2. * PI >= a
            && -2. * PI <= b
            && 2. * PI >= b
            && -2. * PI <= g
            && 2. * PI >= g)
        {
            return Err("Given angles are not within boundaries -2PI, 2PI.");
        }
        let mut data: [f64; 9] = Default::default();
        data[0] = b.cos() * g.cos();
        data[1] = a.sin() * b.sin() * g.cos() - a.cos() * g.sin();
        data[2] = a.cos() * b.sin() * g.cos() + a.sin() * g.sin();
        data[3] = b.cos() * g.sin();
        data[4] = a.sin() * b.sin() * g.sin() + a.cos() * g.cos();
        data[5] = a.cos() * b.sin() * g.sin() - a.sin() * g.cos();
        data[6] = -b.sin();
        data[7] = a.sin() * b.cos();
        data[8] = a.cos() * b.cos();
        Ok(Self { data })
    }
}

pub trait IntoMatrix3<T: Good<T>> {
    fn into_matrix3(&self) -> Result<Matrix3<T>, &'static str>;
}

impl<T: Good<T>> IntoMatrix3<T> for [[T; 3]; 3] {
    fn into_matrix3(&self) -> Result<Matrix3<T>, &'static str> {
        let mut m = Matrix3::<T>::new();
        for i in 0..3 {
            for j in 0..3 {
                m.set(self[i][j], i, j);
            }
        }
        Ok(m)
    }
}

impl<T: Good<T>> IntoMatrix3<T> for [T; 9] {
    fn into_matrix3(&self) -> Result<Matrix3<T>, &'static str> {
        let mut m = Matrix3::<T>::new();
        for i in 0..9 {
            m.set_global(self[i], i);
        }
        Ok(m)
    }
}

impl<T: Good<T>> IntoMatrix3<T> for Vec<Vec<T>> {
    fn into_matrix3(&self) -> Result<Matrix3<T>, &'static str> {
        let mut m = Matrix3::<T>::new();
        if self.len() != 3 {
            return Err("Vector of vectors into_matrix3 must have length 3.");
        }
        if self.iter().any(|v| v.len() != 3) {
            return Err("Vector of vectors into_matrix3 must have inner lengths 3.");
        }
        for i in 0..3 {
            for j in 0..3 {
                m.set(self[i][j], i, j);
            }
        }
        Ok(m)
    }
}

impl<T: Good<T>> IntoMatrix3<T> for Vec<T> {
    fn into_matrix3(&self) -> Result<Matrix3<T>, &'static str> {
        let mut m = Matrix3::<T>::new();
        if self.len() != 9 {
            return Err("Vector into_matrix3 must have length 9.");
        }
        for i in 0..9 {
            m.set_global(self[i], i);
        }
        Ok(m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_set_get() {
        let mut a = Matrix3::new();
        a.set(1, 1, 1);
        a.set_global(-1, 5);
        assert_eq!(a.get_global(4), 1);
        assert_eq!(a.get(1, 2), -1);
        assert_eq!(a.data(), [0, 0, 0, 0, 1, -1, 0, 0, 0]);
    }

    #[test]
    fn get_mut() {
        let mut a = Matrix3::<i32>::new();
        *a.get_mut(0, 0) += 2;
        *a.get_mut(1, 1) += 3;
        assert_eq!(a.data(), [2, 0, 0, 0, 3, 0, 0, 0, 0]);
    }

    #[test]
    fn into_matrix3_invalid_sizes() {
        assert_eq!(
            vec!(vec!(0, 1, 2)).into_matrix3().unwrap_err(),
            "Vector of vectors into_matrix3 must have length 3."
        );
        assert_eq!(
            vec!(vec!(0, 1, 2), vec!(0, 1, 2), vec!(0, 1))
                .into_matrix3()
                .unwrap_err(),
            "Vector of vectors into_matrix3 must have inner lengths 3."
        );
        assert_eq!(
            vec!(0, 1, 2).into_matrix3().unwrap_err(),
            "Vector into_matrix3 must have length 9."
        );
    }

    #[test]
    fn array_into_matrix3() {
        let a = [0, 1, 2, 3, 4, 5, 6, 7, 8].into_matrix3().unwrap();
        assert_eq!(a.data(), [0, 1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn array_array_into_matrix3() {
        let a = [[0, 1, 2], [3, 4, 5], [6, 7, 8]].into_matrix3().unwrap();
        assert_eq!(a.data(), [0, 1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn vec_vec_into_matrix3() {
        let a: Vec<Vec<i32>> = vec![vec![0, 1, 2], vec![3, 4, 5], vec![6, 7, 8]];
        let a = a.into_matrix3().unwrap();
        assert_eq!(a.data(), [0, 1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn vec_into_matrix3() {
        let a: Vec<i32> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8];
        let a = a.into_matrix3().unwrap();
        assert_eq!(a.data(), [0, 1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn row_col() {
        let a = [[0, 1, 2], [3, 4, 5], [6, 7, 8]].into_matrix3().unwrap();
        for i in 0..3 {
            assert_eq!(a.row(i), [a.get(i, 0), a.get(i, 1), a.get(i, 2)]);
        }
        for j in 0..3 {
            assert_eq!(a.col(j), [a.get(0, j), a.get(1, j), a.get(2, j)]);
        }
    }

    #[test]
    fn determinant() {
        let a = [[0, 1, 2], [3, 4, 5], [6, 7, 8]].into_matrix3().unwrap();
        assert_eq!(a.determinant(), 0);
    }

    #[test]
    fn inverse_invalid_det0() {
        let a = [[0, 1, 2], [3, 4, 5], [6, 7, 8]].into_matrix3().unwrap();
        assert_eq!(
            a.inverse().unwrap_err(),
            "Matrix is not invertible. Determinant is 0."
        );
    }

    #[test]
    fn inverse_invalid_det0_float() {
        let a = [[0., 1., 2.], [3., 4., 5.], [6., 7., 8.]]
            .into_matrix3()
            .unwrap();
        assert_eq!(
            a.inverse().unwrap_err(),
            "Matrix is not invertible. Determinant is 0."
        );
    }

    #[test]
    fn inverse() {
        let a: Matrix3<i32> = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].into_matrix3().unwrap();
        assert_eq!(a.inverse().unwrap().data(), [1, 0, 0, 0, 1, 0, 0, 0, 1]);
        let a: Matrix3<f32> = [[1., 2., 3.], [0., 1., 4.], [5., 6., 0.]]
            .into_matrix3()
            .unwrap();
        assert_eq!(
            a.inverse().unwrap().data(),
            [-24., 18., 5., 20., -15., -4., -5., 4., 1.]
        );
    }

    #[test]
    fn transpose() {
        let a = [[0, 1, 2], [3, 4, 5], [6, 7, 8]].into_matrix3().unwrap();
        assert_eq!(a.transpose().data(), [0, 3, 6, 1, 4, 7, 2, 5, 8]);
    }

    #[test]
    fn iter_imut() {
        let a = [[0, 1, 2], [3, 4, 5], [6, 7, 8]].into_matrix3().unwrap();
        let b: Vec<_> = a.into_iter().map(|x| x).collect();
        assert_eq!(b, vec![0, 1, 2, 3, 4, 5, 6, 7, 8]);
    }

    // #[test]
    // fn iter_mut() {
    //     let mut a = [[0, 1, 2], [3, 4, 5], [6, 7, 8]].into_matrix3().unwrap();
    //     a.into_iter().for_each(|x| *x += 1);
    //     assert_eq!(a.data(), [1,2,3,4,5,6,7,8,9]);
    // }

    #[test]
    fn matrix_multiplication() {
        let a = [[0, 1, 2], [3, 4, 5], [6, 7, 8]].into_matrix3().unwrap();
        let b = [[8, 7, 6], [5, 4, 3], [2, 1, 0]].into_matrix3().unwrap();
        let c = a * b;
        assert_eq!(c.data(), [9, 6, 3, 54, 42, 30, 99, 78, 57]);
    }

    #[test]
    fn matrix_multiplication_vector() {
        let a: Matrix3<i32> = [[0, 1, 2], [3, 4, 5], [6, 7, 8]].into_matrix3().unwrap();
        let b: Vector3<i32> = [0, 1, 2].into_vector3().unwrap();
        let c = a * b;
        assert_eq!(c.data, [5, 14, 23]);
    }

    #[test]
    fn from_angles() {
        let a = Matrix3::from_angles(PI / 2., 0., 0.).unwrap();
        let b = Matrix3::from_angles(0., PI / 2., 0.).unwrap();
        let c = Matrix3::from_angles(0., 0., PI / 2.).unwrap();
        let expected_a = [1., 0., 0., 0., 0., -1., 0., 1., 0.];
        let expected_b = [0., 0., 1., 0., 1., 0., -1., 0., 0.];
        let expected_c = [0., -1., 0., 1., 0., 0., 0., 0., 1.];
        for i in 0..9 {
            assert!((a.data[i] - expected_a[i]).abs() < std::f64::EPSILON);
            assert!((b.data[i] - expected_b[i]).abs() < std::f64::EPSILON);
            assert!((c.data[i] - expected_c[i]).abs() < std::f64::EPSILON);
        }
    }

    #[test]
    fn vec3f64_distance() {
        let a = [1., 2., 5.].into_vector3().unwrap();
        let b = [4., 2., 1.].into_vector3().unwrap();
        println!("distance = {}", a.distance(&b));
        assert!((a.distance(&b) - 5.).abs() < std::f64::EPSILON);
    }
}

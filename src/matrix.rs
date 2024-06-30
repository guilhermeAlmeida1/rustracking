use std::{default::Default, ops};

#[derive(Default, Debug)]
pub struct Vector3<T: Default>([T; 3]);

#[derive(Default, Debug)]
pub struct Matrix3<T>
where
    T: Default,
{
    data: [T; 9],
}
impl<T: Default + Copy> Matrix3<T> {
    pub fn new() -> Self {
        Default::default()
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

    // galmeida: missing rows() and cols() iters
}

impl<T: Default> IntoIterator for Matrix3<T> {
    type Item = T;
    type IntoIter = std::array::IntoIter<T, 9>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}
impl<'a, T: Default> IntoIterator for &'a Matrix3<T> {
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
// impl<'a, T: Default> IntoIterator for &'a mut Matrix3<T> {
//     type Item = &'a mut T;
//     type IntoIter = std::array::IntoIter<&'a mut T, 9>;

//     fn into_iter(self) -> Self::IntoIter {
//         unsafe {
//             let arr: [_; 9] = [
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
//             arr.into_iter()
//         }
//     }
// }

impl<T: Default + Copy + ops::AddAssign + ops::Mul<Output = T>> ops::Mul<Matrix3<T>>
    for Vector3<T>
{
    type Output = Vector3<T>;

    fn mul(self, rhs: Matrix3<T>) -> Self::Output {
        let mut r: Self::Output = Self::Output::default();
        for i in 0..3 {
            for j in 0..3 {
                r.0[i] += self.0[j] * rhs.get(j, i);
            }
        }
        r
    }
}

impl<T: Default + Copy + ops::AddAssign + ops::Mul<Output = T>> ops::Mul<Matrix3<T>>
    for Matrix3<T>
{
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

pub trait IntoMatrix3<T: Copy + Default> {
    fn into_matrix3(&self) -> Result<Matrix3<T>, &'static str>;
}

impl<T: Default + Copy> IntoMatrix3<T> for [[T; 3]; 3] {
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

impl<T: Default + Copy> IntoMatrix3<T> for [T; 9] {
    fn into_matrix3(&self) -> Result<Matrix3<T>, &'static str> {
        let mut m = Matrix3::<T>::new();
        for i in 0..9 {
            m.set_global(self[i], i);
        }
        Ok(m)
    }
}

impl<T: Default + Copy> IntoMatrix3<T> for Vec<Vec<T>> {
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

impl<T: Default + Copy> IntoMatrix3<T> for Vec<T> {
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
        let a: Vec<Vec<usize>> = vec![vec![0, 1, 2], vec![3, 4, 5], vec![6, 7, 8]];
        let a = a.into_matrix3().unwrap();
        assert_eq!(a.data(), [0, 1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn vec_into_matrix3() {
        let a: Vec<usize> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8];
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
    fn iter_imut() {
        let a = [[0, 1, 2], [3, 4, 5], [6, 7, 8]].into_matrix3().unwrap();
        let b: Vec<_> = a.into_iter().map(|x|{
            x
        }).collect();
        assert_eq!(b, vec![0,1,2,3,4,5,6,7,8]);
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
        let b: Vector3<i32> = Vector3([0, 1, 2]);
        let c = b * a;
        assert_eq!(c.0, [15, 18, 21]);
    }
}

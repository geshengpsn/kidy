use liealg::{Algebra, Real, Vec6, Vector, SE3};
use nalgebra::SMatrix;

pub struct Jacobian<T, const N: usize> {
    val: SMatrix<T, 6, N>,
}

impl<T: Real, const N: usize> Jacobian<T, N> {
    fn new() -> Self {
        Self {
            val: SMatrix::zeros(),
        }
    }

    fn set_column(&mut self, idx: usize, s: &Vec6<T>) {
        self.val
            .view_mut((0, idx), (6, 1))
            .copy_from_slice(&s.as_array());
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum JocbError {}

pub fn jacobian<T: Real, const N: usize>(
    s_list: &[Vec6<T>; N],
    theta_list: &[T; N],
) -> Result<Jacobian<T, N>, JocbError> {
    // let mut j = Jacobian::new();
    // j.set_column(0, &s_list[0]);
    // let mut trans = SE3::identity();
    // for i in 0..N - 1 {
    //     trans *= (s_list[i] * theta_list[i]).hat().exp();
    //     j.view_mut((0, i + 1), (6, 1))
    //         .copy_from(&(adjoint(&trans) * s_list[i + 1]));
    // }
    // Ok(j)
    todo!()
}

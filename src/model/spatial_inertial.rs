use nalgebra::{Matrix3, Matrix4, Matrix6};
use spatial_alg::prelude::*;

pub(super) fn to_local_spatial_inertial(
    inertia_frame: &Matrix4<f64>,
    inertia: &Matrix3<f64>,
    mass: f64,
) -> Matrix6<f64> {
    let mut i_b = Matrix6::from_diagonal_element(mass);
    i_b.fixed_view_mut::<3, 3>(0, 0).copy_from(inertia);
    let b_t_a = inertia_frame.inv();
    let adj_b_t_a = Matrix6::from_column_slice(b_t_a.lee_adjoint().as_slice());
    adj_b_t_a.transpose() * i_b * adj_b_t_a
}

#[cfg(test)]
mod tests {
    use nalgebra::Matrix4;

    #[test]
    fn to_local_spatial_inertial_test() {
        use nalgebra::Matrix3;

        use crate::model::spatial_inertial::to_local_spatial_inertial;

        let mut inertia_frame = Matrix4::identity();
        inertia_frame[(0, 3)] = 1.0; // translation along x-axis
        let inertia = Matrix3::from_diagonal_element(4.);
        let mass = 5.0;
        let spatial_inertial = to_local_spatial_inertial(&inertia_frame, &inertia, mass);
        println!("spatial_inertial: {:.3}", spatial_inertial);
    }
}

// pub mod base_functions;

// use std::simd::StdFloat;
use liealg::prelude::*;
use liealg::Vec6;
use liealg::SE3;
use nalgebra::SMatrix;
use nalgebra::SVector;
use nalgebra::Vector3;
use nalgebra::Vector6;

pub(crate) struct KiChain<const N:usize> {
    pub(crate) joints_screw: [Vec6<f64>; N],
    pub(crate) zero: SE3<f64>,
}

#[derive(Debug)]
pub(crate) enum IkError {
    PseudoInverse,
    MaxIter,
}

impl<const N:usize> KiChain<N> {
    pub(crate) fn fk(&self, joints: &[f64; N]) -> SE3<f64> {
        let mut pose = SE3::<f64>::identity();
        for (joint, screw) in joints.iter().zip(self.joints_screw.iter()) {
            pose *= (screw.hat() * joint).exp();
        }
        pose * &self.zero
    }

    pub(crate) fn jacobian(&self, joints: &[f64]) -> SMatrix<f64, 6, N> {
        let mut jacobian = SMatrix::<f64, 6, N>::zeros();
        let mut pose = SE3::<f64>::identity();
        for (i, (joint, screw)) in joints.iter().zip(self.joints_screw.iter()).enumerate() {
            let jn = screw;
            let jn = pose.adjoint().act(&jn.hat()).vee();
            jacobian
                .fixed_view_mut::<6, 1>(0, i)
                .copy_from_slice(jn.as_slice());
            pose *= (screw.hat() * joint).exp();
        }
        jacobian
    }

    pub(crate) fn ik(&self, target_pose: &SE3<f64>, init_joints: &[f64; N], r_error: f64, p_error: f64, pinv_eps: f64, max_time: usize) -> Result<([f64; N],usize), IkError>
    where
        nalgebra::Const<6>: nalgebra::DimMin<nalgebra::Const<N>>,
        nalgebra::Const<N>: nalgebra::ToTypenum,
        <nalgebra::Const<6> as nalgebra::DimMin<nalgebra::Const<N>>>::Output: nalgebra::DimSub<nalgebra::Const<1>>,
        nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<<<nalgebra::Const<6> as nalgebra::DimMin<nalgebra::Const<N>>>::Output as nalgebra::DimSub<nalgebra::Const<1>>>::Output>,
        nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<<nalgebra::Const<6> as nalgebra::DimMin<nalgebra::Const<N>>>::Output, nalgebra::Const<N>>,
        nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<nalgebra::Const<6>, <nalgebra::Const<6> as nalgebra::DimMin<nalgebra::Const<N>>>::Output>,
        nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<<nalgebra::Const<6> as nalgebra::DimMin<nalgebra::Const<N>>>::Output>,
    {
        let mut current_joints = *init_joints;
        for t in 0..max_time {
            let current_pose = self.fk(&current_joints);
            println!("current_joints: {current_joints:?}");
            println!("current_pose: {current_pose:.7}");
            // T_c^-1 * T_t = c_T_t
            let t_ct = (current_pose.inv() * target_pose).log();
            let vs = current_pose.adjoint().act(&t_ct).vee();
            println!("{}", vs);
            if Vector3::new(vs.r()[0], vs.r()[1], vs.r()[2]).norm() <= r_error
                && Vector3::new(vs.p()[0], vs.p()[1], vs.p()[2]).norm() <= p_error
            {
                return Ok((current_joints, t));
            }
            let Ok(pinv_j) = self.jacobian(&current_joints).pseudo_inverse(pinv_eps) else {
                return Err(IkError::PseudoInverse);
            };
            let v = Vector6::from_row_slice(vs.as_slice());
            let update: SVector<f64, N> = pinv_j * v;
            println!("pinv: {pinv_j:.5} update:{update}");
            for (j, u) in current_joints.iter_mut().zip(update.iter()) {
                *j += *u;
            }
        }
        Err(IkError::MaxIter)
    }
}

#[cfg(test)]
mod test {
    use approx::assert_relative_eq;
    use liealg::SO3;

    use super::*;
    use core::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn test_fk() {
        let chain = KiChain {
            joints_screw: [
                Vec6::new([0., 0., 1.], [0., 0., 0.]),
                Vec6::new([0., 0., 1.], [0., -1., 0.]),
            ],
            zero: (Vec6::new([0., 0., 0.], [1., 0., 0.]) * 2.).hat().exp(),
        };

        let theta_list = [0., FRAC_PI_2];
        let res = chain.fk(&theta_list);
        assert_relative_eq!(
            res,
            SE3::new(&SO3::from_euler_angles(0., 0., FRAC_PI_2), [1., 1., 0.])
        );
    }

    #[test]
    fn test_jacobian() {
        let chain = KiChain {
            joints_screw: [
                Vec6::new([0., 0., 1.], [0., 0.2, 0.2]),
                Vec6::new([1., 0., 0.], [2., 0., 3.]),
                Vec6::new([0., 1., 0.], [0., 2., 1.]),
                Vec6::new([1., 0., 0.], [0.2, 0.3, 0.4]),
            ],
            zero: SE3::identity(),
        };
        let theta_list = [0.2, 1.1, 0.1, 1.2];
        let res = chain.jacobian(&theta_list);
        let expected = SMatrix::<f64, 6, 4>::from_column_slice(&[
            0.0,
            0.0,
            1.0,
            0.0,
            0.2,
            0.2,
            0.9800665778412416,
            0.19866933079506122,
            0.0,
            1.9521863824506809,
            0.4365413247037721,
            2.960266133840988,
            -0.09011563789485476,
            0.4445543984476258,
            0.8912073600614354,
            -2.216352156896298,
            -2.437125727653321,
            3.2357306532803083,
            0.957494264730031,
            0.28487556541794845,
            -0.04528405057966491,
            -0.5116153729819477,
            2.7753571339551537,
            2.2251244335357394,
        ]);
        assert_relative_eq!(res, expected);
    }

    #[test]
    fn test_ik() {
        let chain = KiChain {
            joints_screw: [
                Vec6::new([0., 0., 1.], [4., 0., 0.]),
                Vec6::new([0., 0., 0.], [0., 1., 0.]),
                Vec6::new([0., 0., -1.], [-6., 0., -0.1]),
            ],
            zero: SE3::new(&SO3::from_euler_angles(0., PI, 0.), [0., 6., 2.]),
        };

        let target = SE3::new(
            &SO3::new_unchecked(&[0., 1., 0., 1., 0., 0., 0., 0., -1.]),
            [-5., 4., 1.6858],
        );

        let init_joints = [1.5, 2.5, 4.];
        let (joints, iter_time) = chain.ik(&target, &init_joints, 5e-5, 5e-5, 1e-10, 100).unwrap();
        println!("{joints:?} {iter_time}");

        let res = chain.fk(&joints);
        assert_relative_eq!(res, target, epsilon = 1e-4);
    }
}

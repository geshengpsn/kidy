use liealg::prelude::*;
use liealg::SE3;
use nalgebra::SMatrix;
use nalgebra::SVector;
use nalgebra::Vector3;
use nalgebra::Vector6;
use crate::KidyChain;

#[derive(Debug)]
pub enum IkError {
    PseudoInverse,
    MaxIter,
}

impl<const N: usize> KidyChain<N> {
    pub fn fk(&self, joints: &[f64; N]) -> SE3<f64> {
        let mut pose = SE3::<f64>::identity();
        for (joint, screw) in joints.iter().zip(self.joints_screw.iter()) {
            pose *= (screw.clone() * joint).exp();
        }
        pose * &self.zero_ee_pose()
    }

    pub fn jacobian(&self, joints: &[f64]) -> SMatrix<f64, 6, N> {
        let mut jacobian = SMatrix::<f64, 6, N>::zeros();
        let mut pose = SE3::<f64>::identity();
        for (i, (joint, screw)) in joints.iter().zip(self.joints_screw.iter()).enumerate() {
            let jn = screw;
            let jn = pose.adjoint().act(jn).vee();
            jacobian
                .fixed_view_mut::<6, 1>(0, i)
                .copy_from_slice(jn.as_slice());
            pose *= (screw.clone() * joint).exp();
        }
        jacobian
    }

    pub fn ik(&self, target_pose: &SE3<f64>, init_joints: &[f64; N], r_error: f64, p_error: f64, pinv_eps: f64, max_time: usize) -> Result<([f64; N],usize), IkError>
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
            // println!("current_joints: {current_joints:?}");
            // println!("current_pose: {current_pose:.7}");
            // T_c^-1 * T_t = c_T_t
            let t_ct = (current_pose.inv() * target_pose).log();
            let vs = current_pose.adjoint().act(&t_ct).vee();
            // println!("{}", vs);
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
            // println!("pinv: {pinv_j:.5} update:{update}");
            for (j, u) in current_joints.iter_mut().zip(update.iter()) {
                *j += *u;
            }
        }
        Err(IkError::MaxIter)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::MultiBody;
    use approx::assert_relative_eq;
    use core::f64::consts::FRAC_PI_2;
    use liealg::SO3;

    #[test]
    fn test_fk() {
        let multi_body =
            MultiBody::from_urdf("urdf/franka_description/robots/franka_panda.urdf").unwrap();
        let chain = multi_body.get_kidy_chain::<7>("panda_link0", "panda_hand");
        let theta_list = [FRAC_PI_2, 0., 0., 0., 0., 0., 0.];
        let res = chain.fk(&theta_list);
        println!("{res:.5}");
    }

    #[test]
    fn test_jacobian() {
        let multi_body =
            MultiBody::from_urdf("urdf/franka_description/robots/franka_panda.urdf").unwrap();
        let chain = multi_body.get_kidy_chain::<7>("panda_link0", "panda_hand");
        let theta_list = [FRAC_PI_2, 0., 0., 0., 0., 0., 0.];
        let res = chain.jacobian(&theta_list);
        println!("{res:.5}");
    }

    #[test]
    fn test_ik() {
        let multi_body =
            MultiBody::from_urdf("urdf/franka_description/robots/franka_panda.urdf").unwrap();
        let chain = multi_body.get_kidy_chain::<7>("panda_link0", "panda_hand");

        let target = SE3::new(
            &SO3::new_unchecked(&[1., 0., 0., 0., -1., 0., 0., 0., -1.]),
            [0.5, 0., 0.6],
        );

        let init_joints = [0., 0., 0., FRAC_PI_2, 0., 0., 0.];
        let (joints, iter_time) = chain
            .ik(&target, &init_joints, 1e-6, 1e-6, 1e-10, 1000)
            .unwrap();
        println!("{joints:?} {iter_time}");

        let res = chain.fk(&joints);
        println!("{res:.7}");
        assert_relative_eq!(res, target, epsilon = 1e-4);
    }
}

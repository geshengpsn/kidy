use std::f64::consts::FRAC_PI_2;

use kidy::MultiBody;
use liealg::{SE3, SO3};

fn main() {
    let multi_body =
        MultiBody::from_urdf("urdf/franka_description/robots/franka_panda.urdf").unwrap();
    let chain = multi_body.get_kidy_chain::<7>("panda_link0", "panda_hand");

    let target = SE3::new(
        &SO3::new_unchecked(&[1., 0., 0., 0., -1., 0., 0., 0., -1.]),
        [0.5, 0., 0.6],
    );

    let init_joints = [0., 0., 0., FRAC_PI_2, 0., 0., 0.];
    // println!("{init_joints:?}");
    let (joints, iter_time) = chain
        .ik(&target, &init_joints, 1e-6, 1e-6, 1e-10, 1000)
        .unwrap();
}

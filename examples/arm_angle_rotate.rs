use kidy::{
    analytical::realman75_6f::{ElbowConfig, Realman75AikSolver, SolutionConfig, SphereConfig},
    visual::RerunVisualizer,
    Model,
};
use nalgebra::{IsometryMatrix3, Matrix3, Rotation3, Translation3};
use std::f64::consts::PI;

fn main() {
    let model = Model::from_urdf("./urdf/rm_75_6f_description/urdf/rm_75_6f_description.urdf")
        .expect("Failed to load URDF model");
    let mut viewer = RerunVisualizer::new("arm angle rotate", "./urdf/", &model);
    let realman = Realman75AikSolver::default();
    let target_pose = IsometryMatrix3::from_parts(
        Translation3::new(0.3, 0., 0.6),
        Rotation3::from_matrix(&Matrix3::new(0.5, 0., 0.866, 0., 1., 0., -0.866, 0., 0.5)),
    );
    let target_pose = target_pose.to_homogeneous();

    let nums = 200;
    for i in 0..nums {
        let theta = 2. * PI * i as f64 / 200.;
        if let Ok(joints) = realman.aik(
            &target_pose,
            theta,
            SolutionConfig {
                shoulder: SphereConfig::A,
                elbow: ElbowConfig::Inner,
                wrist: SphereConfig::A,
            },
        ) {
            let joints = [
                0., joints[0], joints[1], joints[2], joints[3], joints[4], joints[5], joints[6],
            ];
            viewer.set_joints(&joints);
            viewer
                .rerun
                .log("/joints/", &rerun::Scalars::new(joints))
                .unwrap();
        }
    }
}

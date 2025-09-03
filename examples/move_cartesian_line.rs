use kidy::{
    analytical::realman75_6f::{
        ElbowConfig, JointsNorm, Realman75AikSolver, SolutionConfig, SphereConfig,
    },
    visual::RerunVisualizer,
    Model,
};
use nalgebra::{matrix, IsometryMatrix3, Rotation3, Translation3, Vector3};

fn main() {
    let model = Model::from_urdf("./urdf/rm_75_6f_description/urdf/rm_75_6f_description.urdf")
        .expect("Failed to load URDF model");
    let mut viewer = RerunVisualizer::new("move cartesian line", "./urdf/", &model);
    let realman = Realman75AikSolver::default();

    let arm_t_ee_rot = matrix![
        0., 0., 1.;
        1., 0., 0.;
        0., 1., 0.
    ];

    let arm_ee_a = Vector3::new(0.379, -0.2, 0.6);
    let arm_ee_b = Vector3::new(0.379, 0.2, 0.6);
    let steps = 100;

    // let mut last_joints = SVector::<f64, 7>::zeros();
    for i in 0..=steps {
        let t = i as f64 / steps as f64;
        let target_pose = IsometryMatrix3::from_parts(
            Translation3::new(
                arm_ee_a.x * (1. - t) + arm_ee_b.x * t,
                arm_ee_a.y * (1. - t) + arm_ee_b.y * t,
                arm_ee_a.z * (1. - t) + arm_ee_b.z * t,
            ),
            Rotation3::from_matrix(&arm_t_ee_rot),
        )
        .to_homogeneous();
        if let Ok(theta) = realman.search_brent(
            &target_pose,
            SolutionConfig {
                shoulder: SphereConfig::A,
                elbow: ElbowConfig::Inner,
                wrist: SphereConfig::A,
            },
            JointsNorm,
        ) {
            let joints = realman
                .aik(
                    &target_pose,
                    theta,
                    SolutionConfig {
                        shoulder: SphereConfig::A,
                        elbow: ElbowConfig::Inner,
                        wrist: SphereConfig::A,
                    },
                )
                .unwrap();
            // last_joints = joints;
            let joints = [
                0., joints[0], joints[1], joints[2], joints[3], joints[4], joints[5], joints[6],
            ];
            viewer.set_joints(&joints);
            viewer
                .rerun
                .log("/joints/", &rerun::Scalars::new(joints))
                .unwrap();
        } else {
            println!("unreachable");
        }
    }
}

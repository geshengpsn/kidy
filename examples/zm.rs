use kidy::{
    visual::RerunVisualizer,
    Model,
};

fn main() {
    let zm_left = "./urdf/ZPFR/urdf/ZPFR.urdf";
    let model = Model::from_urdf(zm_left).expect("Failed to load URDF model");
    let mut viewer = RerunVisualizer::new("zm robot", "./urdf/", &model);
    let joints = [0.; 9];
    // let joints = [0., -2.3977,  1.2343,  0.5737, -1.1258, -0.3082, -0.3993,  -0.58850, 0.];
    // joints[7] = 1.;
    viewer.set_joints(&joints);
    // viewer
    //     .rerun
    //     .log("/joints/", &rerun::Scalars::new(joints))
    //     .unwrap();
}

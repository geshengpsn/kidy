use kidy::{visual::RerunVisualizer, Model};
use rerun::{external::re_log, Boxes3D};

fn main() {
    re_log::setup_logging();
    let zm = "./urdf/Bot2_ZPF_L2Fgripper_R2Fgripper_A_description/arm_left.urdf";
    let model = Model::from_urdf(zm).expect("Failed to load URDF model");
    let mut viewer = RerunVisualizer::new(
        "bot2 left",
        "./urdf/Bot2_ZPF_L2Fgripper_R2Fgripper_A_description/",
        &model,
    );
    // let joints = [0., -0.603467153, 0.581839460, -0.872919501, -0.257218995, -0.568523196, -0.348943927, -1.648817499, 0.];
    // let a = model.fk(&joints);
    // for j in a {
    //     println!("{:.9}", j);
    // }

    // for l in &model.links {
    //     println!("{:.9}", l.space_spatial_twist);
    // }

    viewer.set_joints(&[
        0.,
        0.,
        0.,
        0.,
        -1.57,
        0.,
        0.,
        0.,
        0.,
    ]);
    viewer.rerun.log(
        "cube",
        &Boxes3D::from_centers_and_sizes([(-0.5, 0.5, 0.2)], [(0.5, 0.5, 0.5)]),
    ).unwrap();
}

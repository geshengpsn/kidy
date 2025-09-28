use std::f64;

use kidy::{visual::RerunVisualizer, Model};
fn main() {
    let bot2 = "./urdf/Bot2_ZPF_L2Fgripper_R2Fgripper_A_description/urdf/Bot2_ZPF_L2Fgripper_R2Fgripper_A_description.urdf";
    let model = Model::from_urdf(bot2).expect("Failed to load URDF model");
    for l in &model.links {
        println!("{} {}", l.name, l.space_spatial_twist);
    }
    let viewer = RerunVisualizer::new(
        "bot2 left",
        "./urdf/",
        &model,
    );
}
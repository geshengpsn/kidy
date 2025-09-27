use kidy::{
    visual::RerunVisualizer,
    Model,
};

fn main() {
    let rm_75 = "./urdf/rm_75_6f_description/urdf/rm_75_6f_description.urdf";
    let model = Model::from_urdf(rm_75).expect("Failed to load URDF model");
    let mut viewer = RerunVisualizer::new("rm 75 robot", "./urdf/", &model);
    let joints = [0.; 8];
    viewer.set_joints(&joints);
}

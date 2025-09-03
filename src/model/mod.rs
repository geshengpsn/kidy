// use liealg::{se3, SE3};
use nalgebra::{Matrix4, Matrix6, Vector6};
use petgraph::prelude::DiGraphMap;

mod bfs;
mod mjcf;
mod spatial_inertial;
mod urdf;
pub(crate) mod utils;
// mod macros;
// let model = model!("RM75-6F.urdf");
// 内肘 外肘
// tube 曲线
// let solution_set = model.ik();
//

#[derive(Debug)]
pub struct Joint {
    pub urdf_joint: Option<urdf_rs::Joint>,
}

#[derive(Debug)]
pub struct Link {
    pub name: String,
    pub space_spatial_twist: Vector6<f64>,
    pub local_spatial_twist: Vector6<f64>,

    pub global_zero_pose: Matrix4<f64>,

    // zero pose relative to parent link
    pub parent_zero_pose: Matrix4<f64>,

    // link frame spatial inertia
    pub local_spatial_inertial: Matrix6<f64>,

    // joint
    pub joint: Joint,
    // original urdf link
    pub visual: Vec<urdf_rs::Visual>,
    pub collision: Vec<urdf_rs::Collision>,
}

#[derive(Debug)]
pub struct Model {
    // rigid body tree
    pub links: Vec<Link>,
    pub link_graph: DiGraphMap<usize, ()>,
    pub bfs: Vec<usize>,
}

pub trait ModelBase {
    fn model(&self) -> &Model;
}
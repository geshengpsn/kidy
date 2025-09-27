use kidy::Model;
use nalgebra::{geometry::Rotation3, Matrix3, UnitQuaternion};
use rand::Rng;

fn main() {
    let model = Model::from_urdf("./urdf/rm_75_6f_description/urdf/rm_75_6f_description.urdf")
        .expect("Failed to load URDF model");
    let mut rng = rand::rng();
    
    for _ in 0..5 {
        // gen random joint angles
        let mut joints = vec![0.; 8];
        let mut i = 1;
        for l in model.links.iter() {
            if let Some(urdf_joint) = &l.joint.urdf_joint {
                if urdf_joint.joint_type != urdf_rs::JointType::Fixed {
                    joints[i] = rng.random_range(urdf_joint.limit.lower..urdf_joint.limit.upper);
                    i += 1;
                }
            }
        }

        let a = model.fk(&joints);
        dbg!(joints);
        let mat = a.last().unwrap();
        println!("{}", mat);
        let linear = mat.fixed_view::<3,1>(0, 3);
        println!("position:\n{}, {}, {}", linear[0], linear[1], linear[2]);
        let rot_mat: Matrix3<f64> = mat.fixed_view::<3,3>(0, 0).into();
        let quat = UnitQuaternion::from_matrix(&rot_mat);
        println!("quaternion:\n{}, {}, {}, {}", quat.w, quat.i, quat.j, quat.k);
    }
}

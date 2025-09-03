use nalgebra::{Matrix3, Vector3, VectorView3, U1, U6};
use urdf_rs::JointType;

use crate::{Link, Model};

pub struct SRSModel {
    pub model: Model,
    pub shoulder_euler_order: EulerOrder,
    pub shoulder_angle_sign: [f64; 3], // Sign of the angles for each axis in the shoulder joint, (-1. or 1.)
    pub wrist_euler_order: EulerOrder,
    pub wrist_angle_sign: [f64; 3], // Sign of the angles for each axis in the wrist joint, (-1. or 1.)
    pub upper_arm: f64,
    pub fore_arm: f64,
    pub wrist_point_under_tool_frame: Vector3<f64>,
    pub shoulder_point: Vector3<f64>,
    pub s_r_se: Matrix3<f64>,
    pub w_r_ew: Matrix3<f64>,
    pub joints_limit: [(f64, f64); 7],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EulerAxis {
    X,
    Y,
    Z,
}

#[derive(Debug)]
pub struct EulerOrder([EulerAxis; 3]);

impl EulerOrder {
    // Validate the Euler axis, rule: near axes must be orthogonal
    fn is_validate(&self) -> bool {
        if self.0[0] == self.0[1] || self.0[1] == self.0[2] {
            return false; // No axis should repeat
        }
        true
    }
}

impl Model {
    pub fn try_srs(&self) -> Result<SRSModel, String> {
        // Check if the model is a 7-DOF robot
        if !is_7dof(self) {
            return Err(format!(
                "Model is not a 7-DOF robot, found {} joints",
                self.links.len() - 1
            ));
        }

        // Check if the model's joints are all revolute
        if !is_all_revolute_joints(self) {
            return Err("Model does not have all revolute joints".to_string());
        }

        // extract shoulder spherical joints
        // shoulder point
        // s_r_se

        // extract wrist spherical joints
        // wrist point under tool frame
        // w_r_ew

        // upper arm
        // fore arm

        todo!()
    }
}

fn is_7dof(m: &Model) -> bool {
    // Check if the model has 7 degrees of freedom
    m.links.len() == 8
}

fn is_all_revolute_joints(m: &Model) -> bool {
    // Check if all joints are revolute
    m.links.iter().skip(1).all(|l| {
        l.joint.urdf_joint.as_ref().is_some_and(|j| {
            j.joint_type == JointType::Revolute || j.joint_type == JointType::Continuous
        })
    })
}

fn extract_sphere_joint(
    l1: &Link,
    l2: &Link,
    l3: &Link,
    epsilon: f64,
) -> Result<EulerOrder, String> {
    // Extract the spherical joint from the space_spatial_twist of links
    // For a spherical joint, three consecutive revolute joints should have
    // their rotation axes intersecting at a common point
    let (p1, axis1) = rot_line(l1);
    let (p2, axis2) = rot_line(l2);
    let (p3, axis3) = rot_line(l3);

    // get the order of the axes
    let (a1, a1_sign) = axis_to_euler_axis(&axis1, epsilon)?;
    let (a2, a2_sign) = axis_to_euler_axis(&axis2, epsilon)?;
    let (a3, a3_sign) = axis_to_euler_axis(&axis3, epsilon)?;
    let euler_order = EulerOrder([a1, a2, a3]);
    let signs = [a1_sign, a2_sign, a3_sign];
    if !euler_order.is_validate() {
        return Err(format!(
            "Euler order is not valid: {euler_order:?}, signs: {signs:?}"
        ));
    }

    // Check if the axes intersect at a common point
    let intersection_point = get_intersection_point(p1, axis1, p2, axis2);

    // Check if the third axis intersects at the same point
    let v = p3 - intersection_point; // Vector from intersection point to p3
    if !approx::relative_eq!(v.cross(&axis2), &Vector3::zeros(), epsilon = epsilon) {
        return Err("Axes do not intersect at a common point".to_string());
    }

    

    todo!()
}

// from link to space line, return (origin, direction)
fn rot_line(l: &Link) -> (Vector3<f64>, Vector3<f64>) {
    (
        l.global_zero_pose.fixed_view::<3, 1>(0, 3).into(), // Origin of the joint
        l.space_spatial_twist.fixed_rows::<3>(0).into(),    // Direction of the joint axis
    )
}

fn axis_to_euler_axis(axis: &Vector3<f64>, epsilon: f64) -> Result<(EulerAxis, f64), String> {
    // Convert a rotation axis to an EulerAxis
    let axis_abs = axis.abs(); // Use absolute values to handle negative axes
    let sign = if axis == &axis_abs { 1.0 } else { -1.0 };
    if approx::relative_eq!(axis_abs, &Vector3::x(), epsilon = epsilon) {
        Ok((EulerAxis::X, sign))
    } else if approx::relative_eq!(axis_abs, &Vector3::y(), epsilon = epsilon) {
        Ok((EulerAxis::Y, sign))
    } else if approx::relative_eq!(axis_abs, &Vector3::z(), epsilon = epsilon) {
        Ok((EulerAxis::Z, sign))
    } else {
        Err(format!(
            "Axis is not aligned with any coordinate axis {axis}"
        ))
    }
}

fn get_intersection_point(
    p1: Vector3<f64>,
    axis1: Vector3<f64>,
    p2: Vector3<f64>,
    axis2: Vector3<f64>,
) -> Vector3<f64> {
    // Check if the axes intersect at a common point
    let v = p2 - p1; // Vector from p1 to p2
    let n = axis1.cross(&axis2); // Normal vector of the plane formed by the two axes
    let t = (v.cross(&axis2)).dot(&n) / n.norm_squared(); // Distance from p1 to the line formed by axis2

    p1 + t * axis1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Model;

    #[test]
    fn test_try_srs() {}

    #[test]
    fn test_is_7dof() {
        let model = Model::from_urdf("./urdf/RM75-6F/urdf/RM75-6F.urdf").unwrap();
        assert!(is_7dof(&model));
    }

    #[test]
    fn test_is_all_revolute_joints() {
        let mut model = Model::from_urdf("./urdf/RM75-6F/urdf/RM75-6F.urdf").unwrap();
        assert!(is_all_revolute_joints(&model));

        model.links[1].joint.urdf_joint.as_mut().unwrap().joint_type =
            urdf_rs::JointType::Continuous;
        assert!(is_all_revolute_joints(&model));

        model.links[1].joint.urdf_joint.as_mut().unwrap().joint_type =
            urdf_rs::JointType::Prismatic;
        assert!(!is_all_revolute_joints(&model));
    }

    #[test]
    fn test_extract_sphere_joint() {
        let model = Model::from_urdf("./urdf/RM75-6F/urdf/RM75-6F.urdf").unwrap();
        let l1 = &model.links[1];
        let l2 = &model.links[2];
        let l3 = &model.links[3];

        // Test extracting spherical joint
        let euler_order = extract_sphere_joint(l1, l2, l3, 1e-5);
        assert!(euler_order.is_ok());
    }

    #[test]
    fn test_rot_line() {
        let model = Model::from_urdf("./urdf/RM75-6F/urdf/RM75-6F.urdf").unwrap();
        let link = &model.links[1];
        // println!("Link: {:#?}", link);
        let (origin, direction) = rot_line(link);

        // println!("Origin: {}", origin);
        // println!("Direction: {}", direction);

        // Check if the origin matches the joint origin
        let joint_origin = link.joint.urdf_joint.as_ref().unwrap().origin.xyz;
        assert_eq!(
            origin,
            Vector3::new(joint_origin[0], joint_origin[1], joint_origin[2])
        );

        // Check if the direction matches the space spatial twist
        assert_eq!(direction, link.space_spatial_twist.fixed_rows::<3>(0));
    }

    #[test]
    fn test_axis_to_euler_axis() {
        let model = Model::from_urdf("./urdf/RM75-6F/urdf/RM75-6F.urdf").unwrap();

        let axis_1 = model.links[1].space_spatial_twist.fixed_rows::<3>(0).into();
        let axis_2 = model.links[2].space_spatial_twist.fixed_rows::<3>(0).into();
        let axis_3 = model.links[3].space_spatial_twist.fixed_rows::<3>(0).into();

        // println!("Axis 1: {}", model.links[1].space_spatial_twist);
        // println!("Axis 2: {}", model.links[2].space_spatial_twist);
        // println!("Axis 3: {}", model.links[3].space_spatial_twist);

        assert_eq!(axis_to_euler_axis(&axis_1, 1e-5), Ok((EulerAxis::Z, -1.)));
        assert_eq!(axis_to_euler_axis(&axis_2, 1e-5), Ok((EulerAxis::Y, -1.)));
        assert_eq!(axis_to_euler_axis(&axis_3, 1e-5), Ok((EulerAxis::Z, -1.)));
    }

    #[test]
    fn test_get_intersection_point() {
        let p1 = Vector3::new(0.0, 0.0, 0.0);
        let axis1 = Vector3::new(0.0, 0.0, 1.0);
        let p2 = Vector3::new(0.0, 1.0, 1.0);
        let axis2 = Vector3::new(0.0, 1.0, 0.0);

        // Calculate the intersection point
        let intersection_point = get_intersection_point(p1, axis1, p2, axis2);

        // Check if the intersection point is correct
        assert_eq!(intersection_point, Vector3::new(0.0, 0.0, 1.0));
    }
}

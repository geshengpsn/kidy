use std::f64::consts::PI;

use nalgebra::{vector, Matrix4, Vector3, VectorView3, U1, U4};
use spatial_alg::{se3, SE3};

use crate::analytical::AikError;

pub struct ZM75AIKSolevr {
    shoulder_length: f64,
    robot_shoulder_point: Vector3<f64>,

    upper_arm_length: f64,
    shoulder_to_j4_length: f64,

    fore_arm_length: f64,
    j4_to_j6_length: f64,

    wrist_length: f64,
    upper_arm_angle: f64,
    fore_arm_angle: f64,

    tool_t_joint7: Matrix4<f64>,
    zero_joint7_t_joint6: Matrix4<f64>,
}

impl Default for ZM75AIKSolevr {
    fn default() -> Self {
        let mut solver = ZM75AIKSolevr {
            shoulder_length: 0.1785,
            robot_shoulder_point: Vector3::zeros(),

            upper_arm_length: 0.28998,
            shoulder_to_j4_length: 0.0,

            fore_arm_length: 0.30879,
            j4_to_j6_length: 0.0,

            wrist_length: 0.11462,
            upper_arm_angle: 0.0,
            fore_arm_angle: 0.0,

            tool_t_joint7: Matrix4::identity(),
            zero_joint7_t_joint6: Matrix4::identity(),
        };
        solver.robot_shoulder_point = vector![0., 0., solver.shoulder_length];
        solver.shoulder_to_j4_length =
            f64::sqrt(solver.upper_arm_length * solver.upper_arm_length + 0.01 * 0.01);
        solver.j4_to_j6_length =
            f64::sqrt(solver.fore_arm_length * solver.fore_arm_length + 0.01 * 0.01);
        solver.upper_arm_angle = f64::atan2(solver.upper_arm_length, 0.01);
        solver.fore_arm_angle = f64::atan2(solver.fore_arm_length, 0.01);

        solver.tool_t_joint7[(2, 3)] = -solver.wrist_length;
        solver.zero_joint7_t_joint6[(2, 3)] = -0.01;
        solver
    }
}

struct ZMAIKResult {}

impl ZM75AIKSolevr {
    fn solve_j4(
        &self,
        robot_joint6_point: VectorView3<f64, U1, U4>,
    ) -> Result<(Option<f64>, Option<f64>), AikError> {
        let shoulder_to_j6_length = (robot_joint6_point - self.robot_shoulder_point).norm();
        if shoulder_to_j6_length > self.shoulder_to_j4_length + self.j4_to_j6_length {
            return Err(AikError::Unreachable);
        }

        // [0, pi)
        let middle_angle = f64::acos(
            (f64::powi(self.shoulder_to_j4_length, 2) + f64::powi(self.j4_to_j6_length, 2)
                - f64::powi(shoulder_to_j6_length, 2))
                / (2. * self.shoulder_to_j4_length * self.j4_to_j6_length),
        );

        // [0, pi)
        let a = 2. * PI - middle_angle - self.fore_arm_angle - self.upper_arm_angle;
        let is_a_in_range = 0.0 < a || a < 130.0;

        // [-pi, +)
        let b = middle_angle - self.fore_arm_angle - self.upper_arm_angle;
        let is_b_in_range = -7.0 < b || b < 0.0;

        Ok((
            if is_a_in_range { Some(a) } else { None },
            if is_b_in_range { Some(b) } else { None },
        ))
    }

    fn aik(&self, target: Matrix4<f64>, joint7_theta: f64) -> Result<ZMAIKResult, AikError> {
        let rot_back = vector![0., joint7_theta, 0., 0., 0., 0.].exp();
        let robot_t_joint6 = target
            * self.tool_t_joint7 // 转到joint7坐标系
            * rot_back                   // 绕y轴旋转j7角度，转为zero_joint7坐标系
            * self.zero_joint7_t_joint6; // 转到joint6坐标系
        let robot_joint6_point = robot_t_joint6.fixed_view::<3, 1>(0, 3);

        // a: [0, pi) / b:[-pi, +)
        let (j4_a, j4_b) = self.solve_j4(robot_joint6_point)?;
        for j4 in [j4_a, j4_b].into_iter().flatten() {
            let joint6_t_robot = robot_t_joint6.inv();
            let shoulder_point_homo = self.robot_shoulder_point.push(1.);
            let joint6_shoulder_point =
                (joint6_t_robot * shoulder_point_homo).fixed_view::<3, 1>(0, 0);

            for i in [1, 2] {}
        }

        Ok(ZMAIKResult {})
    }
}

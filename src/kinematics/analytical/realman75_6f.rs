use std::f64::consts::PI;

use argmin::{
    core::{CostFunction, Executor, Gradient, State},
    solver::{
        brent::BrentOpt, goldensectionsearch::GoldenSectionSearch, particleswarm::ParticleSwarm,
    },
};
use nalgebra::{vector, Matrix3, Matrix4, SVector, Vector3};
use spatial_alg::prelude::*;

use crate::analytical::AikError;

#[derive(Debug)]
pub struct Realman75AikSolver {
    upper_arm: f64,
    fore_arm: f64,
    wrist_point_under_tool_frame: Vector3<f64>,
    shoulder_point: Vector3<f64>,
    min_d: f64,
    max_d: f64,

    joints_limit: [(f64, f64); 7],
    joints_speed_limit: SVector<f64, 7>,
}

impl Default for Realman75AikSolver {
    fn default() -> Self {
        let mut solver = Self {
            upper_arm: 0.256,
            fore_arm: 0.210,
            wrist_point_under_tool_frame: Vector3::new(0.0, 0.0, -0.1613),
            shoulder_point: Vector3::new(0.0, 0.0, 0.2405),
            min_d: 0.,
            max_d: 0.,
            joints_limit: [
                (-178. * PI / 180., 178. * PI / 180.),
                (-130. * PI / 180., 130. * PI / 180.),
                (-178. * PI / 180., 178. * PI / 180.),
                (-135. * PI / 180., 135. * PI / 180.),
                (-178. * PI / 180., 178. * PI / 180.),
                (-128. * PI / 180., 128. * PI / 180.),
                (-360. * PI / 180., 360. * PI / 180.),
            ],
            joints_speed_limit: vector![
                180. * PI / 180., // 180 deg/s
                180. * PI / 180., // 180 deg/s
                225. * PI / 180., // 225 deg/s
                225. * PI / 180., // 225 deg/s
                225. * PI / 180., // 225 deg/s
                225. * PI / 180., // 225 deg/s
                225. * PI / 180., // 225 deg/s
            ],
        };
        let rad = PI - solver.joints_limit[3].1;
        solver.min_d = (solver.upper_arm.powi(2) + solver.fore_arm.powi(2)
            - 2. * solver.upper_arm * solver.fore_arm * rad.cos())
        .sqrt();
        solver.max_d = solver.upper_arm + solver.fore_arm;
        solver
    }
}

impl Realman75AikSolver {
    fn wrist_point(&self, target_pose: &Matrix4<f64>) -> Vector3<f64> {
        let wrist_point_homo = target_pose * self.wrist_point_under_tool_frame.push(1.0);
        Vector3::new(
            wrist_point_homo[0],
            wrist_point_homo[1],
            wrist_point_homo[2],
        )
    }

    fn is_reachable(&self, target_pose: &Matrix4<f64>) -> bool {
        let wrist_point = self.wrist_point(target_pose);
        let sw = wrist_point - self.shoulder_point;
        let d = sw.norm();
        self.min_d <= d && d <= self.max_d
    }

    fn get_middle_result(
        &self,
        target_pose: &Matrix4<f64>,
        config: ElbowConfig,
    ) -> Result<MiddleResult, AikError> {
        // shoulder to wrist point
        let wrist_point = self.wrist_point(target_pose);
        let sw = wrist_point - self.shoulder_point;

        // shoulder to wrist length
        let d = sw.norm();
        if d < self.min_d || d > self.max_d {
            return Err(AikError::Unreachable);
        }

        let d_prime = (self.upper_arm.powi(2) - self.fore_arm.powi(2) + d * d) / (2.0 * d);
        let r = (self.upper_arm.powi(2) - d_prime.powi(2)).sqrt();
        let ep_z = sw / d;
        let ep_x = {
            let ep_x = Vector3::new(0.0, 0.0, 1.0).cross(&ep_z);
            if ep_x.norm().abs() < f64::EPSILON {
                Vector3::new(1.0, 0.0, 0.0) // if ep_z is aligned with z-axis, use x-axis as ep_x
            } else {
                ep_x / ep_x.norm()
            }
        };

        let ep_point = self.shoulder_point + ep_z * d_prime;
        let zero_elbow_point = ep_x * r + ep_point;
        let se = zero_elbow_point - self.shoulder_point;
        let ew = wrist_point - zero_elbow_point;

        let zero_elbow_y = {
            let zero_elbow_y = ep_x.cross(&ep_z);
            zero_elbow_y / zero_elbow_y.norm()
        };

        let se_z = se / se.norm();
        let se_x = {
            let se_x = zero_elbow_y.cross(&se_z);
            se_x / se_x.norm()
        };

        let ew_z = ew / ew.norm();
        let ew_x = {
            let ew_x = zero_elbow_y.cross(&ew_z);
            ew_x / ew_x.norm()
        };

        let inner_elbow_angle = se_z.dot(&ew_z).acos();
        let robot_r_tool = target_pose.fixed_view::<3, 3>(0, 0).into();
        match config {
            ElbowConfig::Inner => Ok(MiddleResult {
                se: Matrix3::from_columns(&[se_x, zero_elbow_y, se_z]),
                ew: Matrix3::from_columns(&[ew_x, zero_elbow_y, ew_z]),
                ep_z,
                inner_elbow_angle,
                robot_r_tool,
            }),
            ElbowConfig::Outer => Ok(MiddleResult {
                se: Matrix3::from_columns(&[-se_x, -zero_elbow_y, se_z]),
                ew: Matrix3::from_columns(&[-ew_x, -zero_elbow_y, ew_z]),
                ep_z,
                inner_elbow_angle: -inner_elbow_angle,
                robot_r_tool,
            }),
        }
    }
}

struct MiddleResult {
    se: Matrix3<f64>,
    ew: Matrix3<f64>,
    ep_z: Vector3<f64>,
    inner_elbow_angle: f64,
    robot_r_tool: Matrix3<f64>,
}

fn zyz(m: &Matrix3<f64>, config: SphereConfig) -> SVector<f64, 3> {
    let r23 = m[(1, 2)];
    let r13 = m[(0, 2)];
    let r31 = m[(2, 0)];
    let r32 = m[(2, 1)];
    let r33 = m[(2, 2)];
    let is_singular = {
        r23.abs() < f64::EPSILON
            && r13.abs() < f64::EPSILON
            && r31.abs() < f64::EPSILON
            && r32.abs() < f64::EPSILON
    };
    let theta1 = r23.atan2(r13);
    let theta2 = (r31 * r31 + r32 * r32).sqrt().atan2(r33);
    let theta3 = r32.atan2(-r31);
    if is_singular {
        SVector::<f64, 3>::new(0., 0., m[(0, 0)].atan2(m[(1, 0)]))
    } else {
        match config {
            SphereConfig::A => SVector::<f64, 3>::new(theta1, theta2, theta3),
            SphereConfig::B => SVector::<f64, 3>::new(
                if theta1 < PI {
                    theta1 + PI
                } else {
                    theta1 - PI
                },
                -theta2,
                if theta3 < PI {
                    theta3 + PI
                } else {
                    theta3 - PI
                },
            ),
        }
    }
}

impl Realman75AikSolver {
    fn aik_from_middle_result(
        &self,
        mr: &MiddleResult,
        theta: f64,
        config: SolutionConfig,
    ) -> SVector<f64, 7> {
        let r = (mr.ep_z * theta).exp();
        let shoulder_r = r * mr.se;
        let wrist_r = (mr.robot_r_tool.transpose() * r * mr.ew).transpose();
        let shoulder = zyz(&shoulder_r, config.shoulder);
        let wrist = zyz(&wrist_r, config.wrist);
        vector![
            shoulder[0],
            shoulder[1],
            shoulder[2],
            mr.inner_elbow_angle,
            wrist[0],
            wrist[1],
            wrist[2],
        ]
    }

    pub fn aik(
        &self,
        target_pose: &Matrix4<f64>,
        theta: f64,
        config: SolutionConfig,
    ) -> Result<SVector<f64, 7>, AikError> {
        let mr = self.get_middle_result(target_pose, config.elbow)?;
        Ok(self.aik_from_middle_result(&mr, theta, config))
    }

    pub fn is_valid_joints(&self, joints: &SVector<f64, 7>) -> bool {
        for (angle, limit) in joints.iter().zip(self.joints_limit) {
            if angle < &limit.0 || angle > &limit.1 {
                return false;
            }
        }
        true
    }

    pub fn is_valid_speed(&self, diff: &SVector<f64, 7>, duration: f64) -> bool {
        let speed = diff.abs() / duration;
        speed < self.joints_speed_limit
    }

    // fn aik_grad_from_middle_result(
    //     &self,
    //     mr: &MiddleResult,
    //     theta: f64,
    //     config: SolutionConfig,
    // ) -> (SVector<f64, 7>, SVector<f64, 7>) {
    //     let r = (mr.ep_z * theta).exp();
    //     let shoulder_r = r * mr.se;
    //     let wrist_r = (mr.robot_r_tool.transpose() * r * mr.ew).transpose();
    //     let shoulder = zyz(&shoulder_r, config.shoulder);
    //     let wrist = zyz(&wrist_r, config.wrist);
    //     let joints = vector![
    //         shoulder[0],
    //         shoulder[1],
    //         shoulder[2],
    //         mr.inner_elbow_angle,
    //         wrist[0],
    //         wrist[1],
    //         wrist[2],
    //     ];

    //     let dshoulder_dtheta = hat(&mr.ep_z) * r * mr.se;
    //     let dwrist_dtheta = -mr.ew.transpose() * r.transpose() * hat(&mr.ep_z) * mr.robot_r_tool;
    //     let (ds1_dshoulder, ds2_dshoulder, ds3_dshoulder) = zyz_jacobian(&shoulder_r);
    //     let (dw1_dwrist, dw2_dwrist, dw3_dwrist) = zyz_jacobian(&wrist_r);
    //     let ele_mul_sum = |m1: Matrix3<f64>, m2: Matrix3<f64>| {
    //         m1.iter().zip(m2.iter()).map(|(a, b)| a * b).sum::<f64>()
    //     };
    //     (
    //         joints,
    //         vector![
    //             ele_mul_sum(ds1_dshoulder, dshoulder_dtheta),
    //             ele_mul_sum(ds2_dshoulder, dshoulder_dtheta),
    //             ele_mul_sum(ds3_dshoulder, dshoulder_dtheta),
    //             0.,
    //             ele_mul_sum(dw1_dwrist, dwrist_dtheta),
    //             ele_mul_sum(dw2_dwrist, dwrist_dtheta),
    //             ele_mul_sum(dw3_dwrist, dwrist_dtheta),
    //         ],
    //     )
    // }

    // pub fn aik_grad(
    //     &self,
    //     target_pose: &Matrix4<f64>,
    //     theta: f64,
    //     config: SolutionConfig,
    // ) -> Result<(SVector<f64, 7>, SVector<f64, 7>), AikError> {
    //     let mr = self.get_middle_result(target_pose, config.elbow)?;
    //     Ok(self.aik_grad_from_middle_result(&mr, theta, config))
    // }
}

// fn zyz_jacobian(m: &Matrix3<f64>) -> (Matrix3<f64>, Matrix3<f64>, Matrix3<f64>) {
//     let r23 = m[(1, 2)];
//     let r13 = m[(0, 2)];
//     let r31 = m[(2, 0)];
//     let r32 = m[(2, 1)];
//     let r33 = m[(2, 2)];
//     let s_square = r23.powi(2) + r13.powi(2);
//     let mut alpha = Matrix3::zeros();
//     let mut beta = Matrix3::zeros();
//     let mut gamma = Matrix3::zeros();
//     if s_square < f64::EPSILON {
//         // m[(0, 0)].atan2(m[(1, 0)])
//         gamma[(0, 0)] = m[(1, 0)];
//         gamma[(1, 0)] = -m[(0, 0)];
//         return (alpha, beta, gamma);
//     }
//     let s = s_square.sqrt();
//     alpha[(0, 2)] = -r23 / s_square;
//     alpha[(1, 2)] = r13 / s_square;
//     beta[(2, 0)] = r33 * r31 / s;
//     beta[(2, 1)] = r33 * r32 / s;
//     beta[(2, 2)] = -s;
//     gamma[(2, 0)] = r32 / s_square;
//     gamma[(2, 1)] = -r31 / s_square;
//     (alpha, beta, gamma)
// }

// #[test]
// fn test_aik_grad() {
//     use nalgebra::matrix;
//     use std::f64::consts::FRAC_PI_2;
//     let solver = Realman75AikSolver::default();
//     let target = matrix![
//         0., 0., 1., 0.379;
//         1., 0., 0., 0.;
//         0., 1., 0., 0.6;
//         0., 0., 0., 1.
//     ];
//     let (res1, grad1) = solver
//         .aik_grad(&target, FRAC_PI_2, ElbowConfig::Inner)
//         .unwrap();
//     let (res2, grad2) = solver
//         .aik_grad(&target, FRAC_PI_2 + 1e-3, ElbowConfig::Inner)
//         .unwrap();
//     let num_grad = (res2 - res1) / 1e-3;
//     println!("num_grad: {:.3}", num_grad);
//     println!("grad1: {:.3}", grad1);
//     println!("grad2: {:.3}", grad2);
// }

impl Realman75AikSolver {
    pub fn search_pso<F: CostFunction<Param = SVector<f64, 7>, Output = f64>>(
        &self,
        target: &Matrix4<f64>,
        config: SolutionConfig,
        f: F,
    ) -> Result<f64, AikError> {
        let pso = ParticleSwarm::new((vec![0.], vec![2. * PI]), 10);
        if !self.is_reachable(target) {
            return Err(AikError::Unreachable);
        }
        let problem = ArmAngleProblem::init(self, *target, config, f).unwrap();
        let res = Executor::new(problem, pso)
            .configure(|state| state.max_iters(100))
            .run()
            .unwrap();
        Ok(res.state.get_best_param().unwrap().position[0])
    }

    pub fn search_brent<F: CostFunction<Param = SVector<f64, 7>, Output = f64>>(
        &self,
        target: &Matrix4<f64>,
        config: SolutionConfig,
        f: F,
    ) -> Result<f64, AikError> {
        let brent = BrentOpt::new(0., 2. * PI);
        let problem = ArmAngleProblem::init(self, *target, config, f).unwrap();
        let res = Executor::new(problem, brent)
            .configure(|state| state.param(0.).max_iters(10))
            .run()
            .unwrap();
        Ok(*res.state.get_best_param().unwrap())
    }

    pub fn golden_section_search<F: CostFunction<Param = SVector<f64, 7>, Output = f64>>(
        &self,
        target: &Matrix4<f64>,
        config: SolutionConfig,
        f: F,
        init: f64,
    ) -> Result<f64, AikError> {
        let golden_section = GoldenSectionSearch::new(0., 2. * PI)
            .unwrap()
            .with_tolerance(1e-4)
            .unwrap();
        let problem = ArmAngleProblem::init(self, *target, config, f).unwrap();
        let res = Executor::new(problem, golden_section)
            .configure(|state| state.param(init).max_iters(100))
            .run()
            .unwrap();
        Ok(*res.state.get_best_param().unwrap())
    }

    // steepest descent method
    // pub fn search_sd<F>(
    //     &self,
    //     target: &Matrix4<f64>,
    //     ik_solver: &Realman75AikSolver,
    //     elbow: ElbowConfig,
    //     f: F,
    // ) -> Result<f64, AikError>
    // where
    //     F: CostFunction<Param = SVector<f64, 7>, Output = f64>
    //         + Gradient<Param = SVector<f64, 7>, Gradient = SVector<f64, 7>>,
    // {
    //     let ls = linesearch::BacktrackingLineSearch::new(
    //         linesearch::condition::ArmijoCondition::new(0.1f64).unwrap(),
    //     );
    //     let sd = SteepestDescent::new(ls);
    //     let problem: ArmAngleProblem<'_, F, f64> =
    //         ArmAngleProblem::init(ik_solver, *target, elbow, f).unwrap();
    //     let res = Executor::new(problem, sd)
    //         .configure(|state| state.param(0.).max_iters(10))
    //         .run()
    //         .unwrap();
    //     Ok(*res.state.get_best_param().unwrap())
    // }
}

#[test]
fn test_aik_and_aik_grad_time() {
    use nalgebra::matrix;
    use std::f64::consts::FRAC_PI_2;
    use std::time::Instant;
    let solver = Realman75AikSolver::default();
    let target = matrix![
        0., 0., 1., 0.379;
        1., 0., 0., 0.;
        0., 1., 0., 0.6;
        0., 0., 0., 1.
    ];
    let start = Instant::now();
    let res = solver
        .aik(
            &target,
            FRAC_PI_2,
            SolutionConfig {
                shoulder: SphereConfig::A,
                elbow: ElbowConfig::Inner,
                wrist: SphereConfig::A,
            },
        )
        .unwrap();

    println!("aik {:?}", start.elapsed());
    println!("{res}");

    // let mr = solver
    //     .get_middle_result(&target, ElbowConfig::Inner)
    //     .unwrap();
    // let start = Instant::now();
    // let res = solver.aik_from_middle_result(&mr, FRAC_PI_2);

    // println!("aik from mr {:?}", start.elapsed());
    // println!("{res}");

    // target[(1, 3)] = 0.2; // change y position to test aik_grad
    // let start = Instant::now();
    // let (res, _grad) = solver
    //     .aik_grad(&target, FRAC_PI_2, ElbowConfig::Inner)
    //     .unwrap();
    // println!("aik_grad {:?}", start.elapsed());
    // println!("{res}");

    // let start = Instant::now();
    // let (res, _grad) = solver.aik_grad_from_middle_result(&mr, FRAC_PI_2);
    // println!("aik_grad from mr {:?}", start.elapsed());
    // println!("{res}");
}

#[test]
fn test_search_pso() {
    use nalgebra::matrix;
    use std::time::Instant;
    let solver = Realman75AikSolver::default();
    let target = matrix![
        0., 0., 1., 0.379;
        1., 0., 0., 0.;
        0., 1., 0., 0.6;
        0., 0., 0., 1.
    ];
    let start = Instant::now();
    let res = solver
        .search_pso(
            &target,
            SolutionConfig {
                shoulder: SphereConfig::A,
                elbow: ElbowConfig::Inner,
                wrist: SphereConfig::A,
            },
            JointsNorm,
        )
        .unwrap();
    println!("res: {}", res);
    println!("time: {:?}", start.elapsed());
}

#[test]
fn test_search_brent() {
    use nalgebra::matrix;
    use std::time::Instant;
    let solver = Realman75AikSolver::default();
    let target = matrix![
        0., 0., 1., 0.379;
        1., 0., 0., 0.;
        0., 1., 0., 0.6;
        0., 0., 0., 1.
    ];
    let start = Instant::now();
    let res = solver
        .search_brent(
            &target,
            SolutionConfig {
                shoulder: SphereConfig::A,
                elbow: ElbowConfig::Inner,
                wrist: SphereConfig::A,
            },
            JointsNorm,
        )
        .unwrap();
    println!("res: {}", res);
    println!("time: {:?}", start.elapsed());
}

#[test]
fn test_golden_section_search() {
    use nalgebra::matrix;
    use std::time::Instant;
    let solver = Realman75AikSolver::default();
    let target = matrix![
        0., 0., 1., 0.379;
        1., 0., 0., 0.;
        0., 1., 0., 0.6;
        0., 0., 0., 1.
    ];
    let start = Instant::now();
    let res = solver
        .golden_section_search(
            &target,
            SolutionConfig {
                shoulder: SphereConfig::A,
                elbow: ElbowConfig::Inner,
                wrist: SphereConfig::A,
            },
            JointsNorm,
            0.,
        )
        .unwrap();
    println!("res: {}", res);
    println!("time: {:?}", start.elapsed());
}

// #[test]
// fn test_steepest_descent() {
//     use nalgebra::matrix;
//     use std::time::Instant;
//     let solver = Realman75AikSolver::default();
//     let target = matrix![
//         0., 0., 1., 0.379;
//         1., 0., 0., 0.;
//         0., 1., 0., 0.6;
//         0., 0., 0., 1.
//     ];
//     let start = Instant::now();
//     let res = solver
//         .search_sd(&target, &solver, ElbowConfig::Inner, JointsNorm)
//         .unwrap();
//     println!("res: {}", res);
//     println!("time: {:?}", start.elapsed());
// }

pub struct JointsNorm;
impl CostFunction for JointsNorm {
    type Param = SVector<f64, 7>;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        Ok(param.norm())
    }
}

impl Gradient for JointsNorm {
    type Param = SVector<f64, 7>;
    type Gradient = SVector<f64, 7>;

    fn gradient(&self, param: &Self::Param) -> Result<Self::Gradient, argmin::core::Error> {
        Ok(param * 2.)
    }
}

struct ArmAngleProblem<'a, F, P> {
    solver: &'a Realman75AikSolver,
    mr: MiddleResult,
    config: SolutionConfig,
    function: F,
    _param: std::marker::PhantomData<P>,
}

impl<'a, F, P> ArmAngleProblem<'a, F, P> {
    fn init(
        solver: &'a Realman75AikSolver,
        target: Matrix4<f64>,
        config: SolutionConfig,
        f: F,
    ) -> Result<Self, AikError> {
        Ok(ArmAngleProblem {
            solver,
            mr: solver.get_middle_result(&target, config.elbow)?,
            config,
            function: f,
            _param: std::marker::PhantomData,
        })
    }
}

impl<F, P> CostFunction for ArmAngleProblem<'_, F, P>
where
    F: CostFunction<Param = SVector<f64, 7>, Output = f64>,
    for<'a> &'a P: Tof64,
{
    type Param = P;

    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        let joints = self
            .solver
            .aik_from_middle_result(&self.mr, param.f64(), self.config);
        self.function.cost(&joints)
    }
}

// impl<F, P> Gradient for ArmAngleProblem<'_, F, P>
// where
//     F: Gradient<Param = SVector<f64, 7>, Gradient = SVector<f64, 7>>,
//     for<'a> &'a P: Tof64,
// {
//     type Param = P;

//     type Gradient = f64;

//     fn gradient(&self, param: &Self::Param) -> Result<Self::Gradient, argmin_math::Error> {
//         let (joints, grad) = self
//             .solver
//             .aik_grad_from_middle_result(&self.mr, param.f64(), self.config);
//         Ok(self.function.gradient(&joints)?.dot(&grad))
//     }
// }

trait Tof64 {
    fn f64(self) -> f64;
}

impl Tof64 for f64 {
    fn f64(self) -> f64 {
        self
    }
}

impl Tof64 for &f64 {
    fn f64(self) -> f64 {
        *self
    }
}

impl Tof64 for &Vec<f64> {
    fn f64(self) -> f64 {
        self[0]
    }
}

// struct ArmAngleProblem<'a, F>
// where
//     F: Fn(&Realman75AikSolver, &SVector<f64, 7>) -> f64,
// {
//     solver: &'a Realman75AikSolver,
//     mr: MiddleResult,
//     target: Matrix4<f64>,
//     cost: F,
// }

// impl<'a, F> ArmAngleProblem<'a, F>
// where
//     F: Fn(&Realman75AikSolver, &SVector<f64, 7>) -> f64,
// {
//     fn init(
//         solver: &'a Realman75AikSolver,
//         target: Matrix4<f64>,
//         elbow: ElbowConfig,
//         f: F,
//     ) -> Result<Self, AikError> {
//         Ok(ArmAngleProblem {
//             solver,
//             mr: solver.get_middle_result(&target, elbow)?,
//             target,
//             cost: f,
//         })
//     }
// }

// impl<F> CostFunction for ArmAngleProblem<'_, F>
// where
//     F: Fn(&Realman75AikSolver, &SVector<f64, 7>) -> f64,
// {
//     type Param = f64;

//     type Output = f64;

//     fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
//         let joints = self
//             .solver
//             .aik_from_middle_result(&self.target, &self.mr, *param);
//         Ok((self.cost)(self.solver, &joints))
//     }
// }

#[derive(Debug, Clone, Copy)]
pub struct SolutionConfig {
    pub shoulder: SphereConfig,
    pub elbow: ElbowConfig,
    pub wrist: SphereConfig,
}

// impl ChooseConfig {
//     pub fn aa() -> Self {
//         Self {
//             shoulder: SphereConfig::A,
//             wrist: SphereConfig::A,
//         }
//     }

//     pub fn ab() -> Self {
//         Self {
//             shoulder: SphereConfig::A,
//             wrist: SphereConfig::B,
//         }
//     }

//     pub fn ba() -> Self {
//         Self {
//             shoulder: SphereConfig::B,
//             wrist: SphereConfig::A,
//         }
//     }

//     pub fn bb() -> Self {
//         Self {
//             shoulder: SphereConfig::B,
//             wrist: SphereConfig::B,
//         }
//     }
// }

#[derive(Debug, Clone, Copy)]
pub enum ElbowConfig {
    Inner,
    Outer,
}

#[derive(Debug, Clone, Copy)]
pub enum SphereConfig {
    A,
    B,
}

use nalgebra::{Matrix4, Matrix6xX, MatrixXx6, Vector3, Vector6};
use spatial_alg::prelude::*;
use std::f64::consts::FRAC_PI_2;

// use crate::model::Model;
use crate::model::Model;

#[cfg(feature = "any_ik")]
pub mod analytical;

#[derive(Debug)]
pub enum IkError {
    PseudoInverse,
    MaxIter,
}

#[derive(Debug)]
pub struct IkSolveParam {
    pub r_error: f64,
    pub p_error: f64,
    /// jacobian inverse eps
    pub pinv_eps: f64,
    pub max_time: usize,
}

impl Model {
    pub fn fk(&self, q: &[f64]) -> Vec<Matrix4<f64>> {
        assert_eq!(q.len(), self.links.len());
        let mut link_poses = vec![Matrix4::identity(); q.len()];
        let log_twists = &self
            .bfs
            .iter()
            .zip(q.iter())
            .map(|(i, q)| (self.links[*i].space_spatial_twist * *q).exp())
            .collect::<Vec<_>>();
        let mut poes = log_twists.clone();
        for link_index in &self.bfs {
            let parent_index = {
                let parent_index = self
                    .link_graph
                    .neighbors_directed(*link_index, petgraph::Direction::Incoming)
                    .next();
                match parent_index {
                    Some(i) => i,
                    None => continue,
                }
            };
            let parent_poe = poes[parent_index];
            poes[*link_index] = parent_poe * log_twists[*link_index];
            link_poses[*link_index] = poes[*link_index] * self.links[*link_index].global_zero_pose;
        }
        link_poses
    }

    pub fn diff_fk(&self, q: &[f64]) -> (Matrix6xX<f64>, Vec<Matrix4<f64>>) {
        assert_eq!(q.len(), self.links.len());
        let mut link_poses = vec![Matrix4::identity(); q.len()];
        let mut jacobian = Matrix6xX::zeros(q.len());
        let log_twists = &self
            .bfs
            .iter()
            .zip(q.iter())
            .map(|(i, q)| (self.links[*i].space_spatial_twist * *q).exp())
            .collect::<Vec<_>>();
        let mut poes = log_twists.clone();
        for link_index in &self.bfs {
            let parent_index = {
                let parent_index = self
                    .link_graph
                    .neighbors_directed(*link_index, petgraph::Direction::Incoming)
                    .next();
                match parent_index {
                    Some(i) => i,
                    None => continue,
                }
            };
            let parent_poe = poes[parent_index];

            let jn = parent_poe.lee_adjoint() * self.links[*link_index].space_spatial_twist;
            jacobian
                .fixed_view_mut::<6, 1>(0, *link_index)
                .copy_from_slice(jn.as_slice());

            poes[*link_index] = parent_poe * log_twists[*link_index];
            link_poses[*link_index] = poes[*link_index] * self.links[*link_index].global_zero_pose;
        }
        (jacobian, link_poses)
    }

    pub fn ik(
        &self,
        link_id: usize,
        target_pose: Matrix4<f64>,
        init_q: &[f64],
        ik_param: &IkSolveParam,
    ) -> Result<(Vec<f64>, usize), IkError> {
        // let target_pose = self.link_poses[link_id].clone();
        let mut q = init_q.to_vec();
        for t in 0..ik_param.max_time {
            let (jacobian, poses) = self.diff_fk(&q);
            let current_pose = &poses[link_id];

            // T_c^-1 * T_t = c_T_t
            let t_ct = (current_pose.inv() * target_pose).log();
            let vs = current_pose.lee_adjoint() * t_ct;
            if Vector3::new(vs[0], vs[1], vs[2]).norm() <= ik_param.r_error
                && Vector3::new(vs[3], vs[4], vs[5]).norm() <= ik_param.p_error
            {
                return Ok((q, t));
            }
            let Ok(pinv_j) = jacobian.clone().pseudo_inverse(ik_param.pinv_eps) else {
                return Err(IkError::PseudoInverse);
            };

            let v = Vector6::from_row_slice(vs.as_slice());
            let update = pinv_j * v;
            for (j, u) in q.iter_mut().zip(update.iter()) {
                *j += *u;
            }
        }
        Err(IkError::MaxIter)
    }

    pub fn dls_ik(
        &self,
        link_id: usize,
        target_pose: Matrix4<f64>,
        init_q: &[f64],
        ik_param: &IkSolveParam,
    ) -> Result<(Vec<f64>, usize), IkError> {
        let mut q = init_q.to_vec();
        for t in 0..ik_param.max_time {
            let (jacobian, poses) = self.diff_fk(&q);
            let current_pose = &poses[link_id];

            // T_c^-1 * T_t = c_T_t
            let t_ct = (current_pose.inv() * target_pose).log();
            let vs = current_pose.lee_adjoint() * (t_ct);
            if Vector3::new(vs[0], vs[1], vs[2]).norm() <= ik_param.r_error
                && Vector3::new(vs[3], vs[4], vs[5]).norm() <= ik_param.p_error
            {
                return Ok((q, t));
            }
            let dls_pinv_j = dls_inverse(jacobian, ik_param.pinv_eps);
            let v = Vector6::from_row_slice(vs.as_slice());
            let mut update = dls_pinv_j * v;

            let clamp_d = FRAC_PI_2;
            if update.norm() > clamp_d {
                update = clamp_d * update.normalize();
            }

            for (j, u) in q.iter_mut().zip(update.iter()) {
                *j += *u;
            }
        }
        Err(IkError::MaxIter)
    }
}

fn dls_inverse(m: Matrix6xX<f64>, eps: f64) -> MatrixXx6<f64> {
    let mut svd = m.svd(true, true);
    let lambda = svd.singular_values.min() / 2.;
    for i in 0..svd.singular_values.len() {
        let val = svd.singular_values[i];
        // let lambda = val;
        if val > eps {
            svd.singular_values[i] = val / (val * val + lambda * lambda);
        } else {
            svd.singular_values[i] = 0.;
        }
    }

    svd.recompose().map(|m| m.adjoint()).unwrap()
}

#[cfg(test)]
mod test {
    use std::f64::consts::FRAC_PI_2;

    use crate::model::Model;

    #[test]
    fn test_fk() {
        let model = Model::from_urdf("./urdf/rm_75_6fb_description/urdf/RM75-6F.urdf").unwrap();
        let q = [0.0, FRAC_PI_2, 0., 0., 0., 0., 0., 0.];
        // first joint is fixed joint
        // state.set_q(&[0.0, FRAC_PI_2, 0., 0., 0., 0., 0., 0.]);
        let poses = model.fk(&q);
        poses.iter().enumerate().for_each(|(i, pose)| {
            println!("link {}: {:.2}", i, pose);
        });
    }

    #[test]
    fn test_diff_fk() {
        let model = Model::from_urdf("./urdf/rm_75_6fb_description/urdf/RM75-6F.urdf").unwrap();
        let q = [0.0, FRAC_PI_2, 0., 0., 0., 0., 0., 0.];
        let (j, _poses) = model.diff_fk(&q);
        println!("{:.2}", j);
    }

    #[test]
    fn test_ik() {}
}

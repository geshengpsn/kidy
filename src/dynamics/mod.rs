use crate::Model;
use nalgebra::Vector6;
use spatial_alg::prelude::*;

impl Model {
    /// RNEA inverse dynamics
    ///
    /// typical gravity: [0., 0., -9.8]
    pub fn id(
        &self,
        q: &[f64],
        dq: &[f64],
        ddq: &[f64],
        external_force: Vector6<f64>,
        gravity: &[f64; 3],
    ) -> Vec<Vector6<f64>> {
        let n_link = self.bfs.len();
        assert_eq!(q.len(), n_link);
        assert_eq!(dq.len(), n_link);
        assert_eq!(ddq.len(), n_link);
        let mut twists = vec![Vector6::zeros(); n_link];
        let mut dtwists = vec![Vector6::zeros(); n_link];

        // // add gravity
        let g = Vector6::new(0.0, 0.0, 0.0, gravity[0], gravity[1], gravity[2]);
        dtwists[0] = g;

        for i in self.bfs.iter().skip(1).cloned() {
            // A_i
            let local_screw = &self.links[i].local_spatial_twist;

            // i^T_i-1: transformation from i to i-1
            let t = (local_screw * -q[i]).exp() * self.links[i].parent_zero_pose.inv();

            // V_i: twist of link i
            let v: Vector6<f64> = t.lee_adjoint() * twists[i - 1] + local_screw * dq[i];

            // dV_i: acceleration of link i
            let dv = t.lee_adjoint() * (dtwists[i - 1])
                + (v.adj() * local_screw) * dq[i]
                + local_screw * ddq[i];

            twists[i] = v;
            dtwists[i] = dv;
        }

        // for t in &twists {
        //     println!("{}", t.vee());
        // }

        // for t in &dtwists {
        //     println!("{}", t.vee());
        // }

        let mut joint_forces = vec![Vector6::zeros(); n_link];
        for i in self.bfs.iter().skip(1).rev().cloned() {
            // println!("{i}");
            let f_child_link = {
                if i == n_link - 1 {
                    // se3::identity()
                    external_force
                } else {
                    let t = (self.links[i + 1].local_spatial_twist * -q[i + 1]).exp()
                        * self.links[i + 1].parent_zero_pose.inv();
                    t.lee_adjoint().transpose() * joint_forces[i + 1]
                }
            };
            // println!("{}", f_child_link.vee());

            // println!("i: {}", self.links[i].local_spatial_inertial);

            joint_forces[i] = f_child_link + self.links[i].local_spatial_inertial * dtwists[i]
                - (twists[i].adj().transpose()
                    * (self.links[i].local_spatial_inertial * twists[i]));
            // println!("f:{}", joint_forces[i].vee());
        }
        joint_forces
    }

    // pub fn mass_matrix(
    //     &self,
    //     q: &[f64],
    //     // external_forces: &[se3<f64>],
    //     gravity: [f64; 3],
    // ) -> DMatrix<f64> {
    //     let n = self.bfs.len();
    //     let a = self.id(q, &vec![0.; n], &vec![0.; n], se3::identity(), gravity);
    // }

    pub fn vel_quadratic_forces(&self, q: &[f64], dq: &[f64]) -> Vec<Vector6<f64>> {
        let n = self.bfs.len();
        self.id(q, dq, &vec![0.; n], Vector6::zeros(), &[0.; 3])
        // self.get_actuator_torque(&joint_forces)
    }

    pub fn gravity_forces(&self, q: &[f64], gravity: &[f64; 3]) -> Vec<Vector6<f64>> {
        let n = self.bfs.len();
        self.id(q, &vec![0.; n], &vec![0.; n], Vector6::zeros(), gravity)
        // self.get_actuator_torque(&joint_forces)
    }

    pub fn get_actuator_torque(&self, joint_forces: &[Vector6<f64>]) -> Vec<f64> {
        assert_eq!(self.links.len(), joint_forces.len());
        let mut torques = Vec::with_capacity(joint_forces.len());
        for (l, f) in self.links.iter().zip(joint_forces.iter()) {
            // println!("{:.3} {:.3} {}", l.local_spatial_screw.vee(), f.vee(), cross_se3(&l.local_spatial_screw, f));
            torques.push(l.local_spatial_twist.dot(f));
        }
        torques
    }
}

// fn ad_mul_se3(m: Matrix6<f64>, w: &Vector6<f64>) -> Vector6<f64> {
//     let vec6 = Vector6::from_column_slice(w.as_slice());
//     let res = m * w;
//     Vector6::new([res[0], res[1], res[2]], [res[3], res[4], res[5]])
// }

// fn cross_se3(v: &se3<f64>, w: &se3<f64>) -> f64 {
//     let vec_v = Vector6::from_column_slice(v.vee().as_slice());
//     let vec_w = Vector6::from_column_slice(w.vee().as_slice());
//     vec_v.dot(&vec_w)
// }

#[cfg(test)]
mod test {
    use nalgebra::Vector6;

    // use super::*;
    use crate::Model;
    // use petgraph::prelude::DiGraphMap;

    #[test]
    fn test_id() {
        let model = Model::from_urdf("./urdf/rm_75_6fb_description/urdf/RM75-6F.urdf").unwrap();
        let q = [0.; 8];
        let dq = [0.; 8];
        let ddq = [0.; 8];
        let joint_forces = model.id(&q, &dq, &ddq, Vector6::zeros(), &[0., 0., -9.8]);
        for joint_f in joint_forces {
            println!("{}", joint_f);
        }
    }

    #[test]
    fn test_gravity_forces() {
        let model = Model::from_urdf("./urdf/rm_75_6fb_description/urdf/RM75-6F.urdf").unwrap();
        let mut q = [0.; 8];
        q[2] = 1.58;
        let gravity = [0., 0., -9.81];
        let forces = model.gravity_forces(&q, &gravity);
        let tau = model.get_actuator_torque(&forces);
        for f in tau {
            println!("{}", f);
        }
    }
}

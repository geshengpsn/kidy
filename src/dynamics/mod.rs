use crate::KidyChain;
use liealg::{prelude::*, se3};
use nalgebra::{Matrix6, Vector3, Vector6};

impl<const N: usize> KidyChain<N> {
    // tau = M(q) * q_dot_dot + C(q, q_dot) * q_dot + G(q)
    pub fn id(
        &self,
        theta: [f64; N],
        dtheta: [f64; N],
        ddtheta: [f64; N],
        gravity: [f64; 3],
        f_tip: se3<f64>,
    ) -> [f64; N] {
        let mut twists = vec![se3::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]); N + 1];
        let mut dtwists = vec![se3::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]); N + 1];

        // add gravity
        let g = se3::new([0.0, 0.0, 0.0], gravity);
        dtwists[0] = g;

        // i^T_i-1
        for i in 1..N + 1 {
            // i^T_i-1: transformation from i to i-1
            let t = (self.local_screw[i - 1].clone() * -theta[i - 1]).exp()
                * self.local_zero_pose[i].clone().inv();

            // V_i: twist of link i
            let v =
                t.adjoint().act(&twists[i - 1]) + self.local_screw[i - 1].clone() * dtheta[i - 1];
            twists[i] = v.clone();

            // dV_i: acceleration of link i
            let dv = t.adjoint().act(&dtwists[i - 1])
                + ad_mul_se(ad(&v), &self.local_screw[i - 1]) * dtheta[i - 1]
                + self.local_screw[i - 1].clone() * ddtheta[i - 1];
            dtwists[i] = dv;
        }

        let mut wrenchs = vec![se3::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]); N + 1];
        wrenchs[N] = f_tip;
        let mut tau = [0.0; N];
        for i in (1..N + 1).rev() {
            let t = (self.local_screw[i].clone() * -theta[i - 1]).exp()
                * self.local_zero_pose[i + 1].clone().inv();
            wrenchs[i] = t.adjoint().transpose().act(&wrenchs[i + 1])
                + inertial_mul_se3(self.local_spatial_inertial[i], &dtwists[i])
                - ad_mul_se(ad(&twists[i]).transpose(), &twists[i]);
            tau[i - 1] = cross_se3(&wrenchs[i], &self.local_screw[i - 1]);
        }
        tau
    }

    // M(q)^-1 * (tau - C(q, q_dot) * q_dot - G(q)) = q_dot_dot
    pub fn fd(&self) {
        println!("KidyChain fd");
    }
}

fn ad(v: &se3<f64>) -> Matrix6<f64> {
    let v = v.vee().as_array();
    let p = Vector3::new(v[0], v[1], v[2]);
    let w = Vector3::new(v[3], v[4], v[5]);
    let a = liealg::hat(&p);
    let b = liealg::hat(&w);
    let mut res = Matrix6::zeros();
    res.view_mut((0, 0), (3, 3)).copy_from(&a);
    res.view_mut((3, 3), (3, 3)).copy_from(&a);
    res.view_mut((3, 0), (3, 3)).copy_from(&b);
    res
}

fn ad_mul_se(m: Matrix6<f64>, w: &se3<f64>) -> se3<f64> {
    let vec6 = Vector6::from_column_slice(w.vee().as_slice());
    let res = m * vec6;
    se3::new([res[0], res[1], res[2]], [res[3], res[4], res[5]])
}

fn inertial_mul_se3(g: Matrix6<f64>, v: &se3<f64>) -> se3<f64> {
    let vec6 = Vector6::from_column_slice(v.vee().as_slice());
    let res = g * vec6;
    se3::new([res[0], res[1], res[2]], [res[3], res[4], res[5]])
}

fn cross_se3(v: &se3<f64>, w: &se3<f64>) -> f64 {
    let vec_v = Vector6::from_column_slice(v.vee().as_slice());
    let vec_w = Vector6::from_column_slice(w.vee().as_slice());
    vec_v.dot(&vec_w)
}

#[cfg(test)]
mod test {
    // use super::*;

    #[test]
    fn test_ad() {
        // let v = se3::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        // let res = ad(v);
        // println!("{}", res);
    }
}

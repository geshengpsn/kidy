use nalgebra::{Matrix4, Vector6};

// fk
fn fk(link_zero_pose: &Matrix4<f64>, twists_se3: &[Matrix4<f64>], twists_se3_sq: &[Matrix4<f64>], joints: &[f64])
{
    let init_pose = link_zero_pose;
    // joints.iter().zip(twists.iter()).skip(1).fold(, |pose, (q, xi)| {
    //     let exp_xi_q = exp_se3(xi, *q);
    //     pose * exp_xi_q
    // })
}


fn fk_diff<J>(joints: J) {}

fn fk_all<J>(joints: J) {}

fn fk_all_diff<J>(joints: J) {}

// ik
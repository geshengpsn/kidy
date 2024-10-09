mod jacobian;

use liealg::{prelude::*, Vec6, SE3};

use self::jacobian::{jacobian, Jacobian};

#[derive(Debug, PartialEq, Eq)]
pub enum FkError {}

pub fn fk<T: Real, const N: usize>(
    m: &SE3<T>,
    s_list: &[Vec6<T>; N],
    theta_list: &[T; N],
) -> Result<SE3<T>, FkError> {
    Ok(s_list
        .iter()
        .rev()
        .cloned()
        .zip(theta_list.iter().rev())
        .fold(m.clone(), |acc, (s, theta)| {
            (s * *theta).hat().exp().mat_mul(&acc)
        }))
}

#[derive(Debug, PartialEq, Eq)]
pub enum IkError {
    FKError(FkError),
    InvalidInput,
}

impl From<FkError> for IkError {
    fn from(e: FkError) -> Self {
        Self::FKError(e)
    }
}

#[derive(Debug)]
pub struct IkSolveParam<T: PartialEq + Real> {
    r_error: T,
    p_error: T,
    max_iter: usize,
}

impl<T: PartialEq + Real> Default for IkSolveParam<T> {
    fn default() -> Self {
        Self {
            r_error: T::epsilon(),
            p_error: T::epsilon(),
            max_iter: 100,
        }
    }
}

fn norm<T: Real>(v: &[T; 3]) -> T
// where for<'a> &'a T: core::ops::Mul<&'a T>
{
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

pub fn ik<T: Real + PartialEq, const N: usize>(
    m: &SE3<T>,
    s_list: &[Vec6<T>; N],
    theta_list: &mut [T; N],
    target: &SE3<T>,
    ik_solve_param: IkSolveParam<T>,
) -> Result<usize, IkError> {
    // f(θ') = f(θ) + J(θ)(θ'-θ)
    // f(θ') - f(θ) = J(θ)(θ'-θ)
    // θ' = θ + J⁻¹(θ)(f(θ') - f(θ))
    let mut i = 0;
    loop {
        // T_s_tip
        let tip = fk(m, s_list, theta_list)?;
        // Te = T_tip_target = T_s_tip⁻¹ * T_s_target
        let e = tip.inv().mat_mul(target);
        let e = e.log().vee();
        if norm(&e.p()) < ik_solve_param.p_error && norm(&e.r()) < ik_solve_param.r_error {
            return Ok(i);
        }
        let j = jacobian(s_list, theta_list);
        i += 1;
    }
}

#[cfg(test)]
mod test {
    use core::f64::consts::FRAC_PI_2;

    use super::*;
    use approx::*;
    use liealg::{Vec6, SO3};

    #[test]
    fn test_fk() {
        let m = (Vec6::new([0., 0., 0.], [1., 0., 0.]) * 2.).hat().exp();
        let s_list = [
            Vec6::new([0., 0., 1.], [0., 0., 0.]),
            Vec6::new([0., 0., 1.], [0., -1., 0.]),
        ];
        let theta_list = [0., FRAC_PI_2];
        let res = fk(&m, &s_list, &theta_list).unwrap();
        let m = SE3::new(&SO3::from_euler_angles(0., 0., FRAC_PI_2), [1., 1., 0.]);
        assert_abs_diff_eq!(m, res);
        println!("{:}", res);
    }

    #[test]
    fn test_ik() {
        let m = (Vec6::new([0., 0., 0.], [1., 0., 0.]) * 2.).hat().exp();
        let s_list = [
            Vec6::new([0., 0., 1.], [0., 0., 0.]),
            Vec6::new([0., 0., 1.], [0., -1., 0.]),
        ];
        let mut theta_list = [0., FRAC_PI_2];
        let target = SE3::new(&SO3::from_euler_angles(0., 0., FRAC_PI_2), [1., 1., 0.]);
        let res = ik(
            &m,
            &s_list,
            &mut theta_list,
            &target,
            IkSolveParam::default(),
        );
        // assert_eq!(res, Ok(()));
    }
}

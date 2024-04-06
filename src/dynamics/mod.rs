use liealg::{prelude::*, Vec6, SE3};

#[derive(Debug)]
pub enum FkError {
    InvalidInput,
}

pub fn fk<T>(m: &SE3<T>, s_list: &[Vec6<T>], theta_list: &[T]) -> Result<SE3<T>, FkError>
where
    T: Real,
{
    if s_list.len() != theta_list.len() {
        return Err(FkError::InvalidInput);
    }
    Ok(s_list
        .iter()
        .rev()
        .cloned()
        .zip(theta_list.iter().rev().copied())
        .fold(m.clone(), |acc, (s, theta)| {
            (s * theta).hat().exp().mat_mul(&acc)
        }))
}

#[cfg(test)]
mod test {
    use core::f64::consts::FRAC_PI_2;

    use super::*;
    use approx::*;
    use liealg::{Vec6, SO3};

    #[test]
    #[should_panic]
    fn test_fk_panic() {
        let m = (Vec6::new([0., 0., 0.], [1., 0., 0.]) * 2.).hat().exp();
        let s_list = [
            Vec6::new([0., 0., 1.], [0., 0., 0.]),
            Vec6::new([0., 0., 1.], [0., -1., 0.]),
        ];
        let theta_list = [0.];
        let _ = fk(&m, &s_list, &theta_list);
    }

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
        println!("{:.2}", res);
    }
}

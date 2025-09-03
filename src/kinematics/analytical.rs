#[cfg(feature = "realman75_6f_aik")]
pub mod realman75_6f;
// mod srsmodel;
pub mod zm75;

#[derive(Debug)]
pub enum AikError {
    Unreachable,
    // GradSingularity,
    // InvalidLimits,
}

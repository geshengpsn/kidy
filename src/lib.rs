//! # kidy is a library for kinematics and dynamics of multi-body.
//! kidy = kinematics + dynamics

// #![cfg_attr(not(test), no_std)]
// #![deny(missing_docs)]
#![deny(unsafe_code)]

mod model;
pub use model::*;

mod kinematics;
pub use kinematics::*;

mod dynamics;
pub mod visual;
// pub use liealg;

mod algo;
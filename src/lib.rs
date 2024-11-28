//! # kidy is a library for kinematics and dynamics of multi-body.
//! kidy = kinematics + dynamics

// #![cfg_attr(not(test), no_std)]
// #![deny(missing_docs)]
#![deny(unsafe_code)]

mod dynamics;
mod kinematics;
mod multi_body;

// re-export liealg
pub use liealg;
pub use multi_body::{KidyChain, Link, MultiBody};

// Entry point for the zkMatrix library
//
//! Zero-knowledge proof for linear algebra 
//!
#![doc = include_str!("../README.md")]
// extern crate sprs;

mod utils;
mod config;
mod test;

// pub mod util;
// pub mod protocols;
pub mod curve;
pub mod commit_mat;
pub mod dirac;
pub mod mat;
// pub mod protocols;
pub mod setup;

pub mod test_data;



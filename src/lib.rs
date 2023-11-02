// Entry point for the zkMatrix library
//
//! Zero-knowledge proof for linear algebra 
//!
#![doc = include_str!("../README.md")]

// extern crate sprs;

pub mod utils;
mod test;

// pub mod util;
// pub mod protocols;
pub mod curve;
pub mod dirac;
pub mod mat;
pub mod test_data;

mod config;
//! Entry point for the zkMatrix library
//!
//! Zero-knowledge proof for linear algebra 
//!
//! 
#![doc = include_str!("../README.md")]
// extern crate sprs;


pub mod commit_mat;
pub mod mat;
pub mod setup;

pub mod utils;
pub mod protocols;
pub mod zkprotocols;

pub mod config;
pub mod experiment_data;





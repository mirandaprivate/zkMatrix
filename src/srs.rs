//! Define the SRS used in the library.
//! 
//! 
use crate::curve::{G1Element, G2Element};

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
struct SRS {
    q: u64,
    g1_base_vec: Vec<G1Element>,
    g2_base_vec: Vec<G2Element>,
}
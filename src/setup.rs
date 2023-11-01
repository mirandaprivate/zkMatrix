//! Define the SRS used in the library.
//! 
//! 
use crate::curve::{G1Element, G2Element};

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
struct SRS {
    q: u64,
    g_hat_vec: Vec<G1Element>,
    h_hat_vec: Vec<G2Element>,
    g_hat_minus: G1Element,
    h_hat_minus: G2Element,
}
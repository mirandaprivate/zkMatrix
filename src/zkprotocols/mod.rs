//! Add zero-knowledge to the protocols
//! 
//! Since we target at logarithmic additional cost to add zero-knowledge
//! to the protocols, we use the following strategy:
//! 
//! 1. Associate a random Zp number to each transcript element.
//! 2. Blinding each transcript element by the random number.
//! 3. Replace the verification check of the original protocols by a group
//!   of zero-knowledge protocols on scalars, e.g. Schnorr's protocol.
//! 
pub mod zk_trans;
pub mod zk_scalars;

pub mod zk_ip_gt;
pub mod zk_scalar_proj;
pub mod zk_left_proj;
pub mod zk_right_proj;
pub mod zk_matmul;
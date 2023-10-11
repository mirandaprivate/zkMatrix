//! Defines the curve used in the library.
//! We use the BLS12-381 curve, which is a pairing-friendly curve.
//! 
//! We redefine the types here to make the code more readable.
//! 
//! Also, we implement the ops traits for these types by 
//!     wrapping the ops from the bls12_381 crate.
//! 
use std::ops::{Add, Mul, Neg, Sub};
use curv::BigInt;
use curv::arithmetic::traits::Converter;
use bls12_381::{Scalar, G1Projective, G2Projective, Gt};

pub type ZpElement = Scalar;
pub type IntElement = BigInt;
pub type G1Element = G1Projective;
pub type G2Element = G2Projective;
pub type GtElement = Gt;

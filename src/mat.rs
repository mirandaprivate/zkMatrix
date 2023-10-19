//! Define the type of the witness matrices used in the library.
//! 
//! # Example
//! 
use std::String;
use sprs::{CsMat, CsVec};

use crate::curve::{ZpElement, G1Element, G2Element, GtElement};

struct Mat {
    pub name: String,
    pub data_path: String,
    pub cache_path: String,
    pub data: CsMat<ZpElement>,
    pub g1_proj: CsVec<G1Element>,
    pub g2_proj: CsVec<G2Element>,
}
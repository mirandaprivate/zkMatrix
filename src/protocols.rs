pub mod hdsip;
pub mod mat_mul;
pub mod sip;

use crate::curve::{GtElement, ZpElement};

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
enum TransElement{
    ZpElement(ZpElement),
    GtElement(GtElement),
}

#[derive(Debug, PartialEq, Eq, Clone)]
struct Protocol {
    name: String,
    transcript: Vec<TransElement>,
}

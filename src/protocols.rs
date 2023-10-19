pub mod hdsip;
pub mod mat_mul;
pub mod sip;

use crate::curve::{GtElement, ZpElement};

enum TransElement{
    ZpElement(ZpElement),
    GtElement(GtElement),
}

struct Protocol {
    name: String,
    transcript: Vec<TransElement>,
}

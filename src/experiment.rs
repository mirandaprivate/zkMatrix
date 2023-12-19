use crate::utils::curve::{ZpElement, G1Element, G2Element, GtElement};
use crate::utils::to_file::FileIO;
use crate::setup::SRS;
use crate::mat::Mat;

use crate::config::{Q, DATA_DIR_PUBLIC, DATA_DIR_PRIVATE};


fn experiment(){
    let srs = SRS::new(Q);
    srs.to_file(String::from("srs"), true).unwrap();

    let srs_read: SRS = SRS::from_file(
        String::from("srs"), true
        ).unwrap();

    assert_eq!(srs, srs_read);
}

mod test{
    use super::*;
 
    #[test]
    fn experiment(){
        experiment();
    }
}
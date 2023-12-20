#![allow(dead_code)]

use std::time::Instant;

use crate::utils::curve::{ZpElement, G1Element, G2Element, GtElement};
use crate::utils::to_file::FileIO;
use crate::commit_mat::CommitMat;
use crate::setup::SRS;
use crate::mat::Mat;

use crate::experiment_data;
use crate::config::{Q, LOG_DIM, SQRT_MATRIX_DIM};


fn experiment_srs_gen(){
    let srs_timer = Instant::now();

    let srs = SRS::new(Q);

    let srs_duration = srs_timer.elapsed();

    println!(" ** SRS generation time: {:?}", srs_duration);

    srs.to_file(String::from("srs"), true).unwrap();

}

fn experiment_gen_matrices(){
    let mat_timer = Instant::now();

    let (c, a, b) = 
        experiment_data::gen_matrices_sparse(SQRT_MATRIX_DIM);

    let mat_duration = mat_timer.elapsed();

    println!(" ** Matrix generation time: {:?}", mat_duration);

    c.to_file(c.id.to_string(), false).unwrap();
    a.to_file(a.id.to_string(), false).unwrap();
    b.to_file(b.id.to_string(), false).unwrap();

}

fn experiment_commit_matrices(){

    let srs = SRS::from_file(String::from("srs"), true).unwrap();

    let a_read: Mat<ZpElement> = Mat::<ZpElement>::from_file(
        format!("a_sprs_2e{:?}", LOG_DIM), false
        ).unwrap();
    
    let commit_a_timer = Instant::now();

    let _a_commit = a_read.commit_rm(&srs);

    let commit_a_duration = commit_a_timer.elapsed();

    println!(" ** Commit matrix a time: {:?}", commit_a_duration);

    let c_read: Mat<ZpElement> = Mat::<ZpElement>::from_file(
        format!("b_sprs_2e{:?}", LOG_DIM), false
        ).unwrap();
    
    let commit_b_timer = Instant::now();

    let _b_commit = c_read.commit_cm(&srs);

    let commit_b_duration = commit_b_timer.elapsed();

    println!(" ** Commit matrix b time: {:?}", commit_b_duration);

    let c_read: Mat<ZpElement> = Mat::<ZpElement>::from_file(
        format!("c_sprs_2e{:?}", LOG_DIM), false
        ).unwrap();
    
    let commit_c_timer = Instant::now();

    let _c_commit = c_read.commit_cm(&srs);

    let commit_c_duration = commit_c_timer.elapsed();

    println!(" ** Commit matrix c time: {:?}", commit_c_duration);

}

mod test{
    use super::*;
 
    #[test]
    fn run_experiment(){
        experiment_srs_gen();
        experiment_gen_matrices();
        experiment_commit_matrices();
    }
}
//! # Main
//! 
//! Run experiment for the zkMatrix protocol.
//! 
//! The experiment parameters are configured in config.rs
//! 
//! Use the following command to run the experiment:
//! ''' cargo run --release '''
//! 
//! Note using the release mode is important for the performance.
//! 
#![allow(dead_code)]

use std::time::Instant;

use zkmatrix::commit_mat::CommitMat;
use zkmatrix::mat::Mat;
use zkmatrix::setup::SRS;

use zkmatrix::utils::curve::{G1Element, G2Element, GtElement};
use zkmatrix::utils::fiat_shamir::TranSeq;
use zkmatrix::utils::to_file::FileIO;

use zkmatrix::protocols::mat_mul::MatMul;

use zkmatrix::experiment_data;
use zkmatrix::config::{Q, LOG_DIM, SQRT_MATRIX_DIM};


fn main(){
    experiment_srs_gen();
    experiment_gen_matrices();
    experiment_commit_matrices();
    experiment_matmul();
}


fn experiment_srs_gen(){
    let srs_timer = Instant::now();

    let srs = SRS::new(Q);

    let srs_duration = srs_timer.elapsed();

    println!(" ** SRS generation time: {:?}", srs_duration);

    srs.to_file(String::from("srs.srs"), true).unwrap();

}

fn experiment_gen_matrices(){
    let mat_timer = Instant::now();

    let (c, a, b) = 
        experiment_data::gen_matrices_sparse(SQRT_MATRIX_DIM);

    let mat_duration = mat_timer.elapsed();

    println!(" ** Matrix generation time: {:?}", mat_duration);

    c.to_file(format!("{}.dat", c.id), false).unwrap();
    a.to_file(format!("{}.dat", a.id), false).unwrap();
    b.to_file(format!("{}.dat", b.id), false).unwrap();

}


fn experiment_commit_matrices(){

    let srs = SRS::from_file(String::from("srs.srs"), true).unwrap();

    let a_read: Mat<i64> = Mat::<i64>::from_file(
        format!("a_sprs_i64_2e{:?}.dat", LOG_DIM), false
        ).unwrap();
    
    let commit_a_timer = Instant::now();

    let a_commit = a_read.commit_rm(&srs);

    let commit_a_duration = commit_a_timer.elapsed();

    println!(" ** Commit matrix a time: {:?}", commit_a_duration);

    let b_read: Mat<i64> = Mat::<i64>::from_file(
        format!("b_sprs_i64_2e{:?}.dat", LOG_DIM), false
        ).unwrap();
    
    let commit_b_timer = Instant::now();

    let b_commit = b_read.commit_cm(&srs);

    let commit_b_duration = commit_b_timer.elapsed();

    println!(" ** Commit matrix b time: {:?}", commit_b_duration);

    let c_read: Mat<i128> = Mat::<i128>::from_file(
        format!("c_sprs_i128_2e{:?}.dat", LOG_DIM), false
        ).unwrap();
    
    let commit_c_timer = Instant::now();

    let c_commit = c_read.commit_cm(&srs);

    let commit_c_duration = commit_c_timer.elapsed();

    println!(" ** Commit matrix c time: {:?}", commit_c_duration);

    let commit_cab = [c_commit, a_commit, b_commit].to_vec();

    commit_cab.to_file(
        format!("commit_abc_sprs_2e{:?}.com", LOG_DIM), 
        true).unwrap()


}


fn experiment_matmul() {

    let srs = SRS::from_file(String::from("srs.srs"), true).unwrap();

    let a_read: Mat<i64> = Mat::<i64>::from_file(
        format!("a_sprs_i64_2e{:?}.dat", LOG_DIM), false
        ).unwrap();

    let b_read: Mat<i64> = Mat::<i64>::from_file(
        format!("b_sprs_i64_2e{:?}.dat", LOG_DIM), false
        ).unwrap();
    
    
    let c_read: Mat<i128> = Mat::<i128>::from_file(
        format!("c_sprs_i128_2e{:?}.dat", LOG_DIM), false
        ).unwrap();

    let a_chache_read: Vec<G2Element> = FileIO::from_file(
        format!("{}_rp.cache", a_read.id), false
        ).unwrap();

    let b_chache_read: Vec<G1Element> = FileIO::from_file(
        format!("{}_lp.cache", b_read.id), false
        ).unwrap();
    
    let c_chache_read: Vec<G1Element> = FileIO::from_file(
        format!("{}_lp.cache", c_read.id), false
        ).unwrap();

    let commit_cab: Vec<GtElement> = FileIO::from_file(
        format!("commit_abc_sprs_2e{:?}.com", LOG_DIM), 
        true
        ).unwrap();

    let timer_prove = Instant::now();
    
    let matmul_protocol = MatMul::new(
        commit_cab[0],
        commit_cab[1],
        commit_cab[2], 
        SQRT_MATRIX_DIM * SQRT_MATRIX_DIM,
        SQRT_MATRIX_DIM * SQRT_MATRIX_DIM,
        SQRT_MATRIX_DIM * SQRT_MATRIX_DIM,
    );
    
    let mut trans = TranSeq::new();

    matmul_protocol.prove::<i128, i64, i64>(
        &srs,
        &mut trans,
        &c_read, &a_read, &b_read,
        &c_chache_read, &a_chache_read, &b_chache_read, 
    );

    srs.to_file(format!("tr_2e{:?}.tr", LOG_DIM), true).unwrap();

    println!(" ** Prover time of MatMul: {:?}", timer_prove.elapsed());

    let timer_verify = Instant::now();

    let result = matmul_protocol.verify(&srs, &mut trans);

    println!(" * Verification of MatMul result: {:?}", result);

    println!(" ** Verifier time of MatMul : {:?}", timer_verify.elapsed());

}


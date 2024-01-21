//! # Main
//! 
//! Run experiments for the zkMatrix protocol.
//! 
//! The experiment parameters are configured in config.rs
//! 
//! Use the following command to run the experiment:
//! ''' cargo bench --release '''
//! 
//! Note using the release mode is important for the performance.
//! 
#![allow(dead_code)]
use std::time::Instant;

use std::fs::OpenOptions;
use std::fs::File;
use std::io::Write;

use zkmatrix::commit_mat::CommitMat;
use zkmatrix::mat::Mat;
use zkmatrix::setup::SRS;

use zkmatrix::utils::curve::ZpElement;
use zkmatrix::utils::curve::{G1Element, G2Element, GtElement};
use zkmatrix::utils::fiat_shamir::TranSeq;
use zkmatrix::utils::to_file::FileIO;

use zkmatrix::protocols::mat_mul::MatMul;

use zkmatrix::zkprotocols::zk_matmul::ZkMatMul;

use zkmatrix::experiment_data;
use zkmatrix::config::{Q, LOG_DIM, SQRT_MATRIX_DIM};

use zkmatrix::zkprotocols::zk_trans::ZkTranSeqProver;

fn main(){
    let mut log_file = OpenOptions::new()
        .create(true)
        .write(true)
        .append(true)
        .open(format!("log/log_2e{:?}", LOG_DIM)).unwrap();

    // experiment_srs_gen(&mut log_file);
    // experiment_gen_matrices(&mut log_file);
    // experiment_commit_matrices(&mut log_file);
    // experiment_matmul(&mut log_file);
    experiment(&mut log_file);
    // experiment_dense(&mut log_file);
}


fn experiment(log_file: &mut File) {

    for t in 2..(LOG_DIM/2+1) {

        let log_dim = 2*t;
        let sqrt_dim = 2usize.pow(t as u32);
        let matrix_dim = 2usize.pow(log_dim as u32);
        let q: usize = 2usize.pow(log_dim as u32) + 2;

        println!(" ** Experiment for zkMatrix, Matrix Dim 2e{:?} times 2e{:?}; Number of non-zero elements: 2e{:?} ** ",
            log_dim,
            log_dim,
            log_dim / 2 * 3,
        );

        let srs_timer = Instant::now();

        let srs = SRS::new(q);

        let srs_duration = srs_timer.elapsed();

        println!(" ** SRS generation time: {:?}", srs_duration);
        writeln!(log_file, " ** SRS generation time: {:?}", srs_duration).unwrap();

        let mat_timer = Instant::now();

        let (c, a, b) = 
            experiment_data::gen_matrices_sparse(sqrt_dim);

        let mat_duration = mat_timer.elapsed();

        println!(" ** Matrix generation time: {:?}", mat_duration);
        writeln!(log_file, " ** Matrix generation time: {:?}", mat_duration).unwrap();

        let commit_a_timer = Instant::now();

        let a_tilde = ZpElement::rand();
        let (a_com, a_cache) = 
            a.commit_rm(&srs);
        let a_blind = a_com + a_tilde * srs.blind_base;

        let commit_a_duration = commit_a_timer.elapsed();

        println!(" ** Commit matrix a time: {:?}", commit_a_duration);
        writeln!(log_file, " ** Commit matrix a time: {:?}", commit_a_duration).unwrap();
        
        let commit_b_timer = Instant::now();

        let b_tilde = ZpElement::rand();
        let (b_com, b_cache) = b.commit_cm(&srs);
        let b_blind = b_com + b_tilde * srs.blind_base;

        let commit_b_duration = commit_b_timer.elapsed();

        println!(" ** Commit matrix b time: {:?}", commit_b_duration);
        writeln!(log_file, " ** Commit matrix b time: {:?}", commit_b_duration).unwrap();

        let commit_c_timer = Instant::now();

        let c_tilde = ZpElement::rand();
        let (c_com, c_cache) = c.commit_cm(&srs);
        let c_blind = c_com + c_tilde * srs.blind_base;

        let commit_c_duration = commit_c_timer.elapsed();

        println!(" ** Commit matrix c time: {:?}", commit_c_duration);
        writeln!(log_file, " ** Commit matrix c time: {:?}", commit_c_duration).unwrap();

        let commit_cab: Vec<GtElement> = [c_blind, a_blind, b_blind].to_vec();

        let timer_prove = Instant::now();
        
        let matmul_protocol = ZkMatMul::new(
            commit_cab[0],
            commit_cab[1],
            commit_cab[2], 
            matrix_dim,
            matrix_dim,
            matrix_dim,
        );
        
        let mut zk_trans = ZkTranSeqProver::new(&srs);

        matmul_protocol.prove::<i128, i64, i64>(
            &srs,
            &mut zk_trans,
            c, a, b,
            &c_cache, &a_cache, &b_cache, 
            c_tilde, a_tilde, b_tilde,
        );

        let trans = zk_trans.publish_trans();

        println!(" ** Prover time of zkMatMul: {:?}", timer_prove.elapsed());
        writeln!(log_file, " ** Prover time of zkMatMul: {:?}", 
            timer_prove.elapsed()
        ).unwrap();

        trans.save_to_file(
            format!("tr_2e{:?}", log_dim)
        ).unwrap();

        let mut trans_read = TranSeq::read_from_file(
            format!("tr_2e{:?}", log_dim)
        );

        let result = matmul_protocol.verify(&srs, &mut trans_read);

        println!(" * Verification of zkMatMul result: {:?}", result);


        let mut verify_time: f64 = 0.;
        let repeat = 10;
        
        for _ in 0..repeat {
            let mut trans_read = TranSeq::read_from_file(
                format!("tr_2e{:?}", log_dim)
            );

            let timer_verify = Instant::now();
            matmul_protocol.verify(&srs, &mut trans_read);
            verify_time += timer_verify.elapsed().as_secs_f64()/repeat as f64;
        }

        verify_time = verify_time * 1000.;

        println!(" ** Verifier time of zkMatMul : {:?}ms", verify_time);
        writeln!(log_file, " ** Verifier time of zkMatMul : {:?}ms", 
            verify_time)
        .unwrap();

    }

}

fn experiment_dense(log_file: &mut File) {

    println!(" ** Experiment for zkMatrix, Matrix Dim 2e{:?} times 2e{:?}; Number of non-zero elements: 2e{:?} ** ",
    LOG_DIM/2,
    LOG_DIM/2,
    LOG_DIM,
    );

    let srs_timer = Instant::now();

    let srs = SRS::new(2usize.pow(LOG_DIM as u32/2) + 2);

    let srs_duration = srs_timer.elapsed();

    println!(" ** SRS generation time: {:?}", srs_duration);
    writeln!(log_file, " ** SRS generation time: {:?}", srs_duration).unwrap();

    let mat_timer = Instant::now();

    let (c, a, b) = 
        experiment_data::gen_matrices_dense(SQRT_MATRIX_DIM);
    
    c.to_file(c.id.to_string(), false).unwrap();
    a.to_file(a.id.to_string(), false).unwrap();
    b.to_file(b.id.to_string(), false).unwrap();

    let mat_duration = mat_timer.elapsed();

    println!(" ** Matrix generation time: {:?}", mat_duration);
    writeln!(log_file, " ** Matrix generation time: {:?}", mat_duration).unwrap();

    let commit_a_timer = Instant::now();

    let a_tilde = ZpElement::rand();
    let (a_com, a_cache) = 
        a.commit_rm(&srs);
    let a_blind = a_com + a_tilde * srs.blind_base;

    let commit_a_duration = commit_a_timer.elapsed();

    println!(" ** Commit matrix a time: {:?}", commit_a_duration);
    writeln!(log_file, " ** Commit matrix a time: {:?}", commit_a_duration).unwrap();
    
    let commit_b_timer = Instant::now();

    let b_tilde = ZpElement::rand();
    let (b_com, b_cache) = b.commit_cm(&srs);
    let b_blind = b_com + b_tilde * srs.blind_base;

    let commit_b_duration = commit_b_timer.elapsed();

    println!(" ** Commit matrix b time: {:?}", commit_b_duration);
    writeln!(log_file, " ** Commit matrix b time: {:?}", commit_b_duration).unwrap();

    let commit_c_timer = Instant::now();

    let c_tilde = ZpElement::rand();
    let (c_com, c_cache) = c.commit_cm(&srs);
    let c_blind = c_com + c_tilde * srs.blind_base;

    let commit_c_duration = commit_c_timer.elapsed();

    println!(" ** Commit matrix c time: {:?}", commit_c_duration);
    writeln!(log_file, " ** Commit matrix c time: {:?}", commit_c_duration).unwrap();

    let commit_cab: Vec<GtElement> = [c_blind, a_blind, b_blind].to_vec();

    let timer_prove = Instant::now();
    
    let matmul_protocol = ZkMatMul::new(
        commit_cab[0],
        commit_cab[1],
        commit_cab[2], 
        SQRT_MATRIX_DIM,
        SQRT_MATRIX_DIM,
        SQRT_MATRIX_DIM,
    );
    
    let mut zk_trans = ZkTranSeqProver::new(&srs);

    matmul_protocol.prove::<i128, i64, i64>(
        &srs,
        &mut zk_trans,
        c, a, b,
        &c_cache, &a_cache, &b_cache, 
        c_tilde, a_tilde, b_tilde,
    );

    let trans = zk_trans.publish_trans();

    println!(" ** Prover time of zkMatMul: {:?}", timer_prove.elapsed());
    writeln!(log_file, " ** Prover time of zkMatMul: {:?}", 
        timer_prove.elapsed()
    ).unwrap();

    trans.save_to_file(
        format!("tr_2e{:?}", LOG_DIM/2)
    ).unwrap();

    let mut trans_read = TranSeq::read_from_file(
        format!("tr_2e{:?}", LOG_DIM/2)
    );

    let timer_verify = Instant::now();

    let mut verify_time: f64 = 0.;
    let repeat = 10;
        
    for _ in 0..repeat {
        let mut trans_read = TranSeq::read_from_file(
            format!("tr_2e{:?}", LOG_DIM/2)
        );

        let timer_verify = Instant::now();
        matmul_protocol.verify(&srs, &mut trans_read);
        verify_time += timer_verify.elapsed().as_secs_f64()/repeat as f64;
    }

    verify_time = verify_time * 1000.;

    let result = matmul_protocol.verify(&srs, &mut trans_read);

    println!(" * Verification of zkMatMul result: {:?}", result);

    println!(" ** Verifier time of zkMatMul : {:?}", verify_time);
    writeln!(log_file, " ** Verifier time of zkMatMul : {:?}", 
        timer_verify.elapsed())
    .unwrap();

}

fn experiment_srs_gen(log_file: &mut File){
    let srs_timer = Instant::now();

    let srs = SRS::new(Q);

    let srs_duration = srs_timer.elapsed();

    println!(" ** SRS generation time: {:?}", srs_duration);
    writeln!(log_file, " ** SRS generation time: {:?}", srs_duration).unwrap();
    

    srs.to_file(
        format!("srs_2e{:?}.srs", LOG_DIM),
        true,
    ).unwrap();

}

fn experiment_gen_matrices(log_file: &mut File){
    let mat_timer = Instant::now();

    let (c, a, b) = 
        experiment_data::gen_matrices_sparse_from_kronecker(SQRT_MATRIX_DIM);

    let mat_duration = mat_timer.elapsed();

    println!(" ** Matrix generation time: {:?}", mat_duration);
    writeln!(log_file, " ** Matrix generation time: {:?}", mat_duration).unwrap();

    c.to_file(format!("{}.dat", c.id), false).unwrap();
    a.to_file(format!("{}.dat", a.id), false).unwrap();
    b.to_file(format!("{}.dat", b.id), false).unwrap();

}


fn experiment_commit_matrices(log_file: &mut File){

    let srs = SRS::from_file(format!(
        "srs_2e{:?}.srs", LOG_DIM), true
    ).unwrap();

    let timer_read = Instant::now();

    let a_read: Mat<i64> = Mat::<i64>::from_file(
        format!("a_sprs_i64_2e{:?}.dat", LOG_DIM), false
        ).unwrap();

    println!(" ** Read matrix a time: {:?}", timer_read.elapsed());
    
    let commit_a_timer = Instant::now();

    let (a_commit,_) = a_read.commit_rm(&srs);

    let commit_a_duration = commit_a_timer.elapsed();

    println!(" ** Commit matrix a time: {:?}", commit_a_duration);
    writeln!(log_file, " ** Commit matrix a time: {:?}", commit_a_duration).unwrap();

    let b_read: Mat<i64> = Mat::<i64>::from_file(
        format!("b_sprs_i64_2e{:?}.dat", LOG_DIM), false
        ).unwrap();
    
    let commit_b_timer = Instant::now();

    let (b_commit,_) = b_read.commit_cm(&srs);

    let commit_b_duration = commit_b_timer.elapsed();

    println!(" ** Commit matrix b time: {:?}", commit_b_duration);
    writeln!(log_file, " ** Commit matrix b time: {:?}", commit_b_duration).unwrap();

    let c_read: Mat<i128> = Mat::<i128>::from_file(
        format!("c_sprs_i128_2e{:?}.dat", LOG_DIM), false
        ).unwrap();
    
    let commit_c_timer = Instant::now();

    let (c_commit,_) = c_read.commit_cm(&srs);

    let commit_c_duration = commit_c_timer.elapsed();

    println!(" ** Commit matrix c time: {:?}", commit_c_duration);
    writeln!(log_file, " ** Commit matrix c time: {:?}", commit_c_duration).unwrap();

    let commit_cab = [c_commit, a_commit, b_commit].to_vec();

    commit_cab.to_file(
        format!("commit_cab_sprs_2e{:?}.com", LOG_DIM), 
        true).unwrap()


}


fn experiment_matmul(log_file: &mut File) {

    let srs = SRS::from_file(
        format!("srs_2e{:?}.srs", LOG_DIM), 
        true,
    ).unwrap();

    let a_read: Mat<i64> = Mat::<i64>::from_file(
        format!("a_sprs_i64_2e{:?}.dat", LOG_DIM), false
        ).unwrap();

    let b_read: Mat<i64> = Mat::<i64>::from_file(
        format!("b_sprs_i64_2e{:?}.dat", LOG_DIM), false
        ).unwrap();
    
    
    let c_read: Mat<i128> = Mat::<i128>::from_file(
        format!("c_sprs_i128_2e{:?}.dat", LOG_DIM), false
        ).unwrap();

    let a_cache_read: Vec<G2Element> = FileIO::from_file(
        format!("{}_rp.cache", a_read.id), false
        ).unwrap();

    let b_cache_read: Vec<G1Element> = FileIO::from_file(
        format!("{}_lp.cache", b_read.id), false
        ).unwrap();
    
    let c_cache_read: Vec<G1Element> = FileIO::from_file(
        format!("{}_lp.cache", c_read.id), false
        ).unwrap();

    let commit_cab: Vec<GtElement> = FileIO::from_file(
        format!("commit_cab_sprs_2e{:?}.com", LOG_DIM), 
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
        c_read, a_read, b_read,
        &c_cache_read, &a_cache_read, &b_cache_read, 
    );


    println!(" ** Prover time of MatMul: {:?}", timer_prove.elapsed());
    writeln!(log_file, " ** Prover time of MatMul: {:?}", timer_prove.elapsed()).unwrap();

    trans.save_to_file(
        format!("tr_2e{:?}", LOG_DIM)
    ).unwrap();

    let mut trans_read = TranSeq::read_from_file(
        format!("tr_2e{:?}", LOG_DIM)
    );

    let timer_verify = Instant::now();

    let result = matmul_protocol.verify(&srs, &mut trans_read);

    println!(" * Verification of MatMul result: {:?}", result);

    println!(" ** Verifier time of MatMul : {:?}", timer_verify.elapsed());
    writeln!(log_file, " ** Verifier time of MatMul : {:?}", timer_verify.elapsed()).unwrap();

}



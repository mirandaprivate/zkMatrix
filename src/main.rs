use std::time::Instant;

use zkmatrix::commit_mat::CommitMat;
use zkmatrix::mat::Mat;
use zkmatrix::setup::SRS;

use zkmatrix::utils::curve::ZpElement;
use zkmatrix::utils::to_file::FileIO;

use zkmatrix::utils::dirac;

use zkmatrix::experiment_data;
use zkmatrix::config::{Q, LOG_DIM, SQRT_MATRIX_DIM};


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

fn experiment_braket(){


    let srs = SRS::from_file(String::from("srs"), true).unwrap();

    let a_left = experiment_data::gen_mat_rand_diag_i64(
        SQRT_MATRIX_DIM,26);
    let a_right = experiment_data::gen_mat_rand_dense_i64(
        SQRT_MATRIX_DIM, 26);
    let a: Mat<i64> = Mat::new_from_data_vec(
        "a_i64_sprs", 
        (SQRT_MATRIX_DIM*SQRT_MATRIX_DIM, SQRT_MATRIX_DIM*SQRT_MATRIX_DIM), 
        experiment_data::diag_kronecker_dense_from_i64_to_i64(
            &a_left, &a_right),
    );

    let a_zp: Mat<ZpElement> = Mat::new_from_data_vec(
        "a_zp_sprs", 
        (SQRT_MATRIX_DIM*SQRT_MATRIX_DIM, SQRT_MATRIX_DIM*SQRT_MATRIX_DIM), 
        experiment_data::diag_kronecker_dense_from_i64_to_zp(
            &a_left, &a_right),
    );

    let timer_bra_opt = Instant::now();
    
    let result_bra_opt = dirac::bra_opt_i64(&a, &srs.g_hat_vec);
    
    println!(" ** Bra opt time: {:?}", timer_bra_opt.elapsed());

    let timer_bra_no_opt = Instant::now();

    let result_bra_no_opt = dirac::proj_left(
        &a_zp, &srs.g_hat_vec);

    println!(" ** Bra no opt time: {:?}", timer_bra_no_opt.elapsed());

    assert_eq!(result_bra_opt, result_bra_no_opt);
    println!("Assert Equal of two method");
    
    let timer_ket_opt = Instant::now();
    
    let result_ket_opt = dirac::ket_opt_i64(
        &a, &srs.h_hat_vec);
    
    println!(" ** Ket opt time: {:?}", timer_ket_opt.elapsed());

    let timer_ket_no_opt = Instant::now();

    let result_ket_no_opt = dirac::proj_right(
        &a_zp, &srs.h_hat_vec);

    println!(" ** Ket no opt time: {:?}", timer_ket_no_opt.elapsed());

    assert_eq!(result_ket_opt, result_ket_no_opt);
    println!("Assert Equal of two method");
    
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


fn main(){
    // experiment_srs_gen();
    // experiment_gen_matrices();
    // experiment_commit_matrices();
    experiment_braket();
}

//! Test data for unit tests and integration tests.
//! 
//! We use kronecker product of a sparse matrix and a dense matrix to generate
//!  the sparse matrices.
//!
//! The dimension of the resulting matrix is:
//!     (SQRT_MATRIX_DIM * SQRT_MATRIX_DIM, SQRT_MATRIX_DIM * SQRT_MATRIX_DIM)
//! 

#![allow(dead_code)]

use std::ops::{Add, Mul};

use curv::arithmetic::Zero;
use rand::Rng;

use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

use crate::mat::Mat;

use crate::utils::curve::ZpElement;
use crate::utils::to_file::FileIO;




pub fn gen_mat_rand_dense_u64(sqrt_dim: usize) -> Vec<Vec<u64>>{
    (0..sqrt_dim).map(|_| {
        (0..sqrt_dim).map(|_| {
            rand::thread_rng().gen_range(0..2u64.pow(26)) as u64
        }).collect::<Vec<u64>>()  
    }).collect::<Vec<Vec<u64>>>()
}

pub fn gen_mat_rand_diag_u64(sqrt_dim: usize) -> Vec<u64>{
    (0..sqrt_dim).map(|_| {
        rand::thread_rng().gen_range(0..2u64.pow(26)) as u64
    }).collect::<Vec<u64>>()  
}

fn mat_mul_dense_u64_to_zp(a: &Vec<Vec<u64>>, b: &Vec<Vec<u64>>) 
    -> Vec<Vec<ZpElement>> {
    let a_rows = a.len();
    let a_cols = a[0].len();
    let b_cols = b[0].len();

    let mut result = 
        vec![vec![ZpElement::zero(); b_cols]; a_rows];

    let pool = ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();

    pool.install(|| {
        result.par_iter_mut().enumerate().for_each(|(i, row)| {
            for j in 0..b_cols {
                for k in 0..a_cols {
                    row[j] += ZpElement::from(a[i][k]) * ZpElement::from(b[k][j]);
                }
            }
        });
    });

    result
}

fn mat_mul_dense_u64_to_u64(a: &Vec<Vec<u64>>, b: &Vec<Vec<u64>>) 
    -> Vec<Vec<u64>> {
    let a_rows = a.len();
    let a_cols = a[0].len();
    let b_cols = b[0].len();

    let mut result = 
        vec![vec![0 as u64; b_cols]; a_rows];

    let pool = ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();

    pool.install(|| {
        result.par_iter_mut().enumerate().for_each(|(i, row)| {
            for j in 0..b_cols {
                for k in 0..a_cols {
                    row[j] += a[i][k] * b[k][j];
                }
            }
        });
    });

    result
}

fn mat_mul_dense_zp_to_zp(a: &Vec<Vec<ZpElement>>, b: &Vec<Vec<ZpElement>>) 
    -> Vec<Vec<ZpElement>> {
    let a_rows = a.len();
    let a_cols = a[0].len();
    let b_cols = b[0].len();

    let mut result = vec![vec![ZpElement::zero(); b_cols]; a_rows];

    for i in 0..a_rows {
        for j in 0..b_cols {
            for k in 0..a_cols {
                result[i][j] += 
                    a[i][k] * b[k][j];
            }
        }
    }

    result
}


fn mat_mul_diag_u64_to_zp(a: &Vec<u64>, b: &Vec<u64>) -> Vec<ZpElement> {
    
    a.iter().zip(b.iter()).map(|(left, right)| {
        ZpElement::from(*left) *  ZpElement::from(*right)
    }).collect()
}

fn mat_mul_diag_u64_to_u64(a: &Vec<u64>, b: &Vec<u64>) -> Vec<u64> {
    
    a.iter().zip(b.iter()).map(|(left, right)| {
        left * right
    }).collect()
}

fn diag_kronecker_dense_from_u64(a: &Vec<u64>, b: &Vec<Vec<u64>>) 
    -> Vec<(usize, usize, ZpElement)> {
    
    let sqrt_dim = a.len();

    let pool = ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();

    pool.install(|| {
        (0..sqrt_dim).into_par_iter().flat_map(|left_ij|{
            (0..sqrt_dim).flat_map( |right_i|{
                (0..sqrt_dim).map( move |right_j|{
                    (
                        left_ij * sqrt_dim + right_i,
                        left_ij * sqrt_dim + right_j,
                        ZpElement::from(a[left_ij]) * ZpElement::from(b[right_i][right_j]) 
                    )
                })
            }).collect::<Vec<(usize, usize, ZpElement)>>()
        }).collect::<Vec<(usize, usize, ZpElement)>>()
    })
    
}

fn diag_kronecker_dense_from_zp(a: &Vec<ZpElement>, b: &Vec<Vec<ZpElement>>) 
    -> Vec<(usize, usize, ZpElement)> {
    
    let sqrt_dim = a.len();

    let pool = ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();

    pool.install(|| {
        (0..sqrt_dim).into_par_iter().flat_map(|left_ij|{
            (0..sqrt_dim).flat_map( |right_i|{
                (0..sqrt_dim).map( move |right_j|{
                    (
                        left_ij * sqrt_dim + right_i,
                        left_ij * sqrt_dim + right_j,
                        a[left_ij] * b[right_i][right_j] 
                    )
                })
            }).collect::<Vec<(usize, usize, ZpElement)>>()
        }).collect::<Vec<(usize, usize, ZpElement)>>()
    })
}


fn gen_matrices_sparse(sqrt_dim: usize) -> (Mat<ZpElement>, Mat<ZpElement>, Mat<ZpElement>) {
    
    let dim = sqrt_dim * sqrt_dim;
    let log_dim = (dim as u64).ilog2() as usize;


    let a_left = gen_mat_rand_diag_u64(sqrt_dim);
    let a_right = gen_mat_rand_dense_u64(sqrt_dim);

    let b_left = gen_mat_rand_diag_u64(sqrt_dim);
    let b_right = gen_mat_rand_dense_u64(sqrt_dim);

    let c_left = mat_mul_diag_u64_to_u64(&a_left, &b_left);
    let c_right = mat_mul_dense_u64_to_u64(&a_right, &b_right);

    let a = Mat::new_from_data_vec(
        "a",
        (sqrt_dim * sqrt_dim, sqrt_dim * sqrt_dim), 
        diag_kronecker_dense_from_u64(&a_left, &a_right)
    );

    let b = Mat::new_from_data_vec(
        "b",
        (sqrt_dim * sqrt_dim, sqrt_dim * sqrt_dim), 
        diag_kronecker_dense_from_u64(&b_left, &b_right)
    );

    let c = Mat::new_from_data_vec(
        "c",
        (sqrt_dim * sqrt_dim, sqrt_dim * sqrt_dim), 
        diag_kronecker_dense_from_u64(&c_left, &c_right)
    );

    a.to_file(format!("a_sprs_2e{:?}", log_dim), false).unwrap();
    b.to_file(format!("b_sprs_2e{:?}", log_dim), false).unwrap();
    c.to_file(format!("c_sprs_2e{:?}", log_dim), false).unwrap();

    (a, b, c)
}

fn sprs_to_dense_zp(sprs: &Mat<ZpElement>) -> Vec<Vec<ZpElement>> {
    let mut dense = 
        vec![vec![ZpElement::zero(); sprs.shape.1]; sprs.shape.0];

    for &(row, col, ref val) in &sprs.data {
        dense[row][col] = *val;
    }

    dense
}

fn dense_to_sprs_zp(id: &str, dense: &Vec<Vec<ZpElement>>) -> Mat<ZpElement> {
    let mut sprs = Mat::new(
        id,
        (dense.len(), dense[0].len())
    );

    for i in 0..dense.len() {
        for j in 0..dense[0].len() {
            if dense[i][j] != ZpElement::zero() {
                sprs.push(i, j, dense[i][j]);
            }
        }
    }

    sprs
}

fn dense_to_sprs_zp_from_u64(id: &str, dense: &Vec<Vec<u64>>) -> Mat<ZpElement> {
    let mut sprs = Mat::new(
        id,
        (dense.len(), dense[0].len())
    );

    for i in 0..dense.len() {
        for j in 0..dense[0].len() {
            if dense[i][j] != 0 {
                sprs.push(i, j, ZpElement::from(dense[i][j]));
            }
        }
    }

    sprs
}

fn gen_matrices_dense(dim: usize) -> (Mat<ZpElement>, Mat<ZpElement>, Mat<ZpElement>) {
    
    let log_dim = (dim as u64).ilog2() as usize;

    let a_u64 = gen_mat_rand_dense_u64(dim);

    let b_u64 = gen_mat_rand_dense_u64(dim);

    let c_zp = mat_mul_dense_u64_to_zp(&a_u64, &b_u64);

    let a = dense_to_sprs_zp_from_u64("a", &a_u64);

    let b = dense_to_sprs_zp_from_u64("b", &b_u64);

    let c = dense_to_sprs_zp("c", &c_zp);

    a.to_file(format!("a_dense_2e{:?}", log_dim), false).unwrap();
    b.to_file(format!("b_dense_2e{:?}", log_dim), false).unwrap();
    c.to_file(format!("c_dense_2e{:?}", log_dim), false).unwrap();

    (a, b, c)
}


#[cfg(test)]
mod tests{

    use super::*;

    use crate::config::SQRT_MATRIX_DIM_TEST;

    #[test]
    fn test_experiment_data(){
        let (a, b, c) = 
            gen_matrices_sparse(SQRT_MATRIX_DIM_TEST);

        let (a_d, b_d, c_d) = 
            gen_matrices_dense(SQRT_MATRIX_DIM_TEST*SQRT_MATRIX_DIM_TEST);

        let a_d_dense = sprs_to_dense_zp(&a_d);
        let b_d_dense = sprs_to_dense_zp(&b_d);
        let c_d_dense = sprs_to_dense_zp(&c_d);

        assert_eq!(mat_mul_dense_zp_to_zp(&a_d_dense, &b_d_dense), c_d_dense);

        // let a_dense = sprs_to_dense_zp(&a);
        // let b_dense = sprs_to_dense_zp(&b);
        // let c_dense = sprs_to_dense_zp(&c);

        // assert_eq!(mat_mul_dense_zp_to_zp(&a_dense, &b_dense), c_dense);
    }
}
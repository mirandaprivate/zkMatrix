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

    for i in 0..a_rows {
        for j in 0..b_cols {
            for k in 0..a_cols {
                result[i][j] += 
                    ZpElement::from(a[i][k] * b[k][j]);
            }
        }
    }

    result
}

fn mat_mul_dense_u64_to_u64(a: &Vec<Vec<u64>>, b: &Vec<Vec<u64>>) 
    -> Vec<Vec<u64>> {
    let a_rows = a.len();
    let a_cols = a[0].len();
    let b_cols = b[0].len();

    let mut result = 
        vec![vec![0 as u64; b_cols]; a_rows];

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
        ZpElement::from(*left *  *right)
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

    (0..sqrt_dim).flat_map(|left_ij|{
        (0..sqrt_dim).flat_map( move |right_i|{
            (0..sqrt_dim).map( move |right_j|{
                (
                    left_ij * sqrt_dim + right_i,
                    left_ij * sqrt_dim + right_j,
                    ZpElement::from(a[left_ij] * b[right_i][right_j]) 
                )
            })
        })
    }).collect()
}

fn diag_kronecker_dense_from_zp(a: &Vec<ZpElement>, b: &Vec<Vec<ZpElement>>) 
    -> Vec<(usize, usize, ZpElement)> {
    
    let sqrt_dim = a.len();

    (0..sqrt_dim).flat_map(|left_ij|{
        (0..sqrt_dim).flat_map( move |right_i|{
            (0..sqrt_dim).map( move |right_j|{
                (
                    left_ij * sqrt_dim + right_i,
                    left_ij * sqrt_dim + right_j,
                    a[left_ij] * b[right_i][right_j] 
                )
            })
        })
    }).collect()
}

fn gen_matrices(sqrt_dim: usize) -> (Mat<ZpElement>, Mat<ZpElement>, Mat<ZpElement>) {
    
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

    a.to_file(String::from("a"), false).unwrap();
    b.to_file(String::from("b"), false).unwrap();
    c.to_file(String::from("c"), false).unwrap();

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

#[cfg(test)]
mod tests{

    use super::*;

    use crate::config::SQRT_MATRIX_DIM_TEST;

    #[test]
    fn test_experiment_data(){
        let (a, b, c) = 
            gen_matrices(SQRT_MATRIX_DIM_TEST);
        let a_dense = sprs_to_dense_zp(&a);
        let b_dense = sprs_to_dense_zp(&b);
        let c_dense = sprs_to_dense_zp(&c);

        // assert_eq!(mat_mul_dense_zp_to_zp(&a_dense, &b_dense), c_dense);
    }
}
//! Test data for unit tests and integration tests.
//! 
//! We use kronecker product of a sparse matrix and a dense matrix to generate
//!  the sparse matrices.
//!
//! The dimension of the resulting matrix is:
//!     (SQRT_MATRIX_DIM * SQRT_MATRIX_DIM, SQRT_MATRIX_DIM * SQRT_MATRIX_DIM)
//! 
//! The product of the test matrices are manually computed for testing.
//! 

#![allow(dead_code)]
use rand::Rng;

use crate::mat::Mat;
use crate::config::SQRT_MATRIX_DIM_TEST;

use crate::utils::curve::{ZpElement, G1Element, G2Element, GtElement};
use crate::utils::kronecker::{kronecker, kronecker_vec};

#[cfg(test)]
pub fn gen_test_matrices() -> (Mat<i128>, Mat<i64>, Mat<i64>){
    crate::experiment_data::gen_matrices_dense(
        crate::config::MATRIX_DIM_TEST
    )
}

pub fn dense_to_sprs_i64(a: &Vec<Vec<i64>>) -> Vec<(usize, usize, i64)> {
    let mut sprs = Vec::new();

    for i in 0..a.len() {
        for j in 0..a[0].len() {
            if a[i][j] != 0 {
                sprs.push((i, j, a[i][j]));
            }
        }
    }

    sprs
}

pub fn dense_to_sprs_i128(a: &Vec<Vec<i128>>) -> Vec<(usize, usize, i128)> {
    let mut sprs = Vec::new();

    for i in 0..a.len() {
        for j in 0..a[0].len() {
            if a[i][j] != 0 {
                sprs.push((i, j, a[i][j]));
            }
        }
    }

    sprs
}

pub fn gen_test_matrices_not_square() -> (Mat<i128>, Mat<i64>, Mat<i64>){

    let a_data = (0..32).map(|_| {
        (0..16).map(|_| {
            rand::thread_rng()
            .gen_range(-2i64.pow(52)..2i64.pow(52)) as i64
        }).collect::<Vec<i64>>()  
    }).collect::<Vec<Vec<i64>>>();

    let b_data = (0..16).map(|_| {
        (0..8).map(|_| {
            rand::thread_rng()
            .gen_range(-2i64.pow(52)..2i64.pow(52)) as i64
        }).collect::<Vec<i64>>()  
    }).collect::<Vec<Vec<i64>>>();

    let c_data = 
        crate::experiment_data::mat_mul_dense_i64_to_i128(&a_data, &b_data);

    let a = Mat::new_from_data_vec(
        "a_no_square", (32, 16), dense_to_sprs_i64(&a_data)
    );

    let b = Mat::new_from_data_vec(
        "b_no_square", (16, 8), dense_to_sprs_i64(&b_data)
    );

    let c = Mat::new_from_data_vec(
        "c_no_square", (32, 8), dense_to_sprs_i128(&c_data)
    );

    (c, a, b)

}

#[cfg(test)]
pub fn gen_vec_v_from_kronecker_test() -> Vec<ZpElement>{

    let vec_v_part: Vec<ZpElement> = (0..SQRT_MATRIX_DIM_TEST).map(|i|
        ZpElement::from((i+1) as u64)
    ).collect();

    kronecker_vec(&vec_v_part, &vec_v_part)
} 

#[cfg(test)]
pub fn gen_vec_v_direct_test() -> Vec<ZpElement> {
    let vec_v: Vec<ZpElement> = 
        (0..SQRT_MATRIX_DIM_TEST).flat_map(|i|{
            (0..SQRT_MATRIX_DIM_TEST).map(move|j|{
                    ZpElement::from(((i+1)*(j+1)) as u64)
                })
        }).collect();
    
    vec_v
}

#[cfg(test)]
pub fn gen_mat_a_zp_from_kronecker_test() -> Mat<ZpElement>{
    
    let mat_a_left: Mat<ZpElement> = Mat::new_from_data_vec(
        "mat_a_left_test", 
        (SQRT_MATRIX_DIM_TEST, SQRT_MATRIX_DIM_TEST), 
        (0..SQRT_MATRIX_DIM_TEST).map(|i|
            (i, i, ZpElement::from((i+1) as u64))
        ).collect()
    );
    
    let mat_a_right: Mat<ZpElement> = Mat::new_from_data_vec(
        "mat_a_right_test",
        (SQRT_MATRIX_DIM_TEST, SQRT_MATRIX_DIM_TEST), 
        (0..SQRT_MATRIX_DIM_TEST).flat_map(|i|
            (0..SQRT_MATRIX_DIM_TEST).map( move |j|
                (i, j, ZpElement::from(((i+1)*(j+1)) as u64))
            )
        ).collect()
    );

    let mat_a = kronecker(&mat_a_left, &mat_a_right);
    println!("mat_a: {:?}", mat_a);

    mat_a
}

#[cfg(test)]
pub fn gen_mat_a_zp_direct_test() -> Mat<ZpElement>{
    let m: Vec<(usize, usize, ZpElement)> =
    (0..SQRT_MATRIX_DIM_TEST).flat_map(|left_ij|{
        (0..SQRT_MATRIX_DIM_TEST).flat_map(move|right_i|{
            (0..SQRT_MATRIX_DIM_TEST).map(move|right_j|(
                left_ij * SQRT_MATRIX_DIM_TEST + right_i,
                left_ij * SQRT_MATRIX_DIM_TEST + right_j,
                ZpElement::from(
                    (left_ij as u64 + 1) * (right_i as u64 + 1) * (right_j as u64 + 1)
                )
            ))
        })
    }).collect();

    let mat_dim = SQRT_MATRIX_DIM_TEST * SQRT_MATRIX_DIM_TEST;

    Mat::new_from_data_vec("a", (mat_dim, mat_dim), m)
}

#[cfg(test)]
pub fn gen_mat_a_u64_direct_test() -> Mat<u64>{
    let m: Vec<(usize, usize, u64)> =
    (0..SQRT_MATRIX_DIM_TEST).flat_map(|left_ij|{
        (0..SQRT_MATRIX_DIM_TEST).flat_map(move|right_i|{
            (0..SQRT_MATRIX_DIM_TEST).map(move|right_j|(
                left_ij * SQRT_MATRIX_DIM_TEST + right_i,
                left_ij * SQRT_MATRIX_DIM_TEST + right_j,
                (left_ij as u64 + 1) * (right_i as u64 + 1) * (right_j as u64 + 1)
            ))
        })
    }).collect();

    let mat_dim = SQRT_MATRIX_DIM_TEST * SQRT_MATRIX_DIM_TEST;

    Mat::new_from_data_vec("a_test", (mat_dim, mat_dim), m)
}

pub fn gen_vec_va_from_kronecker_test() -> Vec<ZpElement>{
    let sum_n_square: usize = (1..=SQRT_MATRIX_DIM_TEST)
        .map(|i| i*i).sum();

    let vec_va_left: Vec<ZpElement> = (0..SQRT_MATRIX_DIM_TEST)
        .map(|i|
        ZpElement::from(((i+1)*(i+1)) as u64)
        ).collect();

    let vec_va_right: Vec<ZpElement> = (0..SQRT_MATRIX_DIM_TEST)
        .map(|i|
         ZpElement::from( (sum_n_square * (i+1)) as u64)
        ).collect();

    kronecker_vec(&vec_va_left, &vec_va_right)
}

fn gen_mat_a_gt_from_kronecker_test() -> Mat<GtElement>{
    let mat_a_left: Mat<G1Element> = Mat::new_from_data_vec(
        "mat_a_left_test", 
        (SQRT_MATRIX_DIM_TEST, SQRT_MATRIX_DIM_TEST), 
        (0..SQRT_MATRIX_DIM_TEST).map(|i|
            (i, i, G1Element::from((i+1) as u64))
        ).collect()
    );
    
    let mat_a_right: Mat<G2Element> = Mat::new_from_data_vec(
        "mat_a_right_test",
        (SQRT_MATRIX_DIM_TEST, SQRT_MATRIX_DIM_TEST), 
        (0..SQRT_MATRIX_DIM_TEST).flat_map(|i|
            (0..SQRT_MATRIX_DIM_TEST).map( move |j|
                (i, j, G2Element::from(((i+1)*(j+1)) as u64))
            )
        ).collect()
    );

    let mat_a: Mat<GtElement> = 
        kronecker(&mat_a_left, &mat_a_right);

    mat_a
}

#[cfg(test)]
fn gen_mat_a_gt_direct_test() -> Mat<GtElement>{
    let m: Vec<(usize, usize, u64)> =
    (0..SQRT_MATRIX_DIM_TEST).flat_map(|left_ij|{
        (0..SQRT_MATRIX_DIM_TEST).flat_map(move|right_i|{
            (0..SQRT_MATRIX_DIM_TEST).map(move|right_j|(
                left_ij * SQRT_MATRIX_DIM_TEST + right_i,
                left_ij * SQRT_MATRIX_DIM_TEST + right_j,
                ((left_ij as u64 + 1) * (right_i as u64 + 1) * (
                    right_j as u64 + 1))
            ))
        })
    }).collect();

    let m_gt = m.into_iter()
        .map(|(row, col, val)|{
        (row, col, GtElement::from(val))
        }).collect();

    let mat_dim = SQRT_MATRIX_DIM_TEST * SQRT_MATRIX_DIM_TEST;

    Mat::new_from_data_vec("a_test", (mat_dim, mat_dim), m_gt)
}

#[cfg(test)]
fn gen_vec_v_gt_from_kronecker_test()-> Vec<GtElement>{

    let vec_v_left: Vec<G1Element> = (0..SQRT_MATRIX_DIM_TEST)
        .map(|i|
            G1Element::from((i+1) as u64)
        ).collect();

    let vec_v_right: Vec<G2Element> = (0..SQRT_MATRIX_DIM_TEST)
        .map(|i|
            G2Element::from((i+1) as u64)
        ).collect();

    kronecker_vec(&vec_v_left, &vec_v_right)
} 

#[cfg(test)]
pub fn gen_vec_v_gt_direct_test() -> Vec<GtElement> {
    let vec_v: Vec<GtElement> = 
        (0..SQRT_MATRIX_DIM_TEST).flat_map(|i|{
            (0..SQRT_MATRIX_DIM_TEST).map(move|j|{
                    GtElement::from(((i+1)*(j+1)) as u64)
                })
        }).collect();
    
    vec_v
}

#[cfg(test)]
pub fn gen_vec_v_g1_direct_test() -> Vec<G1Element> {
    let vec_v: Vec<G1Element> = 
        (0..SQRT_MATRIX_DIM_TEST).flat_map(|i|{
            (0..SQRT_MATRIX_DIM_TEST).map(move|j|{
                    G1Element::from(((i+1)*(j+1)) as u64)
            })
        }).collect();
    
    vec_v
}

#[cfg(test)]
pub fn gen_vec_v_g2_direct_test() -> Vec<G2Element> {
    let vec_v: Vec<G2Element> = 
        (0..SQRT_MATRIX_DIM_TEST).flat_map(|i|{
            (0..SQRT_MATRIX_DIM_TEST).map(move|j|{
                    G2Element::from(((i+1)*(j+1)) as u64)
                })
        }).collect();
    
    vec_v
}

#[cfg(test)]
pub fn gen_vec_va_gt_from_kronecker_test() -> Vec<GtElement>{
    let sum_n_square: usize = (1..=SQRT_MATRIX_DIM_TEST)
        .map(|i| i*i).sum();

    let vec_va_left: Vec<G1Element> = (0..SQRT_MATRIX_DIM_TEST)
        .map(|i|
            G1Element::from(((i+1)*(i+1)) as u64)
        ).collect();

    let vec_va_right: Vec<G2Element> = (0..SQRT_MATRIX_DIM_TEST)
        .map(|i|
             G2Element::from( (sum_n_square * (i+1)) as u64)
        ).collect();

    kronecker_vec(&vec_va_left, &vec_va_right)
}

#[cfg(test)]
pub fn gen_vec_va_g1_from_kronecker_test() -> Vec<G1Element>{
    let sum_n_square: usize = (1..=SQRT_MATRIX_DIM_TEST)
        .map(|i| i*i).sum();

    let vec_va_left: Vec<ZpElement> = (0..SQRT_MATRIX_DIM_TEST)
        .map(|i|
            ZpElement::from(((i+1)*(i+1)) as u64)
        ).collect();

    let vec_va_right: Vec<G1Element> = (0..SQRT_MATRIX_DIM_TEST)
        .map(|i|
             G1Element::from( (sum_n_square * (i+1)) as u64)
        ).collect();

    kronecker_vec(&vec_va_left, &vec_va_right)
}

#[cfg(test)]
pub fn gen_vec_va_g2_from_kronecker_test() -> Vec<G2Element>{
    let sum_n_square: usize = (1..=SQRT_MATRIX_DIM_TEST)
        .map(|i| i*i).sum();

    let vec_va_left: Vec<ZpElement> = (0..SQRT_MATRIX_DIM_TEST)
        .map(|i|
            ZpElement::from(((i+1)*(i+1)) as u64)
        ).collect();

    let vec_va_right: Vec<G2Element> = (0..SQRT_MATRIX_DIM_TEST)
        .map(|i|
             G2Element::from( (sum_n_square * (i+1)) as u64)
        ).collect();

    kronecker_vec(&vec_va_left, &vec_va_right)
}

#[cfg(test)]
pub fn gen_vav_direct_test() -> ZpElement {
    let sum_n_square: usize = (1..=SQRT_MATRIX_DIM_TEST)
        .map(|i| i*i).sum();
    let sum_n_cubic: usize = (1..=SQRT_MATRIX_DIM_TEST)
        .map(|i| i*i*i).sum();

    ZpElement::from((sum_n_cubic * sum_n_square * sum_n_square) as u64)
}

#[cfg(test)]
pub fn gen_vav_gt_direct_test() -> GtElement{
    let sum_n_square: usize = (1..=SQRT_MATRIX_DIM_TEST)
        .map(|i| i*i).sum();
    let sum_n_cubic: usize = (1..=SQRT_MATRIX_DIM_TEST)
        .map(|i| i*i*i).sum();

    GtElement::from((sum_n_cubic * sum_n_square * sum_n_square) as u64)
}

#[cfg(test)]
mod tests{

    use super::*;

    #[test]
    fn test_gen_mat(){
        assert_eq!(
            gen_mat_a_zp_from_kronecker_test(), gen_mat_a_zp_direct_test()
        );
        assert_eq!(
            gen_vec_v_from_kronecker_test(), gen_vec_v_direct_test()
        );
        assert_eq!(
            gen_mat_a_gt_from_kronecker_test(), gen_mat_a_gt_direct_test()
        );
        assert_eq!(
            gen_vec_v_gt_from_kronecker_test(), gen_vec_v_gt_direct_test()
        );
    }
}
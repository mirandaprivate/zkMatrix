//! Test data for unit tests and integration tests.
//! 
//! We use kronecker product of a sparse matrix and a dense matrix to generate
//!  the sparse matrices.
//!
//! The dimension of the resulting matrix is:
//!     (SQRT_MATRIX_DIM * SQRT_MATRIX_DIM, SQRT_MATRIX_DIM * SQRT_MATRIX_DIM)
//! 
use std::ops::{Add, Mul};

use crate::curve::{ZpElement, G1Element, G2Element, GtElement};
use crate::mat::Mat;

const SQRT_MATRIX_DIM: usize = 2usize.pow(3);

pub fn kronecker<T, U, V>(mat_a: &Mat<T>, mat_b: &Mat<U>) -> Mat<V>
where
    T: Clone + Mul<U, Output =V>,
    U: Clone,
    V: Clone + Add + From<u64>,
{

    let result_id = format!("{}{}{}", mat_a.id, "_kronecker_", mat_b.id);
    let result_shape = (mat_a.shape.0 * mat_b.shape.0, mat_a.shape.1 * mat_b.shape.1);
    let mut result = Mat::new(&result_id, result_shape);

    for &(row_a, col_a, ref val_a) in &mat_a.data {
        for &(row_b, col_b, ref val_b) in &mat_b.data {
            let row = row_a * mat_b.shape.0 + row_b;
            let col = col_a * mat_b.shape.1 + col_b;
            let val = val_a.clone() * val_b.clone();
            result.push(row, col, val);
        }
    }

    result
}

pub fn kronecker_vec<T, U, V>(vec_a: &Vec<T>, vec_b: &Vec<U>) -> Vec<V>
where
    T: Clone + Mul<U, Output =V>,
    U: Clone,
    V: Clone + Add + From<u64>,
{
    let result: Vec<V> =
    (0..vec_a.len()).flat_map(|i|{
        (0..vec_b.len()).map(move|j|{
            vec_a[i].clone() * vec_b[j].clone()
        })
    }).collect();

    result
}

pub fn gen_vec_v_from_kronecker() -> Vec<ZpElement>{

    let vec_v_part: Vec<ZpElement> = (0..SQRT_MATRIX_DIM).map(|i|
        ZpElement::from((i+1) as u64)
    ).collect();

    kronecker_vec(&vec_v_part, &vec_v_part)
} 

pub fn gen_vec_v_direct() -> Vec<ZpElement> {
    let vec_v: Vec<ZpElement> = 
        (0..SQRT_MATRIX_DIM).flat_map(|i|{
            (0..SQRT_MATRIX_DIM).map(move|j|{
                    ZpElement::from(((i+1)*(j+1)) as u64)
                })
        }).collect();
    
    vec_v
}

pub fn gen_mat_a_from_kronecker() -> Mat<ZpElement>{
    
    let mat_a_left: Mat<ZpElement> = Mat::new_from_data_vec(
        "mat_a_left", 
        (SQRT_MATRIX_DIM, SQRT_MATRIX_DIM), 
        (0..SQRT_MATRIX_DIM).map(|i|
            (i, i, ZpElement::from((i+1) as u64))
        ).collect()
    );
    
    let mat_a_right: Mat<ZpElement> = Mat::new_from_data_vec(
        "mat_a_right",
        (SQRT_MATRIX_DIM, SQRT_MATRIX_DIM), 
        (0..SQRT_MATRIX_DIM).flat_map(|i|
            (0..SQRT_MATRIX_DIM).map( move |j|
                (i, j, ZpElement::from(((i+1)*(j+1)) as u64))
            )
        ).collect()
    );

    let mat_a = kronecker(&mat_a_left, &mat_a_right);
    println!("mat_a: {:?}", mat_a);

    mat_a
}

pub fn gen_mat_a_direct() -> Mat<ZpElement>{
    let m: Vec<(usize, usize, ZpElement)> =
    (0..SQRT_MATRIX_DIM).flat_map(|left_ij|{
        (0..SQRT_MATRIX_DIM).flat_map(move|right_i|{
            (0..SQRT_MATRIX_DIM).map(move|right_j|(
                left_ij * SQRT_MATRIX_DIM + right_i,
                left_ij * SQRT_MATRIX_DIM + right_j,
                ZpElement::from((left_ij as u64 + 1) * (right_i as u64 + 1) * (right_j as u64 + 1))
            ))
        })
    }).collect();

    let mat_dim = SQRT_MATRIX_DIM * SQRT_MATRIX_DIM;

    Mat::new_from_data_vec("a", (mat_dim, mat_dim), m)
}

pub fn gen_vec_va_from_kronecker() -> Vec<ZpElement>{
    let sum_n_square: usize = (1..=SQRT_MATRIX_DIM).map(|i| i*i).sum();

    let vec_va_left: Vec<ZpElement> = (0..SQRT_MATRIX_DIM).map(|i|
        ZpElement::from(((i+1)*(i+1)) as u64)
    ).collect();

    let vec_va_right: Vec<ZpElement> = (0..SQRT_MATRIX_DIM).map(|i|
         ZpElement::from( (sum_n_square * (i+1)) as u64)
    ).collect();

    kronecker_vec(&vec_va_left, &vec_va_right)
}

fn gen_mat_a_gt_from_kronecker() -> Mat<GtElement>{
    let mat_a_left: Mat<G1Element> = Mat::new_from_u64_vec(
        "mat_a_left", 
        (SQRT_MATRIX_DIM, SQRT_MATRIX_DIM), 
        (0..SQRT_MATRIX_DIM).map(|i|
            (i, i, ((i+1) as u64))
        ).collect()
    );
    
    let mat_a_right: Mat<G2Element> = Mat::new_from_u64_vec(
        "mat_a_right",
        (SQRT_MATRIX_DIM, SQRT_MATRIX_DIM), 
        (0..SQRT_MATRIX_DIM).flat_map(|i|
            (0..SQRT_MATRIX_DIM).map( move |j|
                (i, j, ((i+1)*(j+1)) as u64)
            )
        ).collect()
    );

    let mat_a: Mat<GtElement> = kronecker(&mat_a_left, &mat_a_right);

    mat_a
}

fn gen_mat_a_gt_direct() -> Mat<GtElement>{
    let m: Vec<(usize, usize, u64)> =
    (0..SQRT_MATRIX_DIM).flat_map(|left_ij|{
        (0..SQRT_MATRIX_DIM).flat_map(move|right_i|{
            (0..SQRT_MATRIX_DIM).map(move|right_j|(
                left_ij * SQRT_MATRIX_DIM + right_i,
                left_ij * SQRT_MATRIX_DIM + right_j,
                ((left_ij as u64 + 1) * (right_i as u64 + 1) * (right_j as u64 + 1))
            ))
        })
    }).collect();

    let mat_dim = SQRT_MATRIX_DIM * SQRT_MATRIX_DIM;

    Mat::new_from_u64_vec("a", (mat_dim, mat_dim), m)
}

fn gen_vec_v_gt_from_kronecker()-> Vec<GtElement>{

    let vec_v_left: Vec<G1Element> = (0..SQRT_MATRIX_DIM).map(|i|
        G1Element::from((i+1) as u64)
    ).collect();

    let vec_v_right: Vec<G2Element> = (0..SQRT_MATRIX_DIM).map(|i|
        G2Element::from((i+1) as u64)
    ).collect();

    kronecker_vec(&vec_v_left, &vec_v_right)
} 

pub fn gen_vec_v_gt_direct() -> Vec<GtElement> {
    let vec_v: Vec<GtElement> = 
        (0..SQRT_MATRIX_DIM).flat_map(|i|{
            (0..SQRT_MATRIX_DIM).map(move|j|{
                    GtElement::from(((i+1)*(j+1)) as u64)
                })
        }).collect();
    
    vec_v
}

pub fn gen_vec_v_g1_direct() -> Vec<G1Element> {
    let vec_v: Vec<G1Element> = 
        (0..SQRT_MATRIX_DIM).flat_map(|i|{
            (0..SQRT_MATRIX_DIM).map(move|j|{
                    G1Element::from(((i+1)*(j+1)) as u64)
                })
        }).collect();
    
    vec_v
}

pub fn gen_vec_v_g2_direct() -> Vec<G2Element> {
    let vec_v: Vec<G2Element> = 
        (0..SQRT_MATRIX_DIM).flat_map(|i|{
            (0..SQRT_MATRIX_DIM).map(move|j|{
                    G2Element::from(((i+1)*(j+1)) as u64)
                })
        }).collect();
    
    vec_v
}

pub fn gen_vec_va_gt_from_kronecker() -> Vec<GtElement>{
    let sum_n_square: usize = (1..=SQRT_MATRIX_DIM).map(|i| i*i).sum();

    let vec_va_left: Vec<G1Element> = (0..SQRT_MATRIX_DIM).map(|i|
        G1Element::from(((i+1)*(i+1)) as u64)
    ).collect();

    let vec_va_right: Vec<G2Element> = (0..SQRT_MATRIX_DIM).map(|i|
         G2Element::from( (sum_n_square * (i+1)) as u64)
    ).collect();

    kronecker_vec(&vec_va_left, &vec_va_right)
}

pub fn gen_vec_va_g1_from_kronecker() -> Vec<G1Element>{
    let sum_n_square: usize = (1..=SQRT_MATRIX_DIM).map(|i| i*i).sum();

    let vec_va_left: Vec<ZpElement> = (0..SQRT_MATRIX_DIM).map(|i|
        ZpElement::from(((i+1)*(i+1)) as u64)
    ).collect();

    let vec_va_right: Vec<G1Element> = (0..SQRT_MATRIX_DIM).map(|i|
         G1Element::from( (sum_n_square * (i+1)) as u64)
    ).collect();

    kronecker_vec(&vec_va_left, &vec_va_right)
}

pub fn gen_vav_direct() -> ZpElement {
    let sum_n_square: usize = (1..=SQRT_MATRIX_DIM).map(|i| i*i).sum();
    let sum_n_cubic: usize = (1..=SQRT_MATRIX_DIM).map(|i| i*i*i).sum();

    ZpElement::from((sum_n_cubic * sum_n_square * sum_n_square) as u64)
}

pub fn gen_vav_gt_direct() -> GtElement{
    let sum_n_square: usize = (1..=SQRT_MATRIX_DIM).map(|i| i*i).sum();
    let sum_n_cubic: usize = (1..=SQRT_MATRIX_DIM).map(|i| i*i*i).sum();

    GtElement::from((sum_n_cubic * sum_n_square * sum_n_square) as u64)
}

#[cfg(test)]
mod tests{

    use super::*;

    #[test]
    fn test_gen_mat(){
        assert_eq!(gen_mat_a_from_kronecker(), gen_mat_a_direct());
        assert_eq!(gen_vec_v_from_kronecker(), gen_vec_v_direct());
        assert_eq!(gen_mat_a_gt_from_kronecker(), gen_mat_a_gt_direct());
        assert_eq!(gen_vec_v_gt_from_kronecker(), gen_vec_v_gt_direct());
    }
}
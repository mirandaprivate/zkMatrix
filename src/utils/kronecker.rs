use std::ops::{Add, Mul};

use crate::mat::Mat;

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
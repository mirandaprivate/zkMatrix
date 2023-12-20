use core::fmt::Debug;
use std::marker::{Send, Sync};
use std::sync::{Arc, Mutex};
use std::thread;

use rayon::{ThreadPoolBuilder, prelude::*};


use crate::utils::curve::{Add, Mul, AddAssign, Zero, G1Element, G2Element, GtElement};

use crate::mat::Mat;

use crate::config::NUM_THREADS;


pub fn inner_product<T, U, V>(vec_a: &Vec<T>, vec_b: &Vec<U>) -> V
where
    T: 'static + Clone + Copy + Send + Sync + Mul<U, Output =V>,
    U: 'static + Clone + Copy + Send + Sync,
    V: 'static + Clone + Copy + Send + Sync + Add + Zero,
{
    let length = vec_a.len();

    if length != vec_b.len() {
        panic!("Length of two vectors are not equal!");
    }

    if length <= NUM_THREADS {
        return vec_a.iter()
        .zip(vec_b.iter())
        .fold(
            V::zero(),
            |acc, (&a_val, &b_val)| acc + a_val * b_val,
        );
    }

    let pool = ThreadPoolBuilder::new()
        .num_threads(NUM_THREADS)
        .build()
        .unwrap();

    let inner_product = pool.install(|| {
        vec_a.par_iter()
        .zip(vec_b.par_iter())
        .fold(
            || V::zero(),
            |acc, (&a_val, &b_val)| acc + a_val * b_val,
        )
        .reduce(
            || V::zero(),
            |acc, val| acc + val,
        )
    });

    return inner_product;

}


pub fn proj_left<T, U, V>(mat_a: &Mat<T>, v_base: &Vec<U>) -> Vec<V>
where
    T: 'static + Clone + Copy + Send + Sync + Mul<U, Output = V>,
    U: 'static + Clone + Copy + Send + Sync + Mul<T, Output = V>,
    V: 'static + Clone + Copy + Send + Sync + Debug + Add + AddAssign + Zero,
{
    let a_data  = &mat_a.data;
    let v_base_arc = Arc::new(v_base.clone());
    let v_len = v_base.len();

    let result = Arc::new(Mutex::new(vec![V::zero(); v_len]));

    let chunks: Vec<_> = a_data.chunks(a_data.len() / NUM_THREADS).collect();

    let mut handles = Vec::new();
    for chunk in chunks {
        let chunk = chunk.to_owned();
        let v_base_clone = Arc::clone(&v_base_arc);
        let handle = thread::spawn(
            move || {
                let mut local_result = vec![V::zero(); v_len];
                for &(row, col, val) in &chunk {
                    local_result[col] += v_base_clone[row] * val;
                }
                local_result
            });
        handles.push(handle);
    }

    for handle in handles {
        let thread_result = handle.join().unwrap();
        let mut result = result.lock().unwrap();
        for i in 0..v_len {
            result[i] += thread_result[i];
        }
    }

    Arc::try_unwrap(result).unwrap().into_inner().unwrap()

    // // The following is slower than the current version
    // // Both version works.
    // let a_data  = &mat_a.data;
    // let v_base_arc = Arc::new(v_base.clone());
    
    // let num_nonzero_entries = a_data.len();
    // let v_len = v_base.len();

    // let result = a_data.chunks(
    //     num_nonzero_entries / NUM_THREADS)
    //     .map(|chunk| {
    //         let v_arc_clone = Arc::clone(&v_base_arc);
    //         let chunk = chunk.to_owned();
    //         thread::spawn(
    //             move || 
    //             chunk.iter().fold(
    //                 vec![V::zero(); v_len],
    //                 |mut acc, &(row, col, val)| {
    //                     acc[col] += v_arc_clone[row] * val;
    //                     acc
    //                 }))
    //     })
    //     .fold(
    //         vec![V::zero(); v_len],
    //         |acc, handle| 
    //         handle.join().unwrap().iter().enumerate().map(
    //             |(i, &l)| l + acc[i]
    //         ).collect()
    //     );
    // result
}

pub fn proj_right<T, U, V>(mat_a: &Mat<T>, v_base: &Vec<U>) -> Vec<V>
where
    T: 'static + Clone + Copy + Send + Sync + Mul<U, Output = V>,
    U: 'static + Clone + Copy + Send + Sync + Mul<T, Output = V>,
    V: 'static + Clone + Copy + Send + Sync + Debug + Add + AddAssign + Zero,
{
    let a_data  = &mat_a.data;
    let v_base_arc = Arc::new(v_base.clone());
    let v_len = v_base.len();

    let result = Arc::new(Mutex::new(vec![V::zero(); v_len]));

    let chunks: Vec<_> = a_data.chunks(a_data.len() / NUM_THREADS).collect();

    let mut handles = Vec::new();
    for chunk in chunks {
        let chunk = chunk.to_owned();
        let v_base_clone = Arc::clone(&v_base_arc);
        let handle = thread::spawn(
            move || {
                let mut local_result = vec![V::zero(); v_len];
                for &(row, col, val) in &chunk {
                    local_result[row] += v_base_clone[col] * val;
                }
                local_result
            });
        handles.push(handle);
    }

    for handle in handles {
        let thread_result = handle.join().unwrap();
        let mut result = result.lock().unwrap();
        for i in 0..v_len {
            result[i] += thread_result[i];
        }
    }

    Arc::try_unwrap(result).unwrap().into_inner().unwrap()
}

pub fn dirac<T, U, V, R, W>(
    g_base: &Vec<T>,
    mat_a: &Mat<U>, 
    h_base: &Vec<R>,
) -> W
where
    T: 'static + Clone + Copy + Send + Sync + Mul<U, Output = V>,
    U: 'static + Clone + Copy + Send + Sync + Mul<T, Output = V>,
    V: 'static + Clone + Copy + Send + Sync + Debug + Add + AddAssign + Zero
        + Mul<R, Output = W>,
    R: 'static + Clone + Copy + Send + Sync,
    W: 'static + Clone + Copy + Send + Sync + Add + Zero,
{
    let proj_left = proj_left(mat_a, g_base);
    let result = inner_product(&proj_left, &h_base);
    
    result
}


pub fn bra(mat_a: &Mat<u64>, v_base: &Vec<G1Element>) -> Vec<G1Element> {
    proj_left(&mat_a, &v_base)
}

pub fn ket(mat_a: &Mat<u64>, v_base: &Vec<G2Element>) -> Vec<G2Element> {
    proj_right(&mat_a, &v_base)
}

pub fn braket(
    g_base: &Vec<G1Element>,
    mat_a: &Mat<u64>, 
    h_base: &Vec<G2Element>,
) -> GtElement {
    dirac(&g_base, &mat_a, &h_base)
}


pub fn bra_opt(mat_a: &Mat<u64>, v_base: &Vec<G1Element>) -> Vec<G1Element> {
    proj_left(&mat_a, &v_base)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_data::*;

    use crate::utils::curve::{G1Element, G2Element, GtElement};

    #[test]
    fn test_dirac() {
        // ga = a left-muliplied by g 
        let mat_a: Mat<u64> = gen_mat_a_u64_direct_test();
        let vec_g: Vec<G1Element> = gen_vec_v_g1_direct_test();
        let vec_ga: Vec<G1Element> = gen_vec_va_g1_from_kronecker_test();
        let vec_h: Vec<G2Element> = gen_vec_v_g2_direct_test();
        let gah: GtElement = gen_vav_gt_direct_test();

        let proj_test_left: Vec<G1Element> = proj_left(&mat_a, &vec_g);
        assert_eq!(proj_test_left, vec_ga);
        let gah_test_1 = inner_product(&proj_test_left, &vec_h);
        let bra_test = bra(&mat_a, &vec_g);
        assert_eq!(bra_test, vec_ga);

        let proj_test_right: Vec<G2Element> = proj_right(&mat_a, &vec_h);
        let gah_test_2 = inner_product(&vec_g, &proj_test_right);
        let ket_test = ket(&mat_a, &vec_h);
        assert_eq!(ket_test, proj_test_right);
        
        let gah_test_3 = dirac(&vec_g, &mat_a, &vec_h);
        let gah_test_4 = braket(&vec_g, &mat_a, &vec_h);

        assert_eq!(gah_test_1, gah);
        assert_eq!(gah_test_2, gah);
        assert_eq!(gah_test_3, gah);
        assert_eq!(gah_test_4, gah);

    }
}
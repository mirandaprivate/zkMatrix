use std::ops::{Add, Mul, AddAssign};

use bls12_381::{pairing, G1Affine, G1Projective, G2Affine, G2Projective, Gt, Scalar};
use curv::arithmetic::Zero;
use rand::Rng;
use std::sync::{Arc, Mutex};
use std::marker::{Send, Sync};
use std::thread;
use core::fmt::Debug;

use crate::mat::Mat;

const MATRIX_SIZE: usize = 64;
const NUM_NONZERO_ENTRIES: usize = 256;
const NUM_THREADS: usize = 8;

fn random_scalar<R: Rng + ?Sized>(rng: &mut R) -> Scalar {
    Scalar::from_raw([rng.gen(), rng.gen(), rng.gen(), rng.gen()])
}

pub fn inner_product<T, U, V>(vec_a: &Vec<T>, vec_b: &Vec<U>) -> V
where
    T: Clone + Copy + Mul<U, Output =V>,
    U: Clone + Copy,
    V: Clone + Copy + Add + Zero + From<u64>,
{
    vec_a.iter()
    .zip(vec_b.iter())
    .fold(V::zero(), |acc, (&a_val, &b_val)| acc + a_val * b_val)
}


fn proj_left<T, U, V>(mat_a: &Mat<T>, v_base: &Vec<U>) -> Vec<V>
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
        let handle = thread::spawn(move || {
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
            result[i] += thread_result[i].clone();
        }
    }

    Arc::try_unwrap(result).unwrap().into_inner().unwrap()

}

fn dirac() {
    let mut rng = rand::thread_rng();
    let m = Arc::new((0..NUM_NONZERO_ENTRIES).map(|_| 
        (rng.gen_range(0..MATRIX_SIZE), rng.gen_range(0..MATRIX_SIZE), random_scalar(&mut rng))
    ).collect::<Vec<_>>());
    
    let v1 = Arc::new((0..MATRIX_SIZE).map(|_| G1Projective::generator() * random_scalar(&mut rng)).collect::<Vec<_>>());

    let aggregate_intermediate = m.chunks(NUM_NONZERO_ENTRIES / NUM_THREADS).map(|chunk| {
        let v1 = Arc::clone(&v1);
        let chunk = chunk.to_owned();
        thread::spawn(move || chunk.iter().fold(vec![G1Projective::identity(); MATRIX_SIZE], |mut acc, &(row, col, ref val)| {
            acc[row] += v1[col] * val;
            acc
        }))
    }).fold(vec![G1Projective::identity(); MATRIX_SIZE], |acc, handle| 
        handle.join().unwrap().iter().enumerate().map(|(i, &l)| l + acc[i]).collect()
    );

    let v2 = (0..MATRIX_SIZE).map(|_| G2Projective::generator() * random_scalar(&mut rng)).collect::<Vec<_>>();
    let result = aggregate_intermediate.iter().enumerate().map(|(i, &aggr)| 
        pairing(&G1Affine::from(aggr), &G2Affine::from(&v2[i]))
    ).fold(Gt::identity(), |acc, x| acc + x);

    println!("Result: {:?}", result);
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_data::*;
    use crate::curve::{ZpElement, G1Element, G2Element, GtElement};

    #[test]
    fn test_dirac() {
        // ga = a left-muliplied by g 
        let mat_a: Mat<ZpElement> = gen_mat_a_direct();
        let vec_g: Vec<G1Element> = gen_vec_v_g1_direct();
        let vec_ga: Vec<G1Element> = gen_vec_va_g1_from_kronecker();
        let vec_h: Vec<G2Element> = gen_vec_v_g2_direct();
        let gah: GtElement = gen_vav_gt_direct();

        let proj_test: Vec<G1Element> = proj_left(&mat_a, &vec_g);
        assert_eq!(proj_test, vec_ga);
        let gah_test = inner_product(&proj_test, &vec_h);

        assert_eq!(gah_test, gah);

        // dirac();
    }
}
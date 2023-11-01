use std::ops::{Add, Mul, AddAssign};

use bls12_381::{pairing, G1Affine, G1Projective, G2Affine, G2Projective, Gt, Scalar};
use rand::Rng;
use std::sync::Arc;
use std::thread;

use crate::mat::Mat;
use crate::curve::{ZpElement, G1Element, G2Element, GtElement};

fn random_scalar<R: Rng + ?Sized>(rng: &mut R) -> Scalar {
    Scalar::from_raw([rng.gen(), rng.gen(), rng.gen(), rng.gen()])
}

const MATRIX_SIZE: usize = 64;
const NUM_NONZERO_ENTRIES: usize = 256;
const NUM_THREADS: usize = 8;

fn proj_left<T, U, V>(mat_a: &Mat<T>, vec_l: &Vec<U>) -> Vec<V>
where
    T: Clone + Mul<U, Output = V>,
    U: Clone,
    V: Clone + Add + AddAssign,
{
    let a_data  = &mat_a.data;
    let result = Vec::new();
    result
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
    use crate::{test_data::*, curve::ZpElement};

    #[test]
    fn test_dirac() {
        // ga = a left-muliplied by g 
        let mat_a: Mat<ZpElement> = gen_mat_a_direct();
        let vec_g: Vec<G1Element> = gen_vec_v_g1_direct();
        let vec_ga: Vec<G1Element> = gen_vec_va_g1_from_kronecker();
        dirac();
    }
}
use bls12_381::*;
use bls12_381::Scalar;
use hex;
use curv::BigInt;
use curv::arithmetic::Zero;
use curv::arithmetic::Converter;
use curv::arithmetic::*;
use curv::arithmetic::traits::*;
use ndarray::Array2;

// Validated
pub fn group_order_bls12_384() -> BigInt {
    // BLS12-381 scalar field $\mathbb{F}_q$
    //! according to https://github.com/zkcrypto/bls12_381/blob/main/src/scalar.rs
    let q = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";
    return BigInt::from_str_radix(q, 16).unwrap();
}


// Reversion is erquired according to https://github.com/zkcrypto/bls12_381/blob/main/src/scalar.rs
pub fn BigInt_to_Scalar(input_integer: &BigInt) -> Scalar {
    let mut hex_string = input_integer.to_str_radix(16);
    if hex_string.len() % 2 != 0 {
        hex_string.insert(0, '0'); // Add a leading zero if the string length is odd.
    }
    let bytes = hex::decode(hex_string).unwrap();
    let mut result = [0u8; 64];
    result[64 - bytes.len()..].copy_from_slice(&bytes);
    result.reverse();
    // println!("Result: {:?}", result);
    return Scalar::from_bytes_wide(&result);
}

pub fn g1_point_add(a: &G1Affine, b: &G1Affine) -> G1Affine{
    let sum_projective = G1Projective::from(a) +  G1Projective::from(b);
    G1Affine::from(sum_projective)
}

pub fn g2_point_add(a: &G2Affine, b: &G2Affine) -> G2Affine{
    let sum_projective = G2Projective::from(a) +  G2Projective::from(b);
    G2Affine::from(sum_projective)
}

pub fn max_of_three(a: usize, b: usize, c: usize) -> usize {
    std::cmp::max(std::cmp::max(a, b), c)
}

pub fn pairing_projective(a: &G1Projective, b: &G2Projective) -> Gt {
    let a_affine = G1Affine::from(*a);
    let b_affine = G2Affine::from(*b);
    return pairing(&a_affine, &b_affine)
}

pub fn inner_product_scalar_scalar(a: &Vec<Scalar>, b: &Vec<Scalar>) -> Scalar {
    let mut result = Scalar::zero();
    for i in 0..a.len() {
        result += a[i] * b[i];
    }
    return result;
}

pub fn inner_product_scalar_G1Projective(a: &Vec<Scalar>, b: &Vec<G1Projective>) -> G1Projective{
    let mut result = G1Projective::identity();
    for i in 0..a.len() {
        result += b[i] * a[i];
    }
    return result;
}


pub fn inner_product_scalar_G2Projective(a: &Vec<Scalar>, b: &Vec<G2Projective>) -> G2Projective{
    let mut result = G2Projective::identity();
    for i in 0..a.len() {
        result += b[i] * a[i];
    }
    return result;
}

pub fn inner_product_G1Projective_G2Projective(a: &Vec<G1Projective>, b: &Vec<G2Projective>) -> Gt{
    let mut result = Gt::identity();
    for i in 0..a.len() {
        result += pairing_projective(&a[i], &b[i]);
    }
    return result;
}

pub fn matmul_scalar(a: &Array2<Scalar>, b: &Array2<Scalar>) -> Array2<Scalar> {
    // Perform matrix multiplication
    let (n1, m1) = a.dim();
    let (n2, m2) = b.dim();
    let mut result_vec = vec![];
    for i in 0..n1 {
        for j in 0..m2 {
            let mut sum= Scalar::zero();
            for k in 0..m1 {
                sum += &a[(i, k)] * &b[(k, j)];
            }
            result_vec.push(sum);
        }
    }

    let result = Array2::from_shape_vec((n1, m2), result_vec).unwrap();
    return result;
}
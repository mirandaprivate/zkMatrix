use bls12_381::*;
use bls12_381::Scalar;
use hex;
use curv::BigInt;
use curv::arithmetic::Zero;
use curv::arithmetic::Converter;
use curv::arithmetic::*;
use curv::arithmetic::traits::*;

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

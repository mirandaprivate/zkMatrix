use curv::BigInt;
use curv::arithmetic::traits::Converter;
use bls12_381::{Scalar, G1Projective, G2Projective, Gt};

pub type ZpElement = Scalar;
pub type IntElement = BigInt;
pub type G1Element = G1Projective;
pub type G2Element = G2Projective;
pub type GtElement = Gt;

// Validated
pub fn group_order() -> BigInt {
    // BLS12-381 scalar field $\mathbb{F}_q$
    /// according to 
    /// <https://github.com/zkcrypto/bls12_381/blob/main/src/scalar.rs>
    let q: &str = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";
    return BigInt::from_str_radix(q, 16).unwrap();
}


// Reversion is erquired according to 
// https://github.com/zkcrypto/bls12_381/blob/main/src/scalar.rs
pub fn int_to_scalar(input_integer: &BigInt) -> Scalar {
    let mut hex_string: String = input_integer.to_str_radix(16);
    if hex_string.len() % 2 != 0 {
        hex_string.insert(0, '0'); // Add a leading zero if the string length is odd.
    }
    let bytes: Vec<u8> = hex::decode(hex_string).unwrap();
    let mut result: [u8; 64] = [0u8; 64];
    result[64 - bytes.len()..].copy_from_slice(&bytes);
    result.reverse();
    // println!("Result: {:?}", result);
    return Scalar::from_bytes_wide(&result);
}

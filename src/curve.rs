//! Wrap the curve used in the library.
//! 
//! We use the BLS12-381 curve, which is a pairing-friendly curve.
//! 
//! We wrap the types and ops here to make the code more readable.
//! 
//! Note that we define the multiplication between a G1 element and a G2 element
//! as the pairing operation.
//!  
use core::convert::From;
use std::ops::{Add, Mul, Neg, Sub, AddAssign};

use bls12_381::{Scalar, G1Projective, G2Projective, Gt, G1Affine, G2Affine};
use curv::BigInt;
use curv::arithmetic::{Samplable, Zero};
use curv::arithmetic::traits::Converter;


#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct ZpElement { pub value: Scalar }

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct G1Element { pub value: G1Projective }

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct G2Element { pub value: G2Projective }

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct GtElement { pub value: Gt }


unsafe impl Send for ZpElement {}
unsafe impl Sync for ZpElement {}

unsafe impl Send for G1Element {}
unsafe impl  Sync for G1Element {}
    
unsafe impl Send for G2Element {}
unsafe impl Sync for G2Element {}

unsafe impl Send for GtElement {}
unsafe impl Sync for GtElement {}

fn bigint_to_scalar(input_integer: &BigInt) -> Scalar {
    let mut hex_string = input_integer.to_str_radix(16);
    if hex_string.len() % 2 != 0 {
        // Add a leading zero if the string length is odd.
        hex_string.insert(0, '0'); 
    }
    let bytes = hex::decode(hex_string).unwrap();
    let mut result = [0u8; 64];
    result[64 - bytes.len()..].copy_from_slice(&bytes);
    result.reverse();
    // println!("Result: {:?}", result);
    return Scalar::from_bytes_wide(&result);
}

impl ZpElement {
    pub fn rand(&self) -> Self {
        let rand_in_256bits: BigInt = BigInt::strict_sample(256);
        let scal: Scalar = bigint_to_scalar(&rand_in_256bits);
        ZpElement { value: scal }
    }

    pub fn pow(&self, exponent: u64) -> Self {
        let exponent_256bits:[u64;4] = [exponent.clone(), 0, 0, 0]; 
        ZpElement { value: self.value.pow(&exponent_256bits)}
    }

}

impl G1Element {
    pub fn generator() -> Self {
        G1Element { value: G1Projective::generator() }
    }
}

impl G2Element {
    pub fn generator() -> Self {
        G2Element { value: G2Projective::generator() }
    }
}

impl GtElement {
    pub fn generator() -> Self {
        let g1_gen = G1Affine::from(G1Projective::generator());
        let g2_gen = G2Affine::from(G2Projective::generator());
        GtElement { value: bls12_381::pairing(&g1_gen, &g2_gen) }
    }
}

impl Zero for ZpElement {
    fn zero() -> Self {
        ZpElement { value: Scalar::zero() }
    }

    fn is_zero(&self) -> bool {
        return self.value == Scalar::zero();
    }
}

impl Zero for G1Element {
    fn zero() -> Self {
        G1Element { value: G1Projective::identity() }
    }

    fn is_zero(&self) -> bool {
        return self.value == G1Projective::identity();
    }
}

impl Zero for G2Element{
    fn zero() -> Self {
        G2Element { value: G2Projective::identity() }
    }

    fn is_zero(&self) -> bool {
        return self.value == G2Projective::identity();
    }
}

impl Zero for GtElement{
    fn zero() -> Self {
        GtElement { value: Gt::identity() }
    }

    fn is_zero(&self) -> bool {
        return self.value == Gt::identity();
    }
}


impl From<BigInt> for ZpElement {
    fn from(input_integer: BigInt) -> Self {
        ZpElement { value: bigint_to_scalar(&input_integer) }
    }
}

impl From<u64> for ZpElement {
    fn from(input_integer: u64) -> Self {
        ZpElement { 
            value: bigint_to_scalar( & BigInt::from(input_integer)) 
        }
    }
}

impl From<u64> for G1Element {
    fn from(input_integer: u64) -> Self {
        G1Element { 
            value: G1Projective::generator() 
            * bigint_to_scalar( & BigInt::from(input_integer))
        }
    }
}

impl From<u64> for G2Element {
    fn from(input_integer: u64) -> Self {
        G2Element { 
            value: G2Projective::generator() 
            * bigint_to_scalar( & BigInt::from(input_integer))
        }
    }
}

impl From<u64> for GtElement {
    fn from(input_integer: u64) -> Self {
        let g1_value = G1Affine::from(
            G1Projective::generator()
            * bigint_to_scalar( & BigInt::from(input_integer))
        );
        let g2_gen = G2Affine::from(G2Projective::generator());
        GtElement { 
            value: bls12_381::pairing(&g1_value, &g2_gen) 
        }
    }
}

impl Add for ZpElement {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        ZpElement { value: self.value + &rhs.value }
    }
}

impl Add for G1Element {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        G1Element { value: self.value + &rhs.value }
    }
}

impl Add for G2Element {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        G2Element { value: self.value + &rhs.value }
    }
}

impl Add for GtElement {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        GtElement { value: self.value + &rhs.value }
    }
}

impl AddAssign for ZpElement {
    fn add_assign(&mut self, rhs: Self) {
        self.value += &rhs.value;
    }
}

impl AddAssign for G1Element {
    fn add_assign(&mut self, rhs: Self) {
        self.value += &rhs.value;
    }
}

impl AddAssign for G2Element {
    fn add_assign(&mut self, rhs: Self) {
        self.value += &rhs.value;
    }
}

impl AddAssign for GtElement {
    fn add_assign(&mut self, rhs: Self) {
        self.value += &rhs.value;
    }
}

impl Sub for ZpElement {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        ZpElement { value: self.value - &rhs.value }
    }
}

impl Neg for ZpElement {
    type Output = Self;

    fn neg(self) -> Self::Output {
        ZpElement { value: -self.value }
    }
}

impl Sub for G1Element {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        G1Element { value: self.value - &rhs.value }
    }
}

impl Neg for G1Element {
    type Output = Self;

    fn neg(self) -> Self::Output {
        G1Element { value: -self.value }
    }
}

impl Sub for G2Element {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        G2Element { value: self.value - &rhs.value }
    }
}

impl Neg for G2Element {
    type Output = Self;

    fn neg(self) -> Self::Output {
        G2Element { value: -self.value }
    }
}

impl Sub for GtElement{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        GtElement { value: self.value - &rhs.value }
    }
}

impl Neg for GtElement {
    type Output = Self;

    fn neg(self) -> Self::Output {
        GtElement { value: -self.value }
    }
}

impl Mul for ZpElement {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        ZpElement { value: self.value * &rhs.value }
    }
}

impl Mul<G1Element> for ZpElement {
    type Output = G1Element;

    fn mul(self, rhs: G1Element) -> Self::Output {
        G1Element { value: self.value * &rhs.value }
    }
}

impl Mul<G2Element> for ZpElement {
    type Output = G2Element;

    fn mul(self, rhs: G2Element) -> Self::Output {
        G2Element { value: self.value * &rhs.value}
    }
}

impl Mul<GtElement> for ZpElement {
    type Output = GtElement;

    fn mul(self, rhs: GtElement) -> Self::Output {
        GtElement { value: &rhs.value * self.value}
    }
}

impl Mul<ZpElement> for G1Element {
    type Output = G1Element;

    fn mul(self, rhs: ZpElement) -> Self::Output {
        G1Element { value: &rhs.value * self.value }
    }
}

impl Mul<ZpElement> for G2Element {
    type Output = G2Element;

    fn mul(self, rhs: ZpElement) -> Self::Output {
        G2Element { value: &rhs.value * self.value }
    }
}

impl Mul<ZpElement> for GtElement {
    type Output = GtElement;

    fn mul(self, rhs: ZpElement) -> Self::Output {
        GtElement { value: &self.value * rhs.value }
    }
}

impl Mul<G2Element> for G1Element {
    type Output = GtElement;

    fn mul(self, rhs: G2Element) -> Self::Output {
        let lfs_affine = G1Affine::from(self.value);
        let rhs_affine = G2Affine::from(rhs.value);
        GtElement { value: bls12_381::pairing(&lfs_affine, &rhs_affine) }
    }
}

impl Mul<G1Element> for G2Element {
    type Output = GtElement;

    fn mul(self, rhs: G1Element) -> Self::Output {
        let lfs_affine = G2Affine::from(self.value);
        let rhs_affine = G1Affine::from(rhs.value);
        GtElement { value: bls12_381::pairing(&rhs_affine, &lfs_affine) }
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_curve() {
        let scalar_zero = ZpElement::zero();
        let g1_zero = G1Element::zero();
        let g2_zero = G2Element::zero();
        let gt_zero = GtElement::zero();

        assert_eq!(scalar_zero, ZpElement::from(0));
        assert_eq!(g1_zero, G1Element::from(0));
        assert_eq!(g2_zero, G2Element::from(0));
        assert_eq!(gt_zero, GtElement::from(0));

        let g1_value = G1Element::from(2);
        let g2_value = G2Element::from(3);
        let gt_value = GtElement::from(6);

        assert_eq!(g1_value * g2_value, gt_value);

        let mut g1_mut = G1Element::from(1);
        let mut g2_mut = G2Element::from(2);
        let mut gt_mut = GtElement::from(5);

        g1_mut += G1Element::generator();
        g2_mut += G2Element::generator();
        gt_mut += GtElement::generator();

        assert_eq!(g1_mut, g1_value);
        assert_eq!(g2_mut, g2_value);
        assert_eq!(gt_mut, gt_value);

        let zp_value = ZpElement::from(2);
        let zp_exp = zp_value.pow(3);

        assert_eq!(zp_exp, ZpElement::from(8));

        assert_eq!(scalar_zero + scalar_zero, ZpElement::zero());
        assert_eq!(scalar_zero - scalar_zero, ZpElement::zero());
        assert_eq!(- scalar_zero - scalar_zero, ZpElement::zero());
        assert_eq!(- scalar_zero * scalar_zero, ZpElement::zero());

        assert_eq!(g1_value + g1_zero, g1_value);
        assert_eq!(g1_value - g1_zero, g1_value);
        assert_eq!(- g1_value - g1_zero, - g1_value);
        assert_eq!(- g1_value * scalar_zero, g1_zero);
        assert_eq!(scalar_zero * g1_value, g1_zero);

        assert_eq!(g2_value + g2_zero, g2_value);
        assert_eq!(g2_value - g2_zero, g2_value);
        assert_eq!(- g2_value - g2_zero, - g2_value);
        assert_eq!(- g2_value * scalar_zero, g2_zero);
        assert_eq!(scalar_zero * g2_value, g2_zero);

        assert_eq!(gt_value + gt_zero, gt_value);
        assert_eq!(gt_value - gt_zero, gt_value);
        assert_eq!(- gt_value - gt_zero, - gt_value);
        assert_eq!(- gt_value * scalar_zero, gt_zero);
        assert_eq!(scalar_zero * gt_value, gt_zero);
        assert_eq!(g1_zero * g2_zero, gt_zero);

    }

}
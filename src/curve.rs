//! Defines the curve used in the library.
//! We use the BLS12-381 curve, which is a pairing-friendly curve.
//! 
//! We wrap the types and ops here to make the code more readable.
//! 
//!  
use core::convert::From;
use std::ops::{Add, Mul, Neg, Sub, AddAssign};

use bls12_381::{Scalar, G1Projective, G2Projective, Gt, G1Affine, G2Affine};
use bls12_381::pairing;
use curv::BigInt;
use curv::arithmetic::Zero;
use curv::arithmetic::traits::Converter;



#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct ZpElement { value: Scalar }

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct G1Element { value: G1Projective }

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct G2Element { value: G2Projective }

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct GtElement { value: Gt }

impl ZpElement{
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

impl Zero for ZpElement {
    fn zero() -> Self {
        ZpElement { value: Scalar::zero() }
    }

    fn is_zero(&self) -> bool {
        return self.value == Scalar::zero();
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

impl Add for ZpElement {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        ZpElement { value: self.value + &rhs.value }
    }
}

impl Add for G1Element{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        G1Element { value: self.value + &rhs.value }
    }
}

impl Add for G2Element{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        G2Element { value: self.value + &rhs.value }
    }
}

impl Add for GtElement{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        GtElement { value: self.value + &rhs.value }
    }
}

impl AddAssign for ZpElement{
    fn add_assign(&mut self, rhs: Self) {
        self.value += &rhs.value;
    }
}

impl AddAssign for G1Element{
    fn add_assign(&mut self, rhs: Self) {
        self.value += &rhs.value;
    }
}

impl AddAssign for G2Element{
    fn add_assign(&mut self, rhs: Self) {
        self.value += &rhs.value;
    }
}

impl AddAssign for GtElement{
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



impl Mul for ZpElement{
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

impl Mul<ZpElement> for G1Element{
    type Output = G1Element;

    fn mul(self, rhs: ZpElement) -> Self::Output {
        G1Element { value: &rhs.value * self.value }
    }
}

impl Mul<ZpElement> for G2Element{
    type Output = G2Element;

    fn mul(self, rhs: ZpElement) -> Self::Output {
        G2Element { value: &rhs.value * self.value }
    }
}

impl Mul<ZpElement> for GtElement{
    type Output = GtElement;

    fn mul(self, rhs: ZpElement) -> Self::Output {
        GtElement { value: &self.value * rhs.value }
    }
}

impl Mul<G2Element> for G1Element{
    type Output = GtElement;

    fn mul(self, rhs: G2Element) -> Self::Output {
        let lfs_affine = G1Affine::from(self.value);
        let rhs_affine = G2Affine::from(rhs.value);
        GtElement { value: pairing(&lfs_affine, &rhs_affine) }
    }
}

impl Mul<G1Element> for G2Element{
    type Output = GtElement;

    fn mul(self, rhs: G1Element) -> Self::Output {
        let lfs_affine = G2Affine::from(self.value);
        let rhs_affine = G1Affine::from(rhs.value);
        GtElement { value: pairing(&rhs_affine, &lfs_affine) }
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zeros() {
        let scalar = ZpElement::zero();
        let mut a = ZpElement::from(0);
        // let b = IntElement::from(0);
        // let c = G1Element::identity();
        // let d = G2Element::identity();
        // let e = GtElement::identity();
        a += scalar;
        let b = a.pow(3);
        assert_eq!(b, ZpElement::zero());
        assert_eq!(a + scalar, ZpElement::zero());
        assert_eq!(scalar - scalar, ZpElement::zero());
        assert_eq!(- scalar - scalar, ZpElement::zero());
        assert_eq!(- scalar * scalar, ZpElement::zero());

        let g1 = G1Element { value: G1Projective::identity() };
        let g2 = G2Element { value: G2Projective::identity() };
        let gt = GtElement::zero();
        
        assert_eq!( g1 * g2, gt);
        assert_eq!( a * g1 , g1);
        assert_eq!( a * g2 , g2);
        assert_eq!( a * gt , gt);
        assert_eq!( g2 * g1, gt);
        assert_eq!( g1 * a , g1);
        assert_eq!( g2 * a , g2);
        assert_eq!( gt * a , gt);
        // assert_eq!(b, BigInt::from(0));
        // assert_eq!(c, G1Projective::identity());
        // assert_eq!(d, G2Projective::identity());
        // assert_eq!(e, Gt::identity());
    }

    // #[test]
    // fn test_arithmetic(){
    //     let a = ZpElement::from(2);
    //     let b = ZpElement::from(3);
    //     let c = a + b;
    //     let d = a * b;
    //     let e = a - b;
    //     let f = -a;
    //     assert_eq!(c, Scalar::from(5));
    //     assert_eq!(d, Scalar::from(6));
    //     assert_eq!(e, - Scalar::from(1));
    //     assert_eq!(f, - Scalar::from(2));
    // }
}
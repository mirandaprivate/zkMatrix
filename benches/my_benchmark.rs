#![allow(warnings)]
#![feature(test)]

extern crate test;

use bls12_381::*;
use bls12_381::Scalar;
use hex;
use curv::BigInt;
use curv::arithmetic::Zero;
use curv::arithmetic::Converter;
use curv::arithmetic::*;
use zkmatrix::utils::*;
use test::Bencher;

#[bench]
fn bench_scalar_mul_256bits(b: &mut Bencher){
    let l_in_256bits: BigInt = BigInt::strict_sample(256);
    let r_in_256bits: BigInt = BigInt::strict_sample(256);
    let l_scalar: Scalar = BigInt_to_Scalar(&l_in_256bits);
    let r_scalar: Scalar = BigInt_to_Scalar(&r_in_256bits);
    b.iter(|| l_scalar * r_scalar );
}

#[bench]
fn bench_scalar_add_256bits(b: &mut Bencher){
    let l_in_256bits: BigInt = BigInt::strict_sample(256);
    let r_in_256bits: BigInt = BigInt::strict_sample(256);
    let l_scalar: Scalar = BigInt_to_Scalar(&l_in_256bits);
    let r_scalar: Scalar = BigInt_to_Scalar(&r_in_256bits);
    b.iter(|| l_scalar + r_scalar );
}

#[bench]
fn bench_G1_exp_256bits(b: &mut Bencher) {
    let rand_in_256bits: BigInt = BigInt::strict_sample(256);
    let scal: Scalar = BigInt_to_Scalar(&rand_in_256bits);
    b.iter(|| G1Affine::generator() * scal);
}


#[bench]
fn bench_G1_add_256bits_projective(b: &mut Bencher) {
    let rand_in_256bits: BigInt = BigInt::strict_sample(256);
    let scal: Scalar = BigInt_to_Scalar(&rand_in_256bits);
    let base: G1Projective = G1Projective::generator() * scal;
    b.iter(|| base + base);
}


#[bench]
fn bench_G2_exp_256bits(b: &mut Bencher) {
    let rand_in_256bits: BigInt = BigInt::strict_sample(256);
    let scal: Scalar = BigInt_to_Scalar(&rand_in_256bits);
    b.iter(|| G2Affine::generator() * scal);
}

#[bench]
fn bench_G2_exp_256bits_projective(b: &mut Bencher) {
    let rand_in_256bits: BigInt = BigInt::strict_sample(256);
    let scal: Scalar = BigInt_to_Scalar(&rand_in_256bits);
    let base: G2Projective = G2Projective::generator() * scal;
    b.iter(|| G2Projective::generator() * scal);
}

#[bench]
fn bench_G2_add_256bits_projective(b: &mut Bencher) {
    let rand_in_256bits: BigInt = BigInt::strict_sample(256);
    let scal: Scalar = BigInt_to_Scalar(&rand_in_256bits);
    let base: G2Projective = G2Projective::generator() * scal;
    b.iter(|| base + base);
}


#[bench]
fn bench_GT_exp_256bits(b: &mut Bencher) {
    let x: Scalar = BigInt_to_Scalar(&BigInt::from(12345));
    let y: Scalar = BigInt_to_Scalar(&BigInt::from(12345));
    let g: G1Affine = G1Affine::from(G1Affine::generator() * x);
    let h: G2Affine = G2Affine::from(G2Affine::generator() * y);
    let p: Gt = pairing(&g, &h);
    let rand_in_256bits: BigInt = BigInt::strict_sample(256);
    let scal: Scalar = BigInt_to_Scalar(&rand_in_256bits);
    b.iter(|| p * scal);
}


#[bench]
fn bench_pairing(b: &mut Bencher) {
    let rand1: BigInt = BigInt::strict_sample(256);
    let scal1: Scalar = BigInt_to_Scalar(&rand1);
    let rand2: BigInt = BigInt::strict_sample(256);
    let scal2: Scalar = BigInt_to_Scalar(&rand2);
    let g: G1Affine = G1Affine::from(G1Affine::generator() * scal1);
    let h: G2Affine = G2Affine::from(G2Affine::generator() * scal2);
    b.iter(|| pairing(&g, &h));
}

#[bench]
fn bench_Pairing_Projective(b: &mut Bencher) {
    let rand1: BigInt = BigInt::strict_sample(256);
    let scal1: Scalar = BigInt_to_Scalar(&rand1);
    let rand2: BigInt = BigInt::strict_sample(256);
    let scal2: Scalar = BigInt_to_Scalar(&rand2);
    let g: G1Projective = G1Projective::generator() * scal1;
    let h: G2Projective = G2Projective::generator() * scal2;
    b.iter(|| pairing_projective(&g, &h));
}
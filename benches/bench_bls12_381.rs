#![allow(warnings)]
#![feature(test)]

extern crate test;
extern crate rand;

use std::arch::x86_64::_MM_FROUND_RAISE_EXC;

use rand::Rng;

use zkmatrix::utils::curve::{ZpElement, G1Element, G2Element, GtElement, Double};
use zkmatrix::mat::Mat;
use zkmatrix::commit_mat::CommitMat;
use zkmatrix::test_data::*;

use test::Bencher;


#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_scalar_add_scalar(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_scalar =  ZpElement::from(l_u64);
    let r_scalar =  ZpElement::from(r_u64);
    b.iter(|| l_scalar + r_scalar );
}

#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_g1_add_g1(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_g1 =  G1Element::from(l_u64);
    let r_g1 =  G1Element::from(r_u64);
    b.iter(|| l_g1 + r_g1 );
}

#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_g2_add_g2(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_g2 =  G2Element::from(l_u64);
    let r_g2 =  G2Element::from(r_u64);
    b.iter(|| l_g2 + r_g2 );
}


#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_gt_add_gt(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_gt =  GtElement::from(l_u64);
    let r_gt =  GtElement::from(r_u64);
    b.iter(|| l_gt + r_gt );
}

#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_g1_mul_scalar(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_g1 =  G1Element::from(l_u64);
    let r_scalar =  ZpElement::from(r_u64);
    b.iter(|| l_g1 * r_scalar );
}


#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_g2_mul_scalar(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_g2 =  G2Element::from(l_u64);
    let r_scalar =  ZpElement::from(r_u64);
    b.iter(|| l_g2 * r_scalar );
}

#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_gt_mul_scalar(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_gt =  GtElement::from(l_u64);
    let r_scalar =  ZpElement::from(r_u64);
    b.iter(|| l_gt * r_scalar );
}

#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_g1_mul_g2(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_g1 =  G1Element::from(l_u64);
    let r_g2 =  G2Element::from(r_u64);
    b.iter(|| l_g1 * r_g2 );
}


#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_scalar_mul_scalar(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_scalar =  ZpElement::from(l_u64);
    let r_scalar =  ZpElement::from(r_u64);
    b.iter(|| l_scalar * r_scalar );
}

#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_scalar_mac_scalar(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let a_u64: u64 = rand::thread_rng().gen();
    let l_scalar =  ZpElement::from(l_u64);
    let r_scalar =  ZpElement::from(r_u64);
    let a_scalar =  ZpElement::from(a_u64);
    b.iter(|| a_scalar + l_scalar * r_scalar );
}



#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_u64_mul_u64(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    b.iter(|| l_u64 * r_u64 );
}

#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_u64_add_u64(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    b.iter(|| l_u64 + r_u64 );
}

#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_scalar_double(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let l_scalar =  ZpElement::from(l_u64);
    b.iter(|| l_scalar.double() );
}

#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_g1_double(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let l_g1 =  G1Element::from(l_u64);
    b.iter(|| l_g1.double() );
}

#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_g2_double(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let l_g2 =  G2Element::from(l_u64);
    b.iter(|| l_g2.double() );
}


#[bench]
#[ignore = "only for testing bls12_381"]
fn bench_gt_double(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let l_gt =  GtElement::from(l_u64);
    b.iter(|| l_gt.double() );
}


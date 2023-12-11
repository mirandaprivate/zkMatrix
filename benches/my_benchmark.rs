#![allow(warnings)]
#![feature(test)]

extern crate test;
extern crate rand;

use std::arch::x86_64::_MM_FROUND_RAISE_EXC;

use rand::Rng;

use zkmatrix::curve::{ZpElement, G1Element, G2Element, GtElement, Double};
use zkmatrix::mat::Mat;
use zkmatrix::commit_mat::CommitMat;
use zkmatrix::dirac;
use zkmatrix::test_data::*;

use test::Bencher;


#[bench]
fn bench_scalar_add_scalar(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_scalar =  ZpElement::from(l_u64);
    let r_scalar =  ZpElement::from(r_u64);
    b.iter(|| l_scalar + r_scalar );
}

#[bench]
fn bench_g1_add_g1(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_g1 =  G1Element::from(l_u64);
    let r_g1 =  G1Element::from(r_u64);
    b.iter(|| l_g1 + r_g1 );
}

#[bench]
fn bench_g2_add_g2(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_g2 =  G1Element::from(l_u64);
    let r_g2 =  G1Element::from(r_u64);
    b.iter(|| l_g2 + r_g2 );
}


#[bench]
fn bench_gt_add_gt(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_gt =  G1Element::from(l_u64);
    let r_gt =  G1Element::from(r_u64);
    b.iter(|| l_gt + r_gt );
}

#[bench]
fn bench_g1_mul_scalar(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_g1 =  G1Element::from(l_u64);
    let r_scalar =  ZpElement::from(r_u64);
    b.iter(|| l_g1 * r_scalar );
}


#[bench]
fn bench_g2_mul_scalar(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_g2 =  G2Element::from(l_u64);
    let r_scalar =  ZpElement::from(r_u64);
    b.iter(|| l_g2 * r_scalar );
}

#[bench]
fn bench_gt_mul_scalar(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_gt =  GtElement::from(l_u64);
    let r_scalar =  ZpElement::from(r_u64);
    b.iter(|| l_gt * r_scalar );
}

#[bench]
fn bench_g1_mul_g2(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_g1 =  G1Element::from(l_u64);
    let r_g2 =  G2Element::from(r_u64);
    b.iter(|| l_g1 * r_g2 );
}


#[bench]
fn bench_scalar_mul_scalar(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_scalar =  ZpElement::from(l_u64);
    let r_scalar =  ZpElement::from(r_u64);
    b.iter(|| l_scalar * r_scalar );
}

#[bench]
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
fn bench_u64_mul_u64(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    b.iter(|| l_u64 * r_u64 );
}

#[bench]
fn bench_u64_add_u64(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    b.iter(|| l_u64 + r_u64 );
}

#[bench]
fn bench_scalar_double(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let l_scalar =  ZpElement::from(l_u64);
    b.iter(|| l_scalar.double() );
}

#[bench]
fn bench_g1_double(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let l_g1 =  G1Element::from(l_u64);
    b.iter(|| l_g1.double() );
}

#[bench]
fn bench_g2_double(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let l_g2 =  G2Element::from(l_u64);
    b.iter(|| l_g2.double() );
}


#[bench]
fn bench_gt_double(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let l_gt =  GtElement::from(l_u64);
    b.iter(|| l_gt.double() );
}


#[bench]
#[ignore = "expensive"]
fn bench_dirac(b: &mut Bencher){
    let mat_a: Mat<u64> = gen_mat_a_u64_direct();
    let vec_g: Vec<G1Element> = gen_vec_v_g1_direct();
    let vec_h: Vec<G2Element> = gen_vec_v_g2_direct();
    b.iter(|| dirac::dirac(&vec_g, &mat_a, &vec_h) );
}

#[bench]
fn bench_mat_reduce_scalar(b: &mut Bencher){
    let mat_a: Mat<u64> = gen_mat_a_u64_direct();
    let left: Vec<ZpElement> = gen_vec_v_direct();
    let right: Vec<ZpElement> = gen_vec_v_direct();
    b.iter(||dirac::inner_product(
        &dirac::proj_left(&mat_a, &left),
        &right) 
    );
}
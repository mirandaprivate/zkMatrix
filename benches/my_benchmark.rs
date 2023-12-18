#![allow(warnings)]
#![feature(test)]

extern crate test;
extern crate rand;

use std::arch::x86_64::_MM_FROUND_RAISE_EXC;

use rand::Rng;

use zkmatrix::utils::curve::{ZpElement, G1Element, G2Element, GtElement, Double};
use zkmatrix::mat::Mat;
use zkmatrix::commit_mat::CommitMat;
use zkmatrix::utils::dirac;
use zkmatrix::test_data::*;

use test::Bencher;

#[bench]
// #[ignore = "expensive"]
fn bench_dirac(b: &mut Bencher){
    let mat_a: Mat<u64> = gen_mat_a_u64_direct_test();
    let vec_g: Vec<G1Element> = gen_vec_v_g1_direct_test();
    let vec_h: Vec<G2Element> = gen_vec_v_g2_direct_test();
    b.iter(|| dirac::dirac(&vec_g, &mat_a, &vec_h) );
}

#[bench]
fn bench_mat_reduce_scalar(b: &mut Bencher){
    let mat_a: Mat<u64> = gen_mat_a_u64_direct_test();
    let left: Vec<ZpElement> = gen_vec_v_direct_test();
    let right: Vec<ZpElement> = gen_vec_v_direct_test();
    b.iter(||dirac::inner_product(
        &dirac::proj_left(&mat_a, &left),
        &right) 
    );
}
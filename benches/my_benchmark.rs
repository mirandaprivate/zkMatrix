#![allow(warnings)]
#![feature(test)]

extern crate test;
extern crate rand;

use std::arch::x86_64::_MM_FROUND_RAISE_EXC;

use rand::Rng;

use zkmatrix::utils::curve::{ZpElement, G1Element, G2Element, GtElement, Double};
use zkmatrix::mat::Mat;
use zkmatrix::commit_mat::CommitMat;
use zkmatrix::utils::dirac::{self, BraKet};
use zkmatrix::utils::test_data::*;

use zkmatrix::experiment_data;
use zkmatrix::setup::SRS;

use test::Bencher;

#[bench]
// #[ignore = "expensive"]
fn bench_dirac(b: &mut Bencher){
    let (_,a, _) = experiment_data::gen_matrices_dense(1024);
    let srs = SRS::new(2048);
    b.iter(|| 
        a.braket(&srs.g_hat_vec, &srs.h_hat_vec) 
    );
}


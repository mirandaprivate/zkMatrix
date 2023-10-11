#![allow(warnings)]
#[cfg(not(nightly))]

use bls12_381::*;
use bls12_381::Scalar;
use hex;
use curv::BigInt;
use curv::arithmetic::Zero;
use curv::arithmetic::Converter;
use curv::arithmetic::traits::*;
use curv::arithmetic::traits::Samplable;
use curv::arithmetic::*;
use crate::utils::*;
use std::time::Instant;
use std::os::raw::c_int;

pub struct HDSIP {
    pub tr: Vec<Gt>,
    pub coins: Vec<Scalar>,
    pub m: usize,
    pub n: usize,
    pub C_c: Gt,
    pub C_a: Gt,
    pub b_vec: Vec<BigInt>,
    pub H_hat: G2Projective,
    pub H_hat_vec: Vec<G2Projective>,
    pub G_hat_vec: Vec<G1Projective>,
    a_vec_vec: Vec<Vec<BigInt>>, // n rows * m columns
    A_vec: Vec<G1Projective>,
}

impl HDSIP{

    // Generate the transcipt
    fn prove(&self) {

        let log_m = (self.m as f64).log2() as usize;
        let log_n = (self.n as f64).log2() as usize;
        let log_mn = ((self.m * self.n) as f64).log2() as usize;
        print!("m = {}, n = {}, ", self.m, self.n);
        println!("log_m = {}, log_n = {}, log_mn = {}", log_m, log_n, log_mn);
        println!("G_hat_vec.len():{}", self.G_hat_vec.len());
        println!("H_hat_vec.len():{}", self.H_hat_vec.len());
        println!("A_vec.len():{}", self.A_vec.len());

    }

    // Verify the transcipt
    fn verify(&self) {

    }

}



pub fn init_HDSIP_inputs(m: usize, n: usize) -> HDSIP {

    let tr: Vec<Gt> = Vec::new();

    let coins: Vec<Scalar> = Vec::new();

    let q = group_order_bls12_384();
    assert!(&BigInt::from(m as u32) <= &q && &BigInt::from(n as u32) <= &q);

    let s_hat = BigInt::sample_below(&q);
    let mut G_hat_vec: Vec<G1Projective> = Vec::new();
    for i in 0..=(m-1) {
        let G_scal: BigInt = BigInt::mod_pow(&s_hat, &BigInt::from(i as u32), &q);
        let G_item: G1Projective = G1Projective::generator() * BigInt_to_Scalar(&G_scal);
        G_hat_vec.push(G_item);
    }


    let mut H_hat_vec: Vec<G2Projective> = Vec::new();
    for i in 0..=(n-1) {
        let iq: BigInt = &BigInt::from(i as u32) * &q;
        let H_scal: BigInt = BigInt::mod_pow(&s_hat, &iq, &q);
        let H_item: G2Projective = G2Projective::generator() * BigInt_to_Scalar(&H_scal);
        H_hat_vec.push(H_item);
    }

    let mut a_vec_vec: Vec<Vec<BigInt>> = Vec::new();
    let mut A_vec: Vec<G1Projective> = Vec::new();
    let mut b_vec: Vec<BigInt> = Vec::new();
    for i in 0..n {
        let mut a_vec: Vec<BigInt> = Vec::new();
        let mut A_i: G1Projective = G1Projective::identity();
        for j in 0..m {
            let a_ij: BigInt = BigInt::sample_below(&q);
            A_i += G_hat_vec[j] * BigInt_to_Scalar(&a_ij);
            a_vec.push(a_ij);
        }
        a_vec_vec.push(a_vec);
        A_vec.push(A_i);
        b_vec.push(BigInt::sample_below(&q));

    } 

    let mut C_a: Gt = Gt::identity();
    for i in 0..n {
        let G1: G1Affine = G1Affine::from(&A_vec[i]);
        let G2: G2Affine = G2Affine::from(&H_hat_vec[i]);
        let pairing_Ai_Hi = pairing(&G1, &G2);
        C_a += pairing_Ai_Hi;
    }

    let a1_bn = BigInt::from(12345);
    let a1 = BigInt_to_Scalar(&a1_bn);
    let a2_bn = BigInt::from(45678);
    let a2: Scalar = BigInt_to_Scalar(&a2_bn);
    let G1: G1Affine = G1Affine::from(G1Affine::generator() * a1);
    let H1: G2Affine = G2Affine::from(G2Affine::generator() * a2);
    let C = pairing(&G1, &H1);

    println!("System paramaters for HDSIP are ready now.");

    HDSIP {
        tr,
        coins,
        m,
        n,
        C_c: C,
        C_a,
        b_vec,
        H_hat: G2Projective::generator(),
        H_hat_vec,
        G_hat_vec,
        a_vec_vec,
        A_vec, 
    }
    
}

#[test]
pub fn test_HDSIP() {
    println!("**************** Experiment for HDSIP ****************");
    let m: usize = 16;
    let n: usize = 8;

    let start = Instant::now();
    let HDSIP_inputs = init_HDSIP_inputs(m, n);
    println!("Time elapsed for init_HDSIP_inputs() is: {:?}", start.elapsed());

    let start = Instant::now();
    HDSIP_inputs.prove();
    println!("Time elapsed for prove() is: {:?}", start.elapsed());

    println!("Hi, HDSIP\n");

}

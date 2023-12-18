#![allow(warnings)]
#![feature(test)]

extern crate curve25519_dalek;
extern crate rand;
extern crate test;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::ristretto::{RistrettoPoint, VartimeRistrettoPrecomputation};
use curve25519_dalek::traits::VartimePrecomputedMultiscalarMul;

use rand::Rng;

use test::Bencher;

#[bench]
#[ignore = "only for testing curve25519"]
fn bench_fixed_base_me(b: &mut Bencher){
    let length = 1024;
    let scalars = vec![
        Scalar::from(rand::thread_rng().gen::<u64>()); length];
    let points = vec![
        RistrettoPoint::mul_base(
            &Scalar::from(rand::thread_rng().gen::<u64>())); length];

    let pre_compute = VartimeRistrettoPrecomputation::new(points);  
    
    b.iter(|| pre_compute.vartime_multiscalar_mul(&scalars));
}

#[bench]
#[ignore = "only for testing curve25519"]
fn bench_fixed_base_precompute(b: &mut Bencher){
    let length = 1024;
    let scalars = vec![
        Scalar::from(rand::thread_rng().gen::<u64>()); length];
    let points = vec![
        RistrettoPoint::mul_base(
            &Scalar::from(rand::thread_rng().gen::<u64>())); length];

    b.iter(|| 
        {   let points_clone = points.clone();
            let pre_compute_clone = VartimeRistrettoPrecomputation::new(points_clone);
        }
    );
}

#[bench]
#[ignore = "only for testing curve25519"]
fn bench_scalar_mul(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_scalar =  Scalar::from(l_u64);
    let r_scalar =  Scalar::from(r_u64);
    b.iter(|| l_scalar * r_scalar );
}

#[bench]
#[ignore = "only for testing curve25519"]
fn bench_g_mul_scalar(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_g1 =  RistrettoPoint::mul_base(&Scalar::from(l_u64));
    let r_scalar =  Scalar::from(r_u64);
    b.iter(|| l_g1 * r_scalar );
}


#[bench]
#[ignore = "only for testing curve25519"]
fn bench_g_add_scalar(b: &mut Bencher){
    let l_u64: u64 = rand::thread_rng().gen();
    let r_u64: u64 = rand::thread_rng().gen();
    let l_g1 =  RistrettoPoint::mul_base(&Scalar::from(l_u64));
    let r_g1=  RistrettoPoint::mul_base(&Scalar::from(r_u64));
    b.iter(|| l_g1 + r_g1 );
}
// #![allow(warnings)]
// #![feature(test)]

// extern crate test;
// extern crate rand;

// use test::Bencher;

// use ark_ec::{pairing::Pairing, AffineRepr};
// use ark_ff::{PrimeField, Field};
// // We'll use the BLS12-381 pairing-friendly group for this example.
// use ark_bls12_381::{
//     Bls12_381, G1Projective as G1, 
//     G2Projective as G2, 
//     G1Affine, G2Affine, 
//     Fr as ScalarField, };
// use ark_std::{Zero, UniformRand};


// #[bench]
// fn bench_pairing(b: &mut Bencher){
//     let mut rng = ark_std::rand::thread_rng();
//     // Let's sample uniformly random field elements:
//     let l: G1Affine = G1::rand(&mut rng).into();
//     let r: G2Affine = G2::rand(&mut rng).into();
//     // We can compute the pairing of `a` and `b`:


//     b.iter(|| Bls12_381::pairing(l, r) );
// }
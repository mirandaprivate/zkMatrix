#![allow(warnings)]

use bls12_381::*;
use bls12_381::Scalar;
use hex;
use curv::BigInt;
use curv::arithmetic::Zero;
use curv::arithmetic::Converter;
use curv::arithmetic::*;
use crate::utils::*;

#[test]
pub fn pairing_demo(){
    println!("Pairing demo is started!");

    /***************   Test 1 Exponentiation + Equivalence Verify  **************/
    let a1_bn = BigInt::from(12345);
    let a1 = BigInt_to_Scalar(&a1_bn);

    assert_eq!(a1, Scalar::from(12345));

    let a2_bn = BigInt::from(45678);
    let a2 = BigInt_to_Scalar(&a2_bn);
    let c = a1 * a2;

    let g = G1Affine::from(G1Affine::generator() * a1);
    let h = G2Affine::from(G2Affine::generator() * a2);
    let p = pairing(&g, &h);

    assert!(p != Gt::identity());

    let expected_1 = G1Affine::from(G1Affine::generator() * c);
    assert_eq!(p, pairing(&expected_1, &G2Affine::generator()));

    let expected_2 = G2Affine::from(G2Affine::generator() * c);
    assert_eq!(p, pairing(&G1Affine::generator(), &expected_2));

    assert_eq!(
        p,
        pairing(&G1Affine::generator(), &G2Affine::generator()) * c
    );


    /***************   Test 2-1 Point addition in G1 **************/
    let g = G1Projective::generator() * Scalar::from(12345);

    // Convert the G1Projective point to G1Affine
    let g_affine = G1Affine::from(g);

    // Convert G1Affine points to G1Projective
    let g_projective_1 = G1Projective::from(g_affine);
    let g_projective_2 = G1Projective::from(g_affine);

    // Add the G1Projective points
    let sum_projective = g_projective_1 + g_projective_2;

    // Convert the result back to G1Affine if needed
    let sum_affine = G1Affine::from(sum_projective);
    assert_eq!(G1Projective::identity() + sum_projective, sum_projective);


    /***************   Test 2-2 Point addition in G2 **************/
    let g = G2Projective::generator() * Scalar::from(12345);

    // Convert the G1Projective point to G1Affine
    let g_affine = G2Affine::from(g);

    // Convert G1Affine points to G1Projective
    let g_projective_1 = G2Projective::from(g_affine);
    let g_projective_2 = G2Projective::from(g_affine);

    // Add the G1Projective points
    let sum_projective = g_projective_1 + g_projective_2;

    // Convert the result back to G1Affine if needed
    let sum_affine = G2Affine::from(sum_projective);
    assert_eq!(G2Projective::identity() + sum_projective, sum_projective);

    // println!("Point 1: {:?}", g_affine);
    // println!("Point 2: {:?}", g_affine);
    // println!("Sum: {:?}", sum_affine);

    /***************   Test 3 Validate in GT **************/
    // reminder: p = pairing(a1, a2)
    assert_eq!(
        (p * Scalar::from(123)) * Scalar::from(456),
        (p * Scalar::from(456)) * Scalar::from(123),
    );
    let p2 = &p + &p; // point add of Gt
    assert_eq!(Gt::identity() + p, p);
    println!("Pairing demo is executed successfully!")
}
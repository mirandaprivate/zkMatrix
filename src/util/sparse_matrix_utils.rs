/// Utility Functions to handle sparse matrices.
/// 
/// The element of the sparse matrix is in Zp.

use sprs::{CsMat, CsMatI, CsVec};
use super::zp_utils::{
    int_to_scalar, group_order, ZpElement, IntElement};


pub fn sprs_demo(){
    println!("********   sprs demo is started!  *******************");
    let mat1 = CsMat::new(
        (3, 3), 
        vec![0, 2, 4, 6], 
        vec![0, 2, 1, 2, 0, 1], 
        vec![1., 2., 3., 4., 5., 6.],
    );
    let mat2 = CsMat::new(
        (3, 3), 
        vec![0, 1, 3, 4], 
        vec![1, 0, 2, 2], 
        vec![7., 8., 9., 10.],
    );

    let product = &mat1 * &mat2;

    println!("product: {:?}", product);


    println!("*********************Sparse vector demo is started!**********************");

    let q: IntElement= group_order() - 1;
    println!("q: {:?}", q);
    let q_scalar: ZpElement = int_to_scalar(&q);
    println!("q_scalar: {:?}", q_scalar);
    let vec_data: Vec<ZpElement> = vec![q_scalar, q_scalar, q_scalar];

    let vec_indices : Vec<usize> = vec![0,2,4];

    let sparse_vec: sprs::CsVecBase<Vec<usize>, Vec<ZpElement>, ZpElement> = CsVec::new(5, vec_indices, vec_data);

    let sparse_vec_2: sprs::CsVecBase<Vec<usize>, Vec<ZpElement>, ZpElement> = CsVec::from(sparse_vec.clone());

    let q_3:ZpElement = q_scalar + q_scalar;

    println!("q_3: {:?}", q_3);

    // let sparse_vec_3 = sparse_vec + sparse_vec_2;

    // let sparse_vec_4 = sparse_vec * sparse_vec_2;

    let q_4:ZpElement = q_scalar * q_scalar;

    println!("q_4: {:?}", q_4);

    println!("Sparse vector: {:?}", sparse_vec);

    println!("Sparse vector demo is executed successfully!");
} 

//! Two-tier commitments to the matrices.
//! 
//! These matrices are stored in data/private.
//! Cache the first-tier commitments in data/private.
//! 
//! The commitments are published to data/public.
//! 
//! 
use std::time::Instant;

use crate::mat::Mat;
use crate::setup::SRS;

use crate::utils::curve::{G1Element, G2Element, GtElement};
use crate::utils::dirac;
use crate::utils::dirac::BraKet;
use crate::utils::to_file::FileIO;

pub trait CommitMat {
    fn commit_row_major(
        &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>
    ) -> (GtElement, Vec<G2Element>);

    fn commit_col_major(
        &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>
    ) -> (GtElement, Vec<G1Element>);

    fn commit_rm(&self, srs: &SRS) -> (GtElement, Vec<G2Element>) {
        self.commit_row_major(
            &srs.g_hat_vec, &srs.h_hat_vec)
    }

    fn commit_cm(&self, srs: &SRS) -> (GtElement, Vec<G1Element>) {
        self.commit_col_major(
            &srs.g_hat_vec, &srs.h_hat_vec)
    }
}

impl CommitMat for Mat<u64> {
    fn commit_row_major(
        &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>
    ) -> (GtElement, Vec<G2Element>) {
        
        let timer = Instant::now();

        let right_cache = self.ket(&h_base_vec);
               
        let result = dirac::inner_product(
            &g_base_vec, &right_cache);
    

        println!(" ** Time for matrix commitment: {:?}", timer.elapsed());

        right_cache.to_file(
            format!("{}_rp.cache", self.id), false)
        .unwrap();

        (result, right_cache.clone())

    }

    fn commit_col_major(
            &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>
    ) -> (GtElement, Vec<G1Element>) {

        let timer = Instant::now();

        let left_cache = self.bra(&g_base_vec);
        
        
        let result = dirac::inner_product(
            &left_cache, &h_base_vec);

        println!(" ** Time for matrix commitment: {:?}", timer.elapsed());

        left_cache.to_file(
                format!("{}_lp.cache", self.id), false)
            .unwrap();

        (result, left_cache.clone())
    }
}

impl CommitMat for Mat<i64> {
    fn commit_row_major(
        &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>
    ) -> (GtElement, Vec<G2Element>) {

        let timer = Instant::now();

        let right_cache = self.ket(&h_base_vec);
        
        let result = dirac::inner_product(
            &g_base_vec, &right_cache);
        
        println!(" ** Time for matrix commitment: {:?}", timer.elapsed());

        right_cache.to_file(
            format!("{}_rp.cache", self.id), false)
        .unwrap();
    
        (result, right_cache.clone())
    }

    fn commit_col_major(
            &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>
    ) -> (GtElement, Vec<G1Element>) {

        let timer = Instant::now();

        let left_cache = self.bra(&g_base_vec);
        
        
        let result = dirac::inner_product(
            &left_cache, &h_base_vec);
        
        println!(" ** Time for matrix commitment: {:?}", timer.elapsed());

        left_cache.to_file(
            format!("{}_lp.cache", self.id), false)
        .unwrap();
        (result, left_cache.clone())
    }
}

impl CommitMat for Mat<i128> {
    fn commit_row_major(
        &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>
    ) -> (GtElement, Vec<G2Element>) {

        let timer = Instant::now();

        let right_cache = self.ket(&h_base_vec);
       

        let result = dirac::inner_product(
            &g_base_vec, &right_cache);
        
        println!(" ** Time for matrix commitment: {:?}", timer.elapsed());

        right_cache.to_file(
            format!("{}_rp.cache", self.id), false)
        .unwrap();
        (result, right_cache.clone())

    }

    fn commit_col_major(
            &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>
    ) -> (GtElement, Vec<G1Element>) {

        let timer = Instant::now();

        let left_cache = self.bra(&g_base_vec);
        
        let result = dirac::inner_product(
            &left_cache, &h_base_vec);
        
        println!(" ** Time for matrix commitment: {:?}", timer.elapsed());

        left_cache.to_file(
            format!("{}_lp.cache", self.id), false)
        .unwrap();
        (result, left_cache.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::test_data::*;

    use crate::utils::curve::{G1Element, G2Element, GtElement};

    #[test]
    fn test_commit() {
        let mat_a: Mat<u64> = gen_mat_a_u64_direct_test();
        let vec_g: Vec<G1Element> = gen_vec_v_g1_direct_test();
        let vec_h: Vec<G2Element> = gen_vec_v_g2_direct_test();

        let right_proj: Vec<G2Element> = 
            gen_vec_va_g2_from_kronecker_test();
        let left_proj: Vec<G1Element> = 
            gen_vec_va_g1_from_kronecker_test();
        

        let gah: GtElement = gen_vav_gt_direct_test();

        let (gah_test_1, _) = mat_a
            .commit_row_major(&vec_g, &vec_h);
        let (gah_test_2, _)= mat_a
            .commit_col_major(&vec_g, &vec_h);

        let right_cache: Vec<G2Element> = FileIO::from_file(
                format!("{}_rp.cache", mat_a.id).to_string(),
                false,
            ).unwrap();
        let left_cache: Vec<G1Element> = FileIO::from_file(
               format!("{}_lp.cache", mat_a.id).to_string(), 
               false,
            ).unwrap();

        assert_eq!(gah_test_1, gah);
        assert_eq!(gah_test_2, gah);
        assert_eq!(right_cache, right_proj);
        assert_eq!(left_cache, left_proj);
    }
}
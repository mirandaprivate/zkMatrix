//! Commit the matrices for experiments.
//! 
//! Store the commitments in data/public.
//! 
//! Cache the first-tier commitments in data/private.
//! 

use crate::mat::Mat;
use crate::setup::SRS;

use crate::utils::curve::{G1Element, G2Element, GtElement};
use crate::utils::dirac;
use crate::utils::dirac::BraKet;
use crate::utils::to_file::FileIO;

pub trait CommitMat {
    fn commit_row_major(
        &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>) 
    -> GtElement;

    fn commit_col_major(
        &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>)
    -> GtElement;

    fn commit_rm(&self, srs: &SRS) -> GtElement {
        self.commit_row_major(
            &srs.g_hat_vec, &srs.h_hat_vec)
    }

    fn commit_cm(&self, srs: &SRS) -> GtElement {
        self.commit_col_major(
            &srs.g_hat_vec, &srs.h_hat_vec)
    }
}

impl CommitMat for Mat<u64> {
    fn commit_row_major(
        &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>) 
    -> GtElement {
        let right_cache = self.ket(&h_base_vec);
        right_cache.to_file(
            format!("{}_right_cache", self.id), false)
        .unwrap();

        let result = dirac::inner_product(
            &g_base_vec, &right_cache);
        result
    }

    fn commit_col_major(
            &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>) 
    -> GtElement {
        let left_cache = self.bra(&g_base_vec);
        left_cache.to_file(
            format!("{}_left_cache", self.id), false)
        .unwrap();
        
        let result = dirac::inner_product(
            &left_cache, &h_base_vec);
        result
    }
}

impl CommitMat for Mat<i64> {
    fn commit_row_major(
        &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>) 
    -> GtElement {
        let right_cache = self.ket(&h_base_vec);
        right_cache.to_file(
            format!("{}_right_cache", self.id), false)
        .unwrap();

        let result = dirac::inner_product(
            &g_base_vec, &right_cache);
        result
    }

    fn commit_col_major(
            &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>) 
    -> GtElement {
        let left_cache = self.bra(&g_base_vec);
        left_cache.to_file(
            format!("{}_left_cache", self.id), false)
        .unwrap();
        
        let result = dirac::inner_product(
            &left_cache, &h_base_vec);
        result
    }
}

impl CommitMat for Mat<i128> {
    fn commit_row_major(
        &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>) 
    -> GtElement {
        let right_cache = self.ket(&h_base_vec);
        right_cache.to_file(
            format!("{}_right_cache", self.id), false)
        .unwrap();

        let result = dirac::inner_product(
            &g_base_vec, &right_cache);
        result
    }

    fn commit_col_major(
            &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>) 
    -> GtElement {
        let left_cache = self.bra(&g_base_vec);
        left_cache.to_file(
            format!("{}_left_cache", self.id), false)
        .unwrap();
        
        let result = dirac::inner_product(
            &left_cache, &h_base_vec);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::test_data::*;

    use crate::utils::curve::{G1Element, G2Element, GtElement};

    #[test]
    fn test_commit() {
        let mat_a: Mat<u64> = {
            match FileIO::from_file(
                "a_test".to_string(), false)
            {
                Ok(mat) => mat,
                Err(_) => {
                    let mat_a: Mat<u64> = gen_mat_a_u64_direct_test();
                    mat_a.to_file(
                        mat_a.id.to_string(), false)
                    .unwrap();
                    mat_a
                }
            }
        };

        let vec_g = {
            match FileIO::from_file(
                "vec_g_test".to_string(), true)
            {
                Ok(vec) => vec,
                Err(_) => {
                    let vec_g: Vec<G1Element> = gen_vec_v_g1_direct_test();
                    vec_g.to_file(
                        "vec_g_test".to_string(), true)
                    .unwrap();
                    vec_g
                }
            }
        };

        let vec_h = {
            match FileIO::from_file(
                "vec_h_test".to_string(), true)
            {
                Ok(vec) => vec,
                Err(_) => {
                    let vec_h: Vec<G2Element> = gen_vec_v_g2_direct_test();
                    vec_h.to_file(
                        "vec_h_test".to_string(), true)
                    .unwrap();
                    vec_h
                }
            }
        };
        
        let right_proj = {
            match FileIO::from_file(
                format!("{}_right_test", mat_a.id).to_string(), false)
            {
                Ok(vec) => vec,
                Err(_) => {
                    let right_proj: Vec<G2Element> = gen_vec_va_g2_from_kronecker_test();
                    right_proj.to_file(
                        format!("{}_right_test", mat_a.id).to_string(), false)
                    .unwrap();
                    right_proj
                }
            }
        };

        let left_proj = {
            match FileIO::from_file(
                format!("{}_left_test", mat_a.id).to_string(), false)
            {
                Ok(vec) => vec,
                Err(_) => {
                    let left_proj: Vec<G1Element> 
                        = gen_vec_va_g1_from_kronecker_test();
                    left_proj.to_file(
                        format!("{}_left_test", mat_a.id).to_string(),
                        false,
                    ).unwrap();
                    left_proj
                }
            }
        };

        let gah: GtElement = gen_vav_gt_direct_test();

        let gah_test_1 = mat_a.commit_row_major(&vec_g, &vec_h);
        let gah_test_2= mat_a.commit_col_major(&vec_g, &vec_h);

        let right_cache: Vec<G2Element> = FileIO::from_file(
            format!("{}_right_cache", mat_a.id).to_string(), false)
        .unwrap();
        let left_cache: Vec<G1Element> = FileIO::from_file(
            format!("{}_left_cache", mat_a.id).to_string(), false)
        .unwrap();

        assert_eq!(gah_test_1, gah);
        assert_eq!(gah_test_2, gah);
        assert_eq!(right_cache, right_proj);
        assert_eq!(left_cache, left_proj);
    }
}
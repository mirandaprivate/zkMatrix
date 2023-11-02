//! Commit the matrices for experiments.
//! 
//! Store the commitments in data/public.
//! 
//! Cache the first-tier commitments in data/private.
//! 

use crate::curve::{G1Element, G2Element, GtElement};
use crate::dirac;
use crate::mat::{Mat,ToFile};

pub trait CommitMat {
    fn commit_row_major(
        &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>) 
    -> GtElement;

    fn commit_col_major(
        &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>)
    -> GtElement;
}

impl CommitMat for Mat<u64> {
    fn commit_row_major(
        &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>) 
    -> GtElement {
        let right_cache = dirac:: ket(&self, &h_base_vec);
        right_cache.to_file(
            format!("{}_right_cache", self.id), false)
        .unwrap();

        let result = dirac::inner_product(&g_base_vec, &right_cache);
        result
    }

    fn commit_col_major(
            &self, g_base_vec: &Vec<G1Element>, h_base_vec: &Vec<G2Element>) 
    -> GtElement {
        let left_cache = dirac::bra(&self, &g_base_vec);
        left_cache.to_file(
            format!("{}_left_cache", self.id), false)
        .unwrap();
        
        let result = dirac::inner_product(&left_cache, &h_base_vec);
        result
    }
}
    
#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_data::*;

    use crate::curve::{G1Element, G2Element, GtElement};

    #[test]
    fn test_commit() {
        let mat_a: Mat<u64> = {
            match ToFile::from_file(
                "a".to_string(), false)
            {
                Ok(mat) => mat,
                Err(_) => {
                    let mat_a: Mat<u64> = gen_mat_a_u64_direct();
                    mat_a.to_file(
                        mat_a.id.to_string(), false)
                    .unwrap();
                    mat_a
                }
            }
        };

        let vec_g = {
            match ToFile::from_file(
                "vec_g".to_string(), true)
            {
                Ok(vec) => vec,
                Err(_) => {
                    let vec_g: Vec<G1Element> = gen_vec_v_g1_direct();
                    vec_g.to_file(
                        "vec_g".to_string(), true)
                    .unwrap();
                    vec_g
                }
            }
        };

        let vec_h = {
            match ToFile::from_file(
                "vec_h".to_string(), true)
            {
                Ok(vec) => vec,
                Err(_) => {
                    let vec_h: Vec<G2Element> = gen_vec_v_g2_direct();
                    vec_h.to_file(
                        "vec_h".to_string(), true)
                    .unwrap();
                    vec_h
                }
            }
        };
        
        let right_proj = {
            match ToFile::from_file(
                format!("{}_right", mat_a.id).to_string(), false)
            {
                Ok(vec) => vec,
                Err(_) => {
                    let right_proj: Vec<G2Element> = gen_vec_va_g2_from_kronecker();
                    right_proj.to_file(
                        format!("{}_right", mat_a.id).to_string(), false)
                    .unwrap();
                    right_proj
                }
            }
        };

        let left_proj = {
            match ToFile::from_file(
                format!("{}_left", mat_a.id).to_string(), false)
            {
                Ok(vec) => vec,
                Err(_) => {
                    let left_proj: Vec<G1Element> = gen_vec_va_g1_from_kronecker();
                    left_proj.to_file(
                        format!("{}_left", mat_a.id).to_string(), false)
                    .unwrap();
                    left_proj
                }
            }
        };

        let gah: GtElement = gen_vav_gt_direct();

        let gah_test_1 = mat_a.commit_row_major(&vec_g, &vec_h);
        let gah_test_2= mat_a.commit_col_major(&vec_g, &vec_h);

        let right_cache: Vec<G2Element> = ToFile::from_file(
            format!("{}_right_cache", mat_a.id).to_string(), false)
        .unwrap();
        let left_cache: Vec<G1Element> = ToFile::from_file(
            format!("{}_left_cache", mat_a.id).to_string(), false)
        .unwrap();

        assert_eq!(gah_test_1, gah);
        assert_eq!(gah_test_2, gah);
        assert_eq!(right_cache, right_proj);
        assert_eq!(left_cache, left_proj);
    }
}
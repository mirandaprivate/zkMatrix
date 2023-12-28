//! Implementation of the Scalar Projection protocol
//!
//! Details of this protocol can be found in the DualMatrix paper 
//!
//! To prove that holding a secret matrix \bm{a} 
//! and two Zp elements c_tilde and a_tilde such that
//! 
//! C_a = < \vec{G}, \bm{a}, \vec{H} > + a_tilde * blind_base,
//! C_c = (l^T \bm{a} r) e(\hat{G}, \hat{H}) + c_tilde * blind_base,
// 
use crate::mat::Mat;
use crate::setup::SRS;

use crate::utils::curve::{
    ZpElement, GtElement, G1Element, 
    ConvertToZp,
};
use crate::utils::dirac::{self, BraKetZp};
use crate::utils::fiat_shamir::{TranElem, TranSeq};
use crate::utils::xi;

use crate::protocols::pip::{PipG1, PipG2};

use crate::zkprotocols::zk_trans::ZkTranSeqProver;
use crate::zkprotocols::zk_scalars::{ZkSchnorr, ZkSemiMulScalar};



/// Interface when l_vec and r_vec are arbitrary public vectors
pub struct ZkScalarProj {
    pub c_com: GtElement,
    pub a_com: GtElement,
    pub shape: (usize, usize),
    pub l_vec: Vec<ZpElement>,
    pub r_vec: Vec<ZpElement>,
}

/// Interface when y_l and y_r are vectors generated from a random y
/// In this case, the number of field operations can be optimized
pub struct ZkScalarProjPoly {
    pub c_com: GtElement,
    pub a_com: GtElement,
    pub shape: (usize, usize),
    pub y: ZpElement, 
}

pub trait ZkScalarProjInterface {
    fn get_c_com(&self) -> GtElement;
    fn get_a_com(&self) -> GtElement;
    fn get_shape(&self) -> (usize, usize);
    fn reduce_l(&self, challenges: &Vec<ZpElement>) -> ZpElement;
    fn reduce_r(&self, challenges: &Vec<ZpElement>) -> ZpElement;
    fn get_l_vec(&self) -> Vec<ZpElement>;
    fn get_r_vec(&self) -> Vec<ZpElement>;
}

impl ZkScalarProjInterface for ZkScalarProj {
    fn reduce_l(&self, challenges: &Vec<ZpElement>) -> ZpElement {
        xi::reduce_from_challenges(challenges, &self.get_l_vec())
    }

    fn reduce_r(&self, challenges: &Vec<ZpElement>) -> ZpElement {
        xi::reduce_from_challenges(challenges, &self.get_r_vec())
    }

    fn get_l_vec(&self) -> Vec<ZpElement> {
        self.l_vec.clone()
    }

    fn get_r_vec(&self) -> Vec<ZpElement> {
        self.r_vec.clone()
    }

    fn get_c_com(&self) -> GtElement {
        self.c_com
    }

    fn get_a_com(&self) -> GtElement {
        self.a_com
    }

    fn get_shape(&self) -> (usize, usize) {
        self.shape
    }

}

impl ZkScalarProjInterface for ZkScalarProjPoly {
    fn reduce_l(&self, challenges: &Vec<ZpElement>) -> ZpElement {
        let n = self.get_shape().1;

        xi::phi_s(
            self.y, 
            challenges,
            0 as usize, 
            n as usize,
        )

    }

    fn reduce_r(&self, challenges: &Vec<ZpElement>) -> ZpElement {
        xi::phi_s(
            self.y, challenges,
            0 as usize, 
            1 as usize,
        )
    }

    fn get_l_vec(&self) -> Vec<ZpElement> {
        let y = self.y;
        let m = self.get_shape().0;
        let n = self.get_shape().1;
        let step = y.pow(n as u64);
        
        
        std::iter::successors(
            Some(ZpElement::from(1 as u64)), 
            |&x| Some(x * step)
        ).take(m).collect::<Vec<ZpElement>>()
    }

    fn get_r_vec(&self) -> Vec<ZpElement> {
        let y = self.y;
        let n = self.get_shape().1;
        
        std::iter::successors(
            Some(ZpElement::from(1 as u64)),
            |&x| Some(y * x),
        ).take(n).collect::<Vec<ZpElement>>()
    }

    fn get_c_com(&self) -> GtElement {
        self.c_com
    }

    fn get_a_com(&self) -> GtElement {
        self.a_com
    }

    fn get_shape(&self) -> (usize, usize) {
        self.shape
    }

}

impl ZkScalarProj {
    pub fn new(
        c_com_value: GtElement, 
        a_com_value: GtElement,  
        shape_value: (usize, usize),
        l_vec: &Vec<ZpElement>,
        r_vec: &Vec<ZpElement>,
    ) -> Self {
        Self {
            c_com: c_com_value,
            a_com: a_com_value,
            shape: shape_value,
            l_vec: l_vec.clone(),
            r_vec: r_vec.clone(),
        }
    }
}

impl ZkScalarProjPoly {
    pub fn new(
        c_com_value: GtElement, 
        a_com_value: GtElement,  
        shape_value: (usize, usize),
        y_value: ZpElement,
    ) -> Self {
        Self {
            c_com: c_com_value,
            a_com: a_com_value,
            shape: shape_value,
            y: y_value,
        }
    }
}

impl ZkScalarProjCmProof for ZkScalarProj {}
impl ZkScalarProjCmProof for ZkScalarProjPoly {}

pub trait ZkScalarProjCmProof: ZkScalarProjInterface {

    fn prove_cm<T>(
        &self, srs: &SRS, 
        zk_trans_seq: &mut ZkTranSeqProver, 
        mat_a: &Mat<T>, 
        a_cache: &Vec<G1Element>,
        c_tilde: ZpElement,
        a_tilde: ZpElement,
    ) where 
        T: 'static + ConvertToZp,
        Mat<T>: 'static + BraKetZp, 
    {

        zk_trans_seq.push_without_blinding(TranElem::Gt(self.get_c_com()));
        zk_trans_seq.push_without_blinding(TranElem::Gt(self.get_a_com()));
        zk_trans_seq.push_without_blinding(TranElem::Size(self.get_shape().0));
        zk_trans_seq.push_without_blinding(TranElem::Size(self.get_shape().1));

        let x = zk_trans_seq.gen_challenge();

        let m = self.get_shape().0;
        let n = self.get_shape().1;

        if m & (m - 1) != 0 || n & (n - 1) != 0 
            || m != mat_a.shape.0 || n != mat_a.shape.1 
        {
            panic!("Invalid shape when proving ScalarProjCm");
        }

        let log_m = (m as u64).ilog2() as usize;
        let log_n = (n as u64).ilog2() as usize;

        let u_0 = srs.g_hat * srs.h_hat;
        let l_vec = self.get_l_vec();
        let r_vec = self.get_r_vec();

        let la = mat_a.bra_zp(&l_vec);

        let mut v_current = dirac::vec_convert_to_zp_vec(
            &la);
        let mut capital_a_current = a_cache[0..n].to_vec();
        let mut h_vec_current = srs.h_hat_vec[0..n].to_vec();
        let mut r_current = r_vec[0..n].to_vec();
        
        let mut challenges_n: Vec<ZpElement> = Vec::new();
        let mut challenges_inv_n: Vec<ZpElement> = Vec::new();
 

        for j in 0..log_n {
            let current_len = n / 2usize.pow(j as u32);
            
            let v_left = 
                v_current[0..current_len/2].to_vec();
            let v_right = 
                v_current[current_len/2..current_len].to_vec();

            let capital_a_left = 
                capital_a_current[0..current_len/2].to_vec();
            let capital_a_right = 
                capital_a_current[current_len/2..current_len].to_vec();
            
            let r_left = 
                r_current[0..current_len/2].to_vec();
            let r_right = 
                r_current[current_len/2..current_len].to_vec();
            

            let h_left = 
                h_vec_current[0..current_len/2].to_vec();
            let h_right = 
                h_vec_current[current_len/2..current_len].to_vec();

            let l_tr = 
                dirac::inner_product(&capital_a_left, &h_right)
                + u_0 * x * dirac::inner_product(&v_left, &r_right);
            let r_tr = 
                dirac::inner_product(&capital_a_right, &h_left)
                + u_0 * x * dirac::inner_product(&v_right, &r_left);

            zk_trans_seq.push_gen_blinding(TranElem::Gt(l_tr));
            zk_trans_seq.push_gen_blinding(TranElem::Gt(r_tr));
            
            let x_j = zk_trans_seq.gen_challenge();
            let x_j_inv = x_j.inv();

            challenges_n.push(x_j);
            challenges_inv_n.push(x_j_inv);

            v_current = dirac::vec_addition(
                &v_left,
                &dirac::vec_scalar_mul(
                    &v_right, x_j_inv),
            );

            capital_a_current = dirac::vec_addition(
                &capital_a_left,
                &dirac::vec_scalar_mul(
                    &capital_a_right, x_j_inv),
            );

            h_vec_current = dirac::vec_addition(
                &h_left,
                &dirac::vec_scalar_mul(
                    &h_right, x_j),
            );

            r_current = dirac::vec_addition(
                &r_left,
                &dirac::vec_scalar_mul(
                    &r_right, x_j),
            );

        }

        let xi_n_inv = xi::xi_from_challenges(&challenges_inv_n);
        let a_xi_inv = mat_a.ket_zp(&xi_n_inv);


        let h_reduce = h_vec_current[0];
        let r_reduce = r_current[0];

        let mut a_current = dirac::vec_convert_to_zp_vec(
            &a_xi_inv);
        let mut g_vec_current = srs.g_hat_vec[0..m].to_vec();
        let mut l_current = l_vec[0..m].to_vec();
        
        let mut challenges_m: Vec<ZpElement> = Vec::new();
        let mut challenges_inv_m: Vec<ZpElement> = Vec::new();
        

        for j in 0..log_m {
            let current_len = m / 2usize.pow(j as u32);
            
            let a_left = 
                a_current[0..current_len/2].to_vec();
            let a_right = 
                a_current[current_len/2..current_len].to_vec();
            
            let l_left = 
                l_current[0..current_len/2].to_vec();
            let l_right = 
                l_current[current_len/2..current_len].to_vec();
            

            let g_left = 
                g_vec_current[0..current_len/2].to_vec();
            let g_right = 
                g_vec_current[current_len/2..current_len].to_vec();

            let l_tr = 
                h_reduce * dirac::inner_product(&a_left, &g_right)
                + u_0 * x * r_reduce 
                * dirac::inner_product(&a_left, &l_right);
            let r_tr = 
                h_reduce * dirac::inner_product(&a_right, &g_left)
                + u_0 * x * r_reduce 
                * dirac::inner_product(&a_right, &l_left);

            zk_trans_seq.push_gen_blinding(TranElem::Gt(l_tr));
            zk_trans_seq.push_gen_blinding(TranElem::Gt(r_tr));
            
            let x_j = zk_trans_seq.gen_challenge();
            let x_j_inv = x_j.inv();

            challenges_m.push(x_j);
            challenges_inv_m.push(x_j_inv);

            a_current = dirac::vec_addition(
                &a_left,
                &dirac::vec_scalar_mul(
                    &a_right, x_j_inv),
            );

            g_vec_current = dirac::vec_addition(
                &g_left,
                &dirac::vec_scalar_mul(
                    &g_right, x_j),
            );

            l_current = dirac::vec_addition(
                &l_left,
                &dirac::vec_scalar_mul(
                    &l_right, x_j),
            );

        }

        let a_reduce = a_current[0];

        let g_reduce = g_vec_current[0];


        zk_trans_seq.push_without_blinding(TranElem::G1(g_reduce));
        zk_trans_seq.push_without_blinding(TranElem::G2(h_reduce));

        let pip_g1 = PipG1::new(
            g_reduce, &challenges_m
        );

        let pip_g2 = PipG2::new(
            h_reduce, &challenges_n
        );


        // /////////////////////////////////////////////////////////////
        // Add Zero-Knowledge from now on
        // /////////////////////////////////////////////////////////////
        

        let xi_l = self.reduce_l(&challenges_m);
        let xi_r = self.reduce_r(&challenges_n);

        let (a_reduce_blind, a_reduce_tilde) = 
            zk_trans_seq.push_gen_blinding(TranElem::Zp(a_reduce));
        
            
        let base_rhs =  
            x * xi_l * xi_r * srs.g_hat * srs.h_hat
            + g_reduce * h_reduce;
        
        let rhs_com = 
            base_rhs * a_reduce;

        let (rhs_blind, rhs_tilde) = 
            zk_trans_seq.push_gen_blinding( 
            TranElem::Gt(rhs_com)
            );
        
        let mut lhs_tilde = x * c_tilde + a_tilde;

        let length = zk_trans_seq.blind_seq.len();
        let mut current_index = length - 2 - 2 * log_n - 2 * log_m;

        // assert_eq!(current_index, 0);

        for j in 0..log_n {
            let l_tilde = zk_trans_seq.blind_seq[current_index];
            let r_tilde = zk_trans_seq.blind_seq[current_index+1];
            let x_j = challenges_n[j];
            let x_j_inv = challenges_inv_n[j];
            lhs_tilde = lhs_tilde + l_tilde * x_j + r_tilde * x_j_inv;
            current_index +=2;
        }  

        for j in 0..log_m {
            let l_tilde = zk_trans_seq.blind_seq[current_index];
            let r_tilde = zk_trans_seq.blind_seq[current_index+1];
            let x_j = challenges_m[j];
            let x_j_inv = challenges_inv_m[j];
            lhs_tilde = lhs_tilde + l_tilde * x_j + r_tilde * x_j_inv;
            current_index +=2;
        }  

        let eq_tilde = lhs_tilde - rhs_tilde;
        let blind_base = srs.blind_base;
        let eq_tilde_com = eq_tilde * blind_base;

        let zk_semi_mul = ZkSemiMulScalar::new(
            rhs_blind,
            a_reduce_blind,
            base_rhs,
        );

        let zk_schnorr = ZkSchnorr::new(
            eq_tilde_com,
            blind_base,
        );


        pip_g1.prove(srs, zk_trans_seq.get_mut_trans_seq());
        pip_g2.prove(srs, zk_trans_seq.get_mut_trans_seq());
        zk_semi_mul.prove(
            srs, zk_trans_seq.get_mut_trans_seq(), 
            a_reduce, rhs_tilde, a_reduce_tilde
        );
        zk_schnorr.prove(
            zk_trans_seq.get_mut_trans_seq(), 
            eq_tilde,
        );

    }

    fn verify_as_subprotocol_cm(
        &self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        let pointer_old = trans_seq.pointer;
        
        if (
            TranElem::Gt(self.get_c_com()),
            TranElem::Gt(self.get_a_com()),
            TranElem::Size(self.get_shape().0),
            TranElem::Size(self.get_shape().1)
        ) != (
            trans_seq.data[pointer_old],
            trans_seq.data[pointer_old + 1],
            trans_seq.data[pointer_old + 2], 
            trans_seq.data[pointer_old + 3],
        ) {
            println!("!! Invalid public input when verifying ScalarProjCm");
            return false;
        } 

        let x: ZpElement;
        
        if let TranElem::Coin(x_read) = 
            trans_seq.data[pointer_old + 4] 
        {
            x = x_read;
        } else {
            println!("!! Invalid transcript when verifying ScalarProjCm");
            return false;
        }

        let m = self.get_shape().0;
        let n = self.get_shape().1;
        if m & (m - 1) != 0 || n & (n - 1) != 0 {
            panic!("Invalid shape when verifying ScalarProjCm");
        }
        let log_m = (m as u64).ilog2() as usize;
        let log_n = (n as u64).ilog2() as usize;

        trans_seq.pointer = 
            pointer_old + 5 + 3* log_m + 3* log_n + 4 ;

        let mut current_pointer = pointer_old + 5;
        let mut lhs: GtElement = 
            x * self.get_c_com().clone() + self.get_a_com().clone();
        
        let mut challenges_n: Vec<ZpElement> = Vec::new();
        let mut challenges_inv_n: Vec<ZpElement> = Vec::new();

        for _ in 0..log_n {

            if let (
                TranElem::Gt(l_tr),
                TranElem::Gt(r_tr),
                TranElem::Coin(x_j),
            ) = (
                trans_seq.data[current_pointer],
                trans_seq.data[current_pointer + 1],
                trans_seq.data[current_pointer + 2],
            ) {
                
                let x_j_inv = x_j.inv();
                lhs = lhs + l_tr * x_j + r_tr * x_j_inv;
                challenges_n.push(x_j);
                challenges_inv_n.push(x_j_inv);

            } else {
                println!("!! Invalid transcript when verifying ScalarProjCm");
                return false;
            }

            current_pointer += 3;
        }

        let mut challenges_m: Vec<ZpElement> = Vec::new();
        let mut challenges_inv_m: Vec<ZpElement> = Vec::new();

        for _ in 0..log_m {

            if let (
                TranElem::Gt(l_tr),
                TranElem::Gt(r_tr),
                TranElem::Coin(x_j),
            ) = (
                trans_seq.data[current_pointer],
                trans_seq.data[current_pointer + 1],
                trans_seq.data[current_pointer + 2],
            ) {
                
                let x_j_inv = x_j.inv();
                lhs = lhs + l_tr * x_j + r_tr * x_j_inv;
                challenges_m.push(x_j);
                challenges_inv_m.push(x_j_inv);

            } else {
                println!("!! Invalid transcript when verifying ScalarProjCm");
                return false;
            }

            current_pointer += 3;
        }


        let xi_l = self.reduce_l(&challenges_m);
        let xi_r = self.reduce_r(&challenges_n);

        if let (
            TranElem::G1(g_reduce),
            TranElem::G2(h_reduce),
            TranElem::Gt(a_reduce_blind),
            TranElem::Gt(rhs_blind),
        ) = (
            trans_seq.data[current_pointer],
            trans_seq.data[current_pointer+1],
            trans_seq.data[current_pointer+2],
            trans_seq.data[current_pointer+3],

        ) {

            let eq_tilde_com = lhs - rhs_blind;
            
            let rhs_base = x * xi_l * xi_r * srs.g_hat * srs.h_hat
                + g_reduce * h_reduce;
            let blind_base = srs.blind_base;

            let pip_g1 = PipG1::new(
                g_reduce, &challenges_m
            );
    
            let pip_g2 = PipG2::new(
                h_reduce, &challenges_n
            );


            let zk_semi_mul = ZkSemiMulScalar::new(
                rhs_blind,
                a_reduce_blind,
                rhs_base,
            );
    
            let zk_schnorr = ZkSchnorr::new(
                eq_tilde_com,
                blind_base,
            );

            let check1 = pip_g1.verify_as_subprotocol(srs, trans_seq);
            let check2 = pip_g2.verify_as_subprotocol(srs, trans_seq);
            let check3 = zk_semi_mul.verify_as_subprotocol::<GtElement, ZpElement>(
                srs, trans_seq
            );
            let check4 = zk_schnorr.verify_as_subprotocol(trans_seq);

            // println!(" check 1 {:?}", check1);
            // println!(" check 2 {:?}", check2);
            // println!(" check 3 {:?}", check3);
            // println!(" check 4 {:?}", check4);
            return check1 && check2 && check3 && check4;

        } else {
            println!("!! Invalid transcript when verifying ScalarProjCm");
            return false;
        }

    }

    fn verify_cm(&self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        if trans_seq.check_fiat_shamir() == false {
            println!("!! Fiat shamir check failed when verifying ScalarProjCm");
            return false;
        }

        self.verify_as_subprotocol_cm(srs, trans_seq)
    }

}


#[cfg(test)]
mod tests {
    
    use super::*;

    use rand::Rng;

    use crate::commit_mat::CommitMat;


    #[test]
    fn test_zk_scalar_proj() {

        let srs = SRS::new(64);

        let m = 16 as usize;
        let n = 32 as usize;

        let a_data = (0..m).into_iter()
            .zip((0..n).into_iter())
            .map(|(i, j)| 
                (i, j, rand::thread_rng().gen::<u32>() as i64 )
            ).collect::<Vec<(usize, usize, i64)>>();

        
        let a = Mat::new_from_data_vec(
            "a_test", 
            (m, n), 
            a_data
        );

        let (a_com, a_cache_cm) 
            = a.commit_col_major(
                &srs.g_hat_vec, &srs.h_hat_vec
            );
    

        let l_vec = (0..m).map(|_| 
            ZpElement::rand()
        ).collect::<Vec<ZpElement>>();

        let r_vec = (0..n).map(|_| 
            ZpElement::rand()
        ).collect::<Vec<ZpElement>>();

        let c = a.braket_zp(&l_vec, &r_vec);

        let c_com = c * srs.g_hat * srs.h_hat;

        let a_tilde = ZpElement::rand();
        let c_tilde = ZpElement::rand();

        let c_blind = c_com + c_tilde * srs.blind_base;
        let a_blind = a_com + a_tilde * srs.blind_base;

        let scalar_proj_protocol = ZkScalarProj::new(
            c_blind,
            a_blind, 
            (m, n), 
            &l_vec, 
            &r_vec,
        );

        let mut zk_trans_seq = ZkTranSeqProver::new(&srs);

        scalar_proj_protocol.prove_cm::<i64>(
            &srs, 
            &mut zk_trans_seq, 
            &a,
            &a_cache_cm,
            c_tilde,
            a_tilde,
        );

        let mut trans_seq = zk_trans_seq.publish_trans();

        let result_cm = scalar_proj_protocol.verify_cm(
            &srs, &mut trans_seq
        );

        assert_eq!(result_cm, true);

        println!(" * Verification of ScalarProjCm passed");


    }

    
}

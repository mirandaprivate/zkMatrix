//! Zero-Knowledge implementation of the Inner Product protocol in Gt
//!
//! Details of this protocol can be found in the DualMatrix paper 
//!
//! To prove that holding two vectors vec{a} and vec{b}
//! and c_tilde such that
//! gt_com = <vec{a}, vec{b}> e(\hat{G}, \hat{H}) 
//!         + e(\hat{G}, <vec{a}, vec{H}>) + e(<vec{b}, \hat{G}>, \hat{H})
//!         + c_tilde * blind_base
//!
use crate::setup::SRS;

use crate::utils::curve::{ZpElement, GtElement, ConvertToZp, G2Element, G1Element};
use crate::utils::dirac;
use crate::utils::fiat_shamir::{TranElem, TranSeq};

use crate::protocols::pip::{PipG1, PipG2};

use crate::zkprotocols::zk_trans::ZkTranSeqProver;
use crate::zkprotocols::zk_scalars::{ZkSchnorr, ZkMulScalar, ZkSemiMulScalar};


pub struct ZkIpGt {
    pub gt_com: GtElement,
    pub length: usize,
}

impl ZkIpGt {
    pub fn new(gt_com_value: GtElement, length_value: usize
    ) -> Self {
        Self {
            gt_com: gt_com_value,
            length: length_value,
        }
    }

    pub fn prove<T>(
        &self, srs: &SRS, 
        zk_trans_seq: &mut ZkTranSeqProver, 
        vec_a: &Vec<T>, 
        vec_b: &Vec<T>,
        tilde_c: ZpElement,
    ) where T: 'static + ConvertToZp {

        zk_trans_seq.push_without_blinding(
            TranElem::Gt(self.gt_com),
        );
        zk_trans_seq.push_without_blinding(
            TranElem::Size(self.length),
        );

        if (self.length & (self.length - 1)) != 0 {
            panic!("Length is not a power of 2 when proving IpGt");
        }

        let n = self.length;
        let log_n = (n as u64).ilog2() as usize;

        let mut vec_a_current = dirac::vec_convert_to_zp_vec(
            vec_a);
        let mut vec_b_current = dirac::vec_convert_to_zp_vec(
            vec_b);
        let mut g_vec_current = srs.g_hat_vec[0..n].to_vec();
        let mut h_vec_current = srs.h_hat_vec[0..n].to_vec();
        
        let mut challenges: Vec<ZpElement> = Vec::new();
        let mut challenges_inv: Vec<ZpElement> = Vec::new();
        let u_0 = srs.g_hat * srs.h_hat;
        let g_0 = srs.g_hat;
        let h_0 = srs.h_hat; 

        for j in 0..log_n {
            let current_len = n / 2usize.pow(j as u32);
            
            let vec_a_left = 
                vec_a_current[0..current_len/2].to_vec();
            let vec_a_right = 
                vec_a_current[current_len/2..current_len].to_vec();
            
            let vec_b_left = 
                vec_b_current[0..current_len/2].to_vec();
            let vec_b_right = 
                vec_b_current[current_len/2..current_len].to_vec();
            
            let g_left = 
                g_vec_current[0..current_len/2].to_vec();
            let g_right = 
                g_vec_current[current_len/2..current_len].to_vec();

            let h_left = 
                h_vec_current[0..current_len/2].to_vec();
            let h_right = 
                h_vec_current[current_len/2..current_len].to_vec();

            let l_tr = 
                g_0 * dirac::inner_product(&vec_a_left, &h_right)
                + h_0 * dirac::inner_product(&vec_b_right, &g_left)
                + u_0 * dirac::inner_product(&vec_a_left, &vec_b_right);
            let r_tr = 
                g_0 * dirac::inner_product(&vec_a_right, &h_left)
                + h_0 * dirac::inner_product(&vec_b_left, &g_right)
                + u_0 * dirac::inner_product(&vec_a_right, &vec_b_left);

            zk_trans_seq.push_gen_blinding(TranElem::Gt(l_tr));
            zk_trans_seq.push_gen_blinding(TranElem::Gt(r_tr));
            
            let x_j = zk_trans_seq.gen_challenge();
            let x_j_inv = x_j.inv();

            challenges.push(x_j);
            challenges_inv.push(x_j_inv);

            vec_a_current = dirac::vec_addition(
                &vec_a_left,
                &dirac::vec_scalar_mul(
                    &vec_a_right, x_j_inv),
            );

            vec_b_current = dirac::vec_addition(
                &vec_b_left,
                &dirac::vec_scalar_mul(
                    &vec_b_right, x_j),
            );

            g_vec_current = dirac::vec_addition(
                &g_left,
                &dirac::vec_scalar_mul(
                    &g_right, x_j_inv),
            );

            h_vec_current = dirac::vec_addition(
                &h_left,
                &dirac::vec_scalar_mul(
                    &h_right, x_j),
            );

        }

        let a_reduce = vec_a_current[0];
        let b_reduce = vec_b_current[0];
        let g_reduce = g_vec_current[0];
        let h_reduce = h_vec_current[0];

        let a_h = a_reduce * h_reduce;
        let b_g = b_reduce * g_reduce;
        let a_b = a_reduce * b_reduce;

        let (a_reduce_blind, a_reduce_tilde) =
            zk_trans_seq.push_gen_blinding(TranElem::Zp(a_reduce));
        let (b_reduce_blind, b_reduce_tilde) =
            zk_trans_seq.push_gen_blinding(TranElem::Zp(b_reduce));
        zk_trans_seq.push_without_blinding(TranElem::G1(g_reduce));
        zk_trans_seq.push_without_blinding(TranElem::G2(h_reduce));


        let pip_g1 = PipG1::new(
            g_reduce, &challenges_inv
        );

        let pip_g2 = PipG2::new(
            h_reduce, &challenges
        );


        // ///////////////////////////////////////////////////////////
        // Add zk from now on
        // //////////////////////////////////////////////////////////
        let (a_h_blind, a_h_tilde) = 
            zk_trans_seq.push_gen_blinding(TranElem::G2(a_h));
        let (b_g_blind, b_g_tilde) =
            zk_trans_seq.push_gen_blinding(TranElem::G1(b_g));
        let (a_b_blind, a_b_tilde) =
            zk_trans_seq.push_gen_blinding(TranElem::Zp(a_b));        

        let length = zk_trans_seq.blind_seq.len();
        
        let rhs_tilde = a_b_tilde + a_h_tilde + b_g_tilde; 
           
       

        let mut current_index = length - 5 - 2 * log_n;
        let mut lhs_tilde = tilde_c;
        for j in 0..log_n {
            let l_tilde = zk_trans_seq.blind_seq[current_index];
            let r_tilde = zk_trans_seq.blind_seq[current_index+1];
            let x_j = challenges[j];
            let x_j_inv = challenges_inv[j];
            lhs_tilde = lhs_tilde +  l_tilde * x_j + r_tilde * x_j_inv;
            
            current_index += 2; 
        }

        let eq_tilde = lhs_tilde - rhs_tilde;
        let blind_base = srs.blind_base;
        let eq_tilde_com = eq_tilde * blind_base;

        let zk_semi_mul_1 = ZkSemiMulScalar::new(
            a_h_blind,
            a_reduce_blind,
            h_reduce,
        );

        let zk_semi_mul_2 = ZkSemiMulScalar::new(
            b_g_blind,
            b_reduce_blind,
            g_reduce,
        );

        let zk_mul_scalar = ZkMulScalar::new(
            a_b_blind,
            a_reduce_blind,
            b_reduce_blind,
        );

        pip_g1.prove(srs, zk_trans_seq.get_mut_trans_seq());
        pip_g2.prove(srs, zk_trans_seq.get_mut_trans_seq());

        let zk_schnorr = ZkSchnorr::new(
            eq_tilde_com,
            blind_base,
        );

        zk_semi_mul_1.prove(
            srs,
            zk_trans_seq.get_mut_trans_seq(),
            a_reduce,
            a_h_tilde,
            a_reduce_tilde,
        );

        zk_semi_mul_2.prove(
            srs,
            zk_trans_seq.get_mut_trans_seq(),
            b_reduce,
            b_g_tilde,
            b_reduce_tilde,
        );

        zk_mul_scalar.prove(
            srs, 
            zk_trans_seq.get_mut_trans_seq(), 
            a_reduce, b_reduce, 
            a_b_tilde, a_reduce_tilde, b_reduce_tilde,
        );

        zk_schnorr.prove(
            zk_trans_seq.get_mut_trans_seq(), 
            eq_tilde,
        )

    }

    pub fn verify_as_subprotocol(
        &self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        let pointer_old = trans_seq.pointer;
        
        if (
            TranElem::Gt(self.gt_com),
            TranElem::Size(self.length),
        ) != (
            trans_seq.data[pointer_old],
            trans_seq.data[pointer_old + 1],
        ) {
            println!("{:?}", self.gt_com);
            println!("{:?}", trans_seq.data[pointer_old]);
            println!("!! Invalid public input when verifying IpGt");
            return false;
        } 

        let n = self.length;
        let log_n = (n as u64).ilog2() as usize;

        trans_seq.pointer = pointer_old + 3 * log_n + 9;

        let mut current_pointer = pointer_old + 2;
        let mut lhs: GtElement = self.gt_com.clone();
        
        let mut challenges: Vec<ZpElement> = Vec::new();
        let mut challenges_inv: Vec<ZpElement> = Vec::new();

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
                challenges.push(x_j);
                challenges_inv.push(x_j_inv);

            } else {
                println!("!! Invalid transcript when verifying IpGt");
                return false;
            }

            current_pointer += 3;
        }

        if let (
            TranElem::Gt(a_reduce_blind),
            TranElem::Gt(b_reduce_blind),
            TranElem::G1(g_reduce),
            TranElem::G2(h_reduce),
            TranElem::Gt(a_h_blind),
            TranElem::Gt(b_g_blind),
            TranElem::Gt(a_b_blind),
        ) = (
            trans_seq.data[current_pointer],
            trans_seq.data[current_pointer+1],
            trans_seq.data[current_pointer+2],
            trans_seq.data[current_pointer+3],
            trans_seq.data[current_pointer+4],
            trans_seq.data[current_pointer+5],
            trans_seq.data[current_pointer+6],
        ) {
            let rhs = 
                a_h_blind + b_g_blind + a_b_blind;
            
            let eq_com = lhs - rhs;
                
            let pip_g1 = PipG1::new(
                g_reduce, &challenges_inv
            );
    
            let pip_g2 = PipG2::new(
                h_reduce, &challenges
            );

            let check1 = pip_g1.verify_as_subprotocol(srs, trans_seq);
            let check2 = pip_g2.verify_as_subprotocol(srs, trans_seq);

            // ///////////////////////////////////////////////////////////
            // Add zk from now on
            // /////////////////////////////////////////////////////////////
            let blind_base = srs.blind_base;

            let zk_semi_mul_1 = ZkSemiMulScalar::new(
                a_h_blind,
                a_reduce_blind,
                h_reduce,
            );
    
            let zk_semi_mul_2 = ZkSemiMulScalar::new(
                b_g_blind,
                b_reduce_blind,
                g_reduce,
            );
    
            let zk_mul_scalar = ZkMulScalar::new(
                a_b_blind,
                a_reduce_blind,
                b_reduce_blind,
            );
    
            let zk_schnorr = ZkSchnorr::new(
                eq_com,
                blind_base,
            );
    
            let check3 = zk_semi_mul_1.verify_as_subprotocol::<G2Element, ZpElement>(
                srs,
                trans_seq,
            );
    
            let check4 = zk_semi_mul_2.verify_as_subprotocol::<G1Element, ZpElement>(
                srs,
                trans_seq,
            );
    
            let check5 = zk_mul_scalar.verify_as_subprotocol::<ZpElement, ZpElement, ZpElement>(
                srs, 
                trans_seq,
            );
    
            let check6 = zk_schnorr.verify_as_subprotocol(trans_seq);

            return check1 && check2 && check3 && check4 && check5 && check6;

        } else {
            println!("!! Invalid transcript when verifying IpGt");
            return false;
        }
        
    }

    pub fn verify(&self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        if trans_seq.check_fiat_shamir() == false {
            println!("!! Fiat shamir check failed when verifying IpGt");
            return false;
        }

        self.verify_as_subprotocol(srs, trans_seq)
    }

}


#[cfg(test)]
mod tests {
    
    use super::*;
    use rand::Rng;

    #[test]
    fn test_zk_ipgt() {

        let srs = SRS::new(64);

        let n = 32 as usize;
        let a_vec = (0..n).map(|_| 
            rand::thread_rng().gen_range(0..2i64.pow(26))
        ).collect::<Vec<i64>>();

        let b_vec = (0..n).map(|_| 
            rand::thread_rng().gen_range(0..2i64.pow(26))
        ).collect::<Vec<i64>>();


        let a_com = srs.g_hat * dirac::inner_product(
            &a_vec, &srs.h_hat_vec
        );

        let b_com = srs.h_hat * dirac::inner_product(
            &b_vec, &srs.g_hat_vec
        );

        let c = dirac::inner_product(&a_vec, &b_vec);

        let blind_base = srs.blind_base;

        let c_tilde = ZpElement::rand();

        let gt_com = 
            c * srs.g_hat * srs.h_hat + a_com + b_com
            + c_tilde * blind_base;

        let gt_ip = ZkIpGt::new(
            gt_com, n
        );

        let mut trans_seq_prover = ZkTranSeqProver::new(&srs);

        gt_ip.prove::<i64>(
            &srs, 
            &mut trans_seq_prover, 
            &a_vec, 
            &b_vec,
            c_tilde,
        );

        let mut trans_seq = trans_seq_prover.publish_trans();

        println!(" * Proof of ZkIpGt generated {:?}", trans_seq.data.len());

        let result = gt_ip.verify(&srs, &mut trans_seq);

        assert_eq!(result, true);
        assert_eq!(trans_seq.data.len(), trans_seq.pointer);

        println!(" * Verification of ZkIpGt passed");

    }

    
}

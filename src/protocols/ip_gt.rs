//! Implementation of the Inner Product protocol in Gt
//!
//! Details of this protocol can be found in the DualMatrix paper 
//!
//! To prove that holding two vectors vec{a} and vec{b} such that
//! gt_com = <vec{a}, vec{b}> e(\hat{G}, \hat{H}) 
//!         + e(\hat{G}, <vec{a}, vec{H}>) + e(<vec{b}, \hat{G}>, \hat{H})

use crate::setup::SRS;

use crate::utils::curve::{ZpElement, GtElement, ConvertToZp};
use crate::utils::dirac;
use crate::utils::fiat_shamir::{TranElem, TranSeq};

use crate::protocols::pip::{PipG1, PipG2};


pub struct IpGt {
    pub gt_com: GtElement,
    pub length: usize,
}

impl IpGt {
    pub fn new(gt_com_value: GtElement, length_value: usize
    ) -> Self {
        Self {
            gt_com: gt_com_value,
            length: length_value,
        }
    }

    pub fn prove<T>(
        &self, srs: &SRS, 
        trans_seq: &mut TranSeq, 
        vec_a: &Vec<T>, 
        vec_b: &Vec<T>,
    ) where T: 'static + ConvertToZp {

        trans_seq.push(TranElem::Gt(self.gt_com));
        trans_seq.push(TranElem::Size(self.length));

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

            trans_seq.push(TranElem::Gt(l_tr));
            trans_seq.push(TranElem::Gt(r_tr));
            
            let x_j = trans_seq.gen_challenge();
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

        

        trans_seq.push(TranElem::Zp(a_reduce));
        trans_seq.push(TranElem::Zp(b_reduce));
        trans_seq.push(TranElem::G1(g_reduce));
        trans_seq.push(TranElem::G2(h_reduce));


        let pip_g1 = PipG1::new(
            g_reduce, &challenges_inv
        );

        let pip_g2 = PipG2::new(
            h_reduce, &challenges
        );


        pip_g1.prove(srs, trans_seq);
        pip_g2.prove(srs, trans_seq);

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
            println!("!! Invalid public input when verifying IpGt");
            return false;
        } 

        let n = self.length;
        let log_n = (n as u64).ilog2() as usize;

        trans_seq.pointer = pointer_old + 3 * log_n + 6;

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
            TranElem::Zp(a_reduce),
            TranElem::Zp(b_reduce),
            TranElem::G1(g_reduce),
            TranElem::G2(h_reduce),
        ) = (
            trans_seq.data[current_pointer],
            trans_seq.data[current_pointer+1],
            trans_seq.data[current_pointer+2],
            trans_seq.data[current_pointer+3],
        ) {
            let rhs = 
                srs.g_hat * a_reduce * h_reduce
                + g_reduce * b_reduce * srs.h_hat
                + a_reduce * b_reduce * srs.g_hat * srs.h_hat;
            if lhs == rhs {
                
                let pip_g1 = PipG1::new(
                    g_reduce, &challenges_inv
                );
        
                let pip_g2 = PipG2::new(
                    h_reduce, &challenges
                );

                let check_1 = pip_g1.verify_as_subprotocol(srs, trans_seq);
                let check_2 = pip_g2.verify_as_subprotocol(srs, trans_seq);
                return check_1 && check_2;

            }
        } else {
            println!("!! Invalid transcript when verifying IpGt");
            return false;
        }
    
        println!("!! Verification of IpGt failed");
        return false;
        
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
    fn test_ipgt() {

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

        let gt_com = c * srs.g_hat * srs.h_hat + a_com + b_com;

        let gt_ip = IpGt::new(
            gt_com, n
        );

        let mut trans_seq = TranSeq::new();

        gt_ip.prove::<i64>(&srs, &mut trans_seq, &a_vec, &b_vec);

        let result = gt_ip.verify(&srs, &mut trans_seq);

        assert_eq!(result, true);
        assert_eq!(trans_seq.data.len(), 3 * 5 + 6 + 4 + 4);

        println!(" * Verification of IpGt passed");

    }

    
}

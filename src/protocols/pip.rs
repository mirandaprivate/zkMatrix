///! Implementation of the Public Inner Product protocol
///!
///! Details of this protocol can be found in the DualMatrix paper 
/// 
use crate::setup::SRS;

use crate::utils::curve::{ZpElement, G1Element, G2Element};
use crate::utils::dirac;
use crate::utils::fiat_shamir::{TranElem, TranSeq};
use crate::utils::xi;


pub struct PIP_G1 {
    pub v: G1Element,
    pub challenges: Vec<ZpElement>,
}

impl PIP_G1 {
    pub fn new(v_value: G1Element, challenges_value: &Vec<ZpElement>
    ) -> Self {
        Self {
            v: v_value,
            challenges: challenges_value.clone(),
        }
    }

    pub fn prove(&mut self, srs: &SRS, trans_seq: &mut TranSeq) -> TranSeq {

        trans_seq.push(TranElem::G1(self.v));

        let s = trans_seq.gen_challenge();

        let v_prime_tr = xi::reduce_from_challenges(
            &self.challenges, &srs.h_hat_prime_vec
        );

        let log_n = self.challenges.len() as usize;
        let n = 2usize.pow(log_n as u32);

        let xi = xi::xi_from_challenges(&self.challenges);
        let psi = xi::psi_from_xi(&xi, s);

        let mut g_hat_1: Vec<G1Element> = vec![srs.g_hat,];
        g_hat_1.append(
            &mut srs.g_hat_vec[..(n-1)].to_vec()
        );

        let w_tr = dirac::inner_product(&psi, &g_hat_1);

        trans_seq.push(TranElem::G2(v_prime_tr));
        trans_seq.push(TranElem::G1(w_tr));

        trans_seq.clone()

    }

    pub fn verify_as_subprotocol(&mut self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        let pointer_old = trans_seq.pointer;
        trans_seq.pointer = pointer_old + 4;

        if TranElem::G1(self.v) != trans_seq.data[pointer_old] {
            return false;
        }

        if let (
            TranElem::Zp(s),
            TranElem::G2(v_prime_tr),
            TranElem::G1(w_tr),
         ) = (
            trans_seq.data[pointer_old + 1],
            trans_seq.data[pointer_old + 2],
            trans_seq.data[pointer_old + 3],
        ) {
            
            let s_hat_h_hat = srs.h_hat_prime_vec[0];
            
            let phi_s = xi::phi_s(
                s, &self.challenges, 1, 1);

            if ( w_tr * (s * srs.h_hat - s_hat_h_hat)
                == (phi_s * srs.g_hat - self.v) * (srs.h_hat) 
             ) && (
                self.v * srs.h_hat == srs.g_hat * v_prime_tr
             )
            {
                return true;
            } 

        } 

        return false;
    }

    pub fn verify(&mut self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        trans_seq.pointer = 0;

        if trans_seq.check_fiat_shamir() == false {
            return false;
        }

        self.verify_as_subprotocol(srs, trans_seq)
    }

}

#[cfg(test)]
mod tests {
    
    use super::*;

    #[test]
    fn test_pip() {
        
    }

    
}
//! Implementation of the Public Inner Product protocol
//!
//! Details of this protocol can be found in the DualMatrix paper 
// 
use crate::setup::SRS;

use crate::utils::curve::{ZpElement, G1Element, G2Element};
use crate::utils::dirac;
use crate::utils::fiat_shamir::{TranElem, TranSeq};
use crate::utils::xi;


pub struct PipG1 {
    pub v: G1Element,
    pub challenges: Vec<ZpElement>,
}

pub struct PipG2 {
    pub v: G2Element,
    pub challenges: Vec<ZpElement>,
}

impl PipG1 {
    pub fn new(v_value: G1Element, challenges_value: &Vec<ZpElement>
    ) -> Self {
        Self {
            v: v_value,
            challenges: challenges_value.clone(),
        }
    }

    pub fn prove(&self, srs: &SRS, trans_seq: &mut TranSeq)  {

        // println!("Proving the PIPs...");

        trans_seq.push(TranElem::G1(self.v));

        let s = trans_seq.gen_challenge();

        let v_prime_tr = xi::reduce_from_challenges(
            &self.challenges, &srs.h_hat_prime_vec
        );

        let log_n = self.challenges.len() as usize;
        let n = 2usize.pow(log_n as u32);
        

        let xi = xi::xi_from_challenges(&self.challenges);
        let psi = xi::psi_from_xi(&xi, s);

        // println!("srs_g_hat: {:?}", srs.g_hat_vec.len());

        let g_hat_1: Vec<G1Element> = vec![srs.g_hat.clone(),]
            .iter().chain(srs.g_hat_vec.iter().take(n-1))
            .cloned().collect::<Vec<G1Element>>();
      
        // println!("g_hat_1: {:?}", g_hat_1);

        let w_tr = dirac::inner_product(&g_hat_1, &psi);

        trans_seq.push(TranElem::G2(v_prime_tr));
        trans_seq.push(TranElem::G1(w_tr));


    }

    pub fn verify_as_subprotocol(
        &self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        let pointer_old = trans_seq.pointer;
        trans_seq.pointer = pointer_old + 4;

        if TranElem::G1(self.v) != trans_seq.data[pointer_old] {
            return false;
        }

        if let (
            TranElem::Coin(s),
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

            println!("!! Pip equation check failed when verifying Pip");

        } 

        println!("!! Type check for transcript elements failed when verifying Pip");
        return false;
    }

    pub fn verify(&self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        if trans_seq.check_fiat_shamir() == false {
            println!("!! Fiat shamir check failed when verifying Pip");
            return false;
        }

        self.verify_as_subprotocol(srs, trans_seq)
    }

}




impl PipG2 {
    pub fn new(v_value: G2Element, challenges_value: &Vec<ZpElement>
    ) -> Self {
        Self {
            v: v_value,
            challenges: challenges_value.clone(),
        }
    }

    pub fn prove(&self, srs: &SRS, trans_seq: &mut TranSeq)  {

        trans_seq.push(TranElem::G2(self.v));

        let s = trans_seq.gen_challenge();

        let s_pow_q = s.pow(srs.q as u64);

        let v_prime_tr = xi::reduce_from_challenges(
            &self.challenges, &srs.g_hat_prime_vec
        );

        let log_n = self.challenges.len() as usize;
        let n = 2usize.pow(log_n as u32);
        

        let xi = xi::xi_from_challenges(&self.challenges);
        let psi = xi::psi_from_xi(&xi, s_pow_q);

        // println!("srs_g_hat: {:?}", srs.g_hat_vec.len());

        let h_hat_1: Vec<G2Element> = vec![srs.h_hat.clone(),]
            .iter().chain(srs.h_hat_vec.iter().take(n-1))
            .cloned().collect::<Vec<G2Element>>();
      
        // println!("g_hat_1: {:?}", g_hat_1);

        let w_tr = dirac::inner_product(&h_hat_1, &psi);

        trans_seq.push(TranElem::G1(v_prime_tr));
        trans_seq.push(TranElem::G2(w_tr));


    }

    pub fn verify_as_subprotocol(
        &self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        let pointer_old = trans_seq.pointer;
        trans_seq.pointer = pointer_old + 4;

        if TranElem::G2(self.v) != trans_seq.data[pointer_old] {
            return false;
        }

        if let (
            TranElem::Coin(s),
            TranElem::G1(v_prime_tr),
            TranElem::G2(w_tr),
         ) = (
            trans_seq.data[pointer_old + 1],
            trans_seq.data[pointer_old + 2],
            trans_seq.data[pointer_old + 3],
        ) {
            
            let s_hat_pow_q_g_hat = srs.g_hat_prime_vec[0];
            
            let q = srs.q as usize;

            let phi_s = xi::phi_s(
                s, &self.challenges, q, q);

            let s_pow_q = s.pow(q as u64);

            if ( (s_pow_q * srs.g_hat - s_hat_pow_q_g_hat) * w_tr 
                == (srs.g_hat) * (phi_s * srs.h_hat - self.v)  
             ) && (
                srs.g_hat * self.v == v_prime_tr * srs.h_hat
             )
            {
                return true;
            } 

            println!("!! Pip equation check failed when verifying Pip");

        } 

        println!("!! Type check for transcript elements failed when verifying Pip");
        return false;
    }

    pub fn verify(&self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {


        if trans_seq.check_fiat_shamir() == false {
            println!("!! Fiat shamir check failed when verifying Pip");
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

        let srs = SRS::new(64);

        let challenges = (0..5).map(|_| 
            ZpElement::rand()
        ).collect::<Vec<ZpElement>>();

        let v1_test = xi::reduce_from_challenges(
            &challenges, &srs.g_hat_vec
        );

        let v2_test = xi::reduce_from_challenges(
            &challenges, &srs.h_hat_vec
        );
        
        let pip_g1_protocol = PipG1::new(
            v1_test, &challenges
        );

        let pip_g2_protocol = PipG2::new(
            v2_test, &challenges
        );


        let mut trans_seq = TranSeq::new();

        pip_g1_protocol.prove(&srs, &mut trans_seq);

        pip_g2_protocol.prove(&srs, &mut trans_seq);

  
        let result1 = pip_g1_protocol.verify(
            &srs, &mut trans_seq
        );

        assert_eq!(result1, true);

        println!(" * Verifications of PipG1 passed");

        let result2 = pip_g2_protocol.verify(&srs, &mut trans_seq);

        assert_eq!(result2, true);

        println!(" * Verifications of PipG1 passed");

    }

    
}
///! Implementation of the Inner Product protocol
///!
///! Details of this protocol can be found in the zkMatrix paper 
///!
///! To prove that holding a vector vec{a} such that
///! <vec{a}, vec{G}> = vec_com,

use crate::setup::SRS;

use crate::utils::curve::{ZpElement, G1Element, ConvertToZp};
use crate::utils::dirac;
use crate::utils::fiat_shamir::{TranElem, TranSeq};

use crate::protocols::pip::PipG1;

pub struct VecComG1 {
    pub vec_com: G1Element,
    pub length: usize,
}

impl VecComG1 {
    pub fn new(vec_com_value: G1Element, length_value: usize
    ) -> Self {
        Self {
            vec_com: vec_com_value,
            length: length_value,
        }
    }

    pub fn prove<T>(
        &self, srs: &SRS, trans_seq: &mut TranSeq, vec_witness: &Vec<T>
    ) where T: 'static + ConvertToZp {

        trans_seq.push(TranElem::G1(self.vec_com));
        trans_seq.push(TranElem::Size(self.length));

        if (self.length & (self.length - 1)) != 0 {
            panic!("Length is not a power of 2 when proving VecCom");
        }

        let n = self.length;
        let log_n = (n as u64).ilog2() as usize;

        let mut vec_current = dirac::vec_convert_to_zp_vec(
            vec_witness);
        let mut g_vec_current = srs.g_hat_vec[0..n].to_vec();
        
        let mut challenges: Vec<ZpElement> = Vec::new();

        for j in 0..log_n {
            let current_len = n / 2usize.pow(j as u32);
            let vec_left = 
                vec_current[0..current_len/2].to_vec();
            let vec_right = 
                vec_current[current_len/2..current_len].to_vec();
            let g_left = 
                g_vec_current[0..current_len/2].to_vec();
            let g_right = 
                g_vec_current[current_len/2..current_len].to_vec();

            let l_tr = dirac::inner_product(
                &vec_left, &g_right
            );
            let r_tr = dirac::inner_product(
                &vec_right, &g_left
            );

            trans_seq.push(TranElem::G1(l_tr));
            trans_seq.push(TranElem::G1(r_tr));
            
            let x_j = trans_seq.gen_challenge();
            let x_j_inv = x_j.inv();

            challenges.push(x_j);

            vec_current = dirac::vec_addition(
                &vec_left,
                &dirac::vec_scalar_mul(
                    &vec_right, x_j_inv),
            );

            g_vec_current = dirac::vec_addition(
                &g_left,
                &dirac::vec_scalar_mul(
                    &g_right, x_j),
            );

        }

        let a_reduce = vec_current[0];
        let g_reduce = g_vec_current[0];

        trans_seq.push(TranElem::Zp(a_reduce));
        trans_seq.push(TranElem::G1(g_reduce));

        let pip_g1 = PipG1::new(
            g_reduce, &challenges
        );

        pip_g1.prove(srs, trans_seq)

    }

    pub fn verify_as_subprotocol(
        &self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        let pointer_old = trans_seq.pointer;
        
        if (
            TranElem::G1(self.vec_com),
            TranElem::Size(self.length),
        ) != (
            trans_seq.data[pointer_old],
            trans_seq.data[pointer_old + 1],
        ) {
            println!("!! Invalid public input when verifying VecCom");
            return false;
        } 

        let n = self.length;
        let log_n = (n as u64).ilog2() as usize;

        trans_seq.pointer = pointer_old + 3 * log_n + 4;

        let mut current_pointer = pointer_old + 2;
        let mut lhs: G1Element = self.vec_com.clone();
        
        let mut challenges: Vec<ZpElement> = Vec::new();

        for _ in 0..log_n {

            if let (
                TranElem::G1(l_tr),
                TranElem::G1(r_tr),
                TranElem::Coin(x_j),
            ) = (
                trans_seq.data[current_pointer],
                trans_seq.data[current_pointer + 1],
                trans_seq.data[current_pointer + 2],
            ) {
                let x_j_inv = x_j.inv();
                lhs = lhs + l_tr * x_j + r_tr * x_j_inv;
                challenges.push(x_j);
            } else {
                println!("!! Invalid transcript when verifying VecCom");
                return false;
            }

            current_pointer += 3;
        }

        if let (
            TranElem::Zp(a_reduce),
            TranElem::G1(g_reduce),
        ) = (
            trans_seq.data[current_pointer],
            trans_seq.data[current_pointer+1],
        ) {
            let rhs = g_reduce * a_reduce;
            if lhs == rhs {
                let pip_g1 = PipG1::new(
                    g_reduce, &challenges
                );

                return pip_g1.verify(srs, trans_seq);
            }
        } else {
            println!("!! Invalid transcript when verifying VecCom");
            return false;
        }
    
        println!("!! Verification of VecCom failed");
        return false;
        
    }

    pub fn verify(&self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        if trans_seq.check_fiat_shamir() == false {
            println!("!! Fiat shamir check failed when verifying VecCom");
            return false;
        }

        self.verify_as_subprotocol(srs, trans_seq)
    }

}


#[cfg(test)]
mod tests {
    
    use super::*;

    #[test]
    fn test_vec_com() {

        let srs = SRS::new(64);

        let n = 32 as usize;
        let a_vec = (0..n).map(|_| 
            ZpElement::rand()
        ).collect::<Vec<ZpElement>>();

        let a_com = dirac::inner_product(
            &a_vec, &srs.g_hat_vec
        );

        let vec_com_ip = VecComG1::new(
            a_com, n
        );


        let mut trans_seq = TranSeq::new();

        vec_com_ip.prove::<ZpElement>(&srs, &mut trans_seq, &a_vec);

        let result = vec_com_ip.verify(&srs, &mut trans_seq);

        assert_eq!(result, true);
        assert_eq!(trans_seq.data.len(), 3 * 5 + 4 + 4);

        println!(" * Verification of VecCom passed");

    }

    
}


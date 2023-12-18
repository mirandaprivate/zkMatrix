//! Define the SRS used in the library.
//! 
//! 
use rand::Rng;
use serde::{Serialize, Deserialize};

use crate::utils::curve::{G1Element, G2Element, ZpElement};

use crate::utils::to_file::FileIO;

#[derive(Debug, PartialEq, Eq, Clone, Serialize, Deserialize)]
pub struct SRS {
    q: usize,
    g_hat: G1Element,
    h_hat: G2Element,
    g_hat_minus: G1Element,
    h_hat_minus: G2Element,
    g_hat_vec: Vec<G1Element>,
    h_hat_vec: Vec<G2Element>,
    g_hat_prime_vec: Vec<G1Element>,
    h_hat_prime_vec: Vec<G2Element>,
}

impl SRS {
    pub fn new(q: usize) -> Self {
        
        let s_hat = ZpElement::from(
            rand::thread_rng().gen::<u64>());
        // let nu = ZpElement::from(
        //     rand::thread_rng().gen::<u64>());

        let g_hat = G1Element::generator();
        let h_hat = G2Element::generator();

        let s_vec: Vec<ZpElement> = std::iter::successors(
            Some(s_hat), 
            |&x| Some(x * s_hat)
        ).take(q-1).collect();

        let s_hat_pow_q = s_hat * (*s_vec.last().unwrap());

        let s_prime_vec: Vec<ZpElement> = std::iter::successors(
            Some(s_hat_pow_q), 
            |&x| Some(x * s_hat_pow_q)
        ).take(q-1).collect();

        let s_hat_pow_q_squre = 
            s_hat_pow_q * (*s_prime_vec.last().unwrap());
        
        let g_hat_minus = s_hat_pow_q_squre * g_hat;
        let h_hat_minus = s_hat_pow_q_squre * h_hat;

        let g_hat_vec = s_vec.iter()
            .map(|&x| x * g_hat)
            .collect();

        let g_hat_prime_vec = s_prime_vec.iter()
            .map(|&x| x * g_hat)
            .collect();

        let h_hat_vec = s_prime_vec.iter()
            .map(|&x| x * h_hat)
            .collect();

        let h_hat_prime_vec = s_vec.iter()
           .map(|&x| x * h_hat)
            .collect();

        Self {
            q,
            g_hat,
            h_hat,
            g_hat_minus,
            h_hat_minus,
            g_hat_vec,
            h_hat_vec,
            g_hat_prime_vec,
            h_hat_prime_vec,
         }

    }

}


impl FileIO for SRS {}

#[cfg(test)]
mod tests {

    use super::*;

    use crate::{config::Q_TEST, utils::to_file::FileIO};

    #[test]
    fn test_srs() {
        let srs = SRS::new(Q_TEST);
        srs.to_file(String::from("srs"), true).unwrap();
        let srs_read = FileIO::from_file(String::from("srs"), true).unwrap();
        println!("srs: {:?}", srs);
        assert_eq!(srs, srs_read);
    }
}
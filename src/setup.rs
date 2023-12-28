//! Generate the SRS required by the protocols.
//! 
//! Here, s_hat is the toxic waste which should be discarded securely
//! in real applications.
//! 
use serde::{Serialize, Deserialize};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

use crate::utils::curve::{G1Element, G2Element, ZpElement, GtElement};

use crate::utils::to_file::FileIO;

use crate::config::NUM_THREADS;

#[derive(Debug, PartialEq, Eq, Clone, Serialize, Deserialize)]
pub struct SRS {
    pub q: usize,
    pub g_hat: G1Element,
    pub h_hat: G2Element,
    pub blind_base: GtElement,
    pub g_hat_vec: Vec<G1Element>,
    pub h_hat_vec: Vec<G2Element>,
    pub g_hat_prime_vec: Vec<G1Element>,
    pub h_hat_prime_vec: Vec<G2Element>,
}

impl SRS {
    pub fn new(q: usize) -> Self {
        
        let s_hat = ZpElement::rand();
        
        let g_hat = G1Element::generator().clone();
        let h_hat = G2Element::generator().clone();

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
        
        let blind_base = s_hat_pow_q_squre * g_hat * h_hat;

        let pool = ThreadPoolBuilder::new()
            .num_threads(NUM_THREADS)
            .build()
            .unwrap();

        let g_hat_vec = 
            pool.install(|| {
                s_vec.par_iter()
                .map(|&x| x * g_hat)
                .collect()
                }
            );

        let g_hat_prime_vec = 
            pool.install(|| {
                s_prime_vec.par_iter()
                .map(|&x| x * g_hat)
                .collect()
                }
            );

        let h_hat_vec = 
            pool.install(|| {
                s_prime_vec.par_iter()
                .map(|&x| x * h_hat)
                .collect()
                }
            );

        let h_hat_prime_vec = 
            pool.install(|| {
                s_vec.iter()
            .map(|&x| x * h_hat)
                .collect()
                }
            );

        let srs = Self {
            q,
            g_hat,
            h_hat,
            blind_base,
            g_hat_vec,
            h_hat_vec,
            g_hat_prime_vec,
            h_hat_prime_vec,
        };

        srs.to_file("srs.srs".to_string(), true).unwrap();

        srs

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
        let srs_read = FileIO::from_file(
            String::from("srs.srs"), true
        ).unwrap();
        println!("srs: {:?}", srs);
        assert_eq!(srs, srs_read);
    }
}
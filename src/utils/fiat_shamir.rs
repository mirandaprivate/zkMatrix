//! Utility functions for the adaptive Fiat-Shamir transformation
//! 
use serde::{Serialize, Deserialize};
use sha2::{Sha256, Digest};

use crate::utils::curve::{ZpElement, G1Element, G2Element, GtElement, Zero};
use crate::utils::to_file::FileIO;


#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TranElem{
    G1(G1Element),
    G2(G2Element),
    Gt(GtElement),
    Zp(ZpElement),
    Size(usize),
    Coin(ZpElement),
}

#[derive(Debug, Clone)]
pub struct TranSeq {
    pub hasher: Sha256,
    pub data: Vec<TranElem>,
    pub pointer: usize,

}


impl TranSeq{
 
    pub fn new() -> Self {
        let mut hasher_val = Sha256::new();
        hasher_val.update(b"0");

        Self { 
            hasher: hasher_val,
            data: Vec::new(),
            pointer: 0,
        }
    }

    pub fn reset_hash(&mut self) {
        let mut hasher_val = Sha256::new();
        hasher_val.update(b"0");

        self.hasher = hasher_val.clone();
        self.pointer = 0;
    }

    pub fn push(&mut self, elem: TranElem) {
        let bytes = 
            bincode::serialize(&elem).unwrap();
        self.hasher.update(bytes);
        self.data.push(elem);
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn gen_challenge(&mut self) -> ZpElement {
        // println!("Gen challenge");
        // println!("{:?}", &self.hasher.clone().finalize());
        let hash_value: [u8; 32] = self.hasher.clone().finalize().into();
        // println!("Hash value: {:?}", hash_value);
        let mut challenge = ZpElement::from_u8(
            &hash_value);
        // println!("Challenge: {:?}", challenge);

        while challenge == ZpElement::zero() {
            self.hasher.update(b"0");
            let hash_value: [u8; 32] = self.hasher.clone().finalize().into();
            challenge = ZpElement::from_u8(
                &hash_value);
        }

        self.data.push(TranElem::Coin(challenge));

        challenge
    }

    pub fn check_fiat_shamir(&self) -> bool {
        let n = self.data.len();
        let mut new_hasher = Sha256::new();
        new_hasher.update(b"0");

        for i in 0..n {
            let current = self.data[i].clone();
    
            match current{
                TranElem::G1(_) | TranElem::G2(_) | TranElem::Gt(_) | 
                    TranElem::Zp(_) | TranElem::Size(_)
                => {
                    let bytes 
                        = bincode::serialize(&current).unwrap();
                    new_hasher.update(bytes);
                },
    
                TranElem::Coin(coin_value) => {
                    let hash_value: [u8; 32] = new_hasher.clone().finalize().into();
                    let mut challenge = ZpElement::from_u8(
                        &hash_value);
                    while challenge == ZpElement::zero() {
                        new_hasher.update(b"0");
                        let hash_value: [u8; 32] = new_hasher.clone().finalize().into();
                        challenge = ZpElement::from_u8(
                            &hash_value);
                    }
                    if challenge != coin_value{
                        return false;
                    }
                },

            }
        }
        return true;
    }

    pub fn save_to_file(&self, file_name: String) -> std::io::Result<()> {
        self.data.to_file(format!("{}.tr", file_name), true)
    }

    pub fn read_from_file(file_name: String) -> Self {
        
        let data_read: Vec<TranElem> = 
            FileIO::from_file(format!("{}.tr", file_name), true)
                .unwrap();

        let mut hasher_val = Sha256::new();
        hasher_val.update(b"0");
        
        Self { 
            hasher: hasher_val,
            data: data_read,
            pointer: 0,
        }
    }

}


#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_fiat_shamir(){
        
        let mut trans = TranSeq::new();

        trans.push(TranElem::Zp(ZpElement::from(1 as u64)) );

        trans.gen_challenge();

        trans.push(TranElem::Zp(ZpElement::from(2 as u64)) );

        trans.gen_challenge();

        trans.save_to_file(String::from("tr_test")).unwrap();

        let trans_read = TranSeq::read_from_file(
            String::from("tr_test")
        );

        assert_eq!(trans_read.data.len(), 4);
        assert_eq!(trans_read.check_fiat_shamir(), true);


    }
}
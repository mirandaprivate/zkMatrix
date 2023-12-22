use serde::{Serialize, Deserialize};
use sha2::{Sha256, Digest};

use crate::utils::curve::{ZpElement, G1Element, G2Element, GtElement, Zero};


#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TranElem{
    G1(G1Element),
    G2(G2Element),
    Gt(GtElement),
    Zp(ZpElement),
    Coin(ZpElement),
}


pub struct TranSeq {
    pub hasher: Sha256,
    pub data: Vec<TranElem>,

}


impl TranSeq{
 
    pub fn new() -> Self {
        Self { 
            hasher: Sha256::new(),
            data: Vec::new(),
        }
    }

    pub fn push(&mut self, elem: TranElem) {
        self.data.push(elem);
    }

    pub fn pop_begining(&mut self) -> TranElem {
        self.data.remove(0)
    }
    pub fn len(&self) -> usize {
        self.data.len()
    }

}

pub fn adaptive_fiat_shamir(trans_seq: &mut TranSeq) -> TranSeq
{
    let n = trans_seq.len();
    let mut output_seq = TranSeq::new();
    let mut hasher = Sha256::new();

    for _ in 0..n {
        let current = trans_seq.pop_begining();

        match current{
            TranElem::G1(_) | TranElem::G2(_) | TranElem::Gt(_) | TranElem::Zp(_) 
            => {
                let bytes 
                    = bincode::serialize(&current).unwrap();
                hasher.update(bytes);
                output_seq.push(current);
            },

            TranElem::Coin(_) => {
                let mut challenge = ZpElement::from_bytes(
                    &hasher.clone().finalize().into());
                while challenge == ZpElement::zero() {
                    hasher.update(b"0");
                    challenge = ZpElement::from_bytes(
                        &hasher.clone().finalize().into());
                }
                output_seq.push(TranElem::Coin(challenge));
            },
        }
    }
    output_seq
}
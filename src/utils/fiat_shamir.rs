use serde::ser::Serializer;
use serde::de::Deserializer;
use serde::{Serialize, Deserialize};
use sha2::{Sha256, Digest};

use crate::utils::curve::{ZpElement, G1Element, G2Element, GtElement};


#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TranElem{
    G1(G1Element),
    G2(G2Element),
    Gt(GtElement),
    Zp(ZpElement),
    Coin(ZpElement),
}


pub struct TranSeq {
    pub data: Vec<TranElem>,

}


impl TranSeq{
 
    pub fn new() -> Self {
        Self { data: Vec::new::<TranElem>() }
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


pub fn adaptive_fiat_shamir(trans_seq: &TranSeq) -> TranSeq
{
    let mut hasher = Sha256::new();

    let n = trans_seq.len();
    let output_seq = TranSeq::new();

    for _ in 0..n {
        let current = trans_seq.pop_begining();

        match current{
            TranElem::G1(_) | TranElem::G2(_) | TranElem::Gt(_) | TranElement::Zp(_) 
            => {
                let mut bytes = Vec::new();
                current.serialize(&mut Serializer::new(&mut bytes)).unwrap();
                hasher.update(bytes);
                output_seq.push(current);
            },

            TranElem::Coin(_) => {
                let mut challenge = ZpElement::from_bytes(
                    &hasher.clone().finalize());
                while challenge == ZpElement::zero() {
                    hasher.update(b"0");
                    challenge = ZpElement::from_bytes(
                        &hasher.clone().finalize());
                }
                output_seq.push(TranElem::Coin(challenge));
            },
        }
    }
    output_seq
}
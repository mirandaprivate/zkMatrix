//! Define the transcript sequence for zero-knowledge protocols
//! 
//! In essence, in the prover algorithm of a zero-knowledge protocol,
//! the prover can hide each transcript element by accompanying it with
//! a blinding factor.
//!
use crate::setup::SRS;

use crate::utils::curve::{ZpElement, G1Element, G2Element, GtElement};
use crate::utils::fiat_shamir::{TranElem, TranSeq};


/// Turn each transcript element into a hidden commitment
/// and recording the blinding factors
/// 
pub struct  ZkTranSeqProver {
    pub trans_seq: TranSeq,
    pub blind_seq: Vec<ZpElement>,
    pub blind_base: GtElement,
    pub g_hat: G1Element,
    pub h_hat: G2Element,
} 

impl ZkTranSeqProver {

    pub fn new(srs: &SRS) -> Self {
        let trans_seq = TranSeq::new();
        let blind_seq: Vec<ZpElement> = Vec::new();
        Self { 
            trans_seq: trans_seq, 
            blind_seq: blind_seq, 
            blind_base: srs.blind_base,
            g_hat: srs.g_hat,
            h_hat: srs.h_hat,
        }
    }

    pub fn gen_challenge(&mut self) -> ZpElement {
        self.trans_seq.gen_challenge()
    }

    pub fn get_mut_trans_seq(&mut self) -> &mut TranSeq {
        &mut self.trans_seq
    }

    pub fn push_without_blinding(&mut self, tr_elem: TranElem) {
        self.trans_seq.push(tr_elem);
    }

    pub fn push_with_blinding  (
        &mut self, tr_elem: TranElem, blinding_factor: ZpElement
    ) -> GtElement {

        let tr_gt: GtElement;

        match tr_elem {
            
            TranElem::Zp(tr_value)
                => { tr_gt = tr_value * self.g_hat * self.h_hat; }
            TranElem::G1(tr_value)
                => { tr_gt = tr_value * self.h_hat; }
            TranElem::G2(tr_value)
                => { tr_gt = tr_value * self.g_hat; }
            TranElem::Gt(tr_value)
                => { tr_gt = tr_value; }
            _ => { panic!("Unable to blind"); }
        }

        let tr_blind = 
            tr_gt + blinding_factor * self.blind_base;
    
        self.trans_seq.push(TranElem::Gt(tr_blind));
        self.blind_seq.push(blinding_factor);
        tr_blind
    }

    pub fn push_gen_blinding(&mut self, tr_elem: TranElem
    ) -> (GtElement, ZpElement) {
        
        let blinding_factor = ZpElement::rand();

        let gt_val = self.push_with_blinding(
            tr_elem, blinding_factor
        );

        (gt_val, blinding_factor)
    }
    
    pub fn publish_trans(&mut self) -> TranSeq {
        self.trans_seq.clone()
    }

}
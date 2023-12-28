//! Zero-knowledge protocols for scalar arithmetic
//! 
use crate::setup::SRS;

use crate::utils::curve::{
    G1Element, G2Element, ZpElement, GtElement,
    Mul, Add, Zero
};
use crate::utils::fiat_shamir::{TranElem, TranSeq};

/// Prove know a scalar a such that a * base = commit
pub struct ZkSchnorr
{
    pub commit: GtElement,
    pub base: GtElement,
}

/// Prove holding scalar c, a, b and \tilde{c}, \tilde{a}, \tilde{b} such that:
/// c_com = c.toGt() + \tilde{c} * blind_base
/// a_com = a.toGt() + \tilde{a} * blind_base
/// b_com = b.toGt() + \tilde{b} * blind_base
/// c = a * b
/// 
pub struct ZkMulScalar
{
    pub c_com: GtElement,
    pub a_com: GtElement,
    pub b_com: GtElement,
}

/// Prove holding scalar c, a and \tilde{c}, \tilde{a} such that:
/// c_com = c.toGt() + \tilde{c} * blind_base
/// a_com = a.toGt() + \tilde{a} * blind_base
/// c = a * b
/// Here, b is a public scalar
/// 
pub struct ZkSemiMulScalar<V> 
where
    V: 'static + Copy + Clone + Add  + ToTranElem,
{
    pub c_com: GtElement,
    pub a_com: GtElement,
    pub b: V,
}

impl ZkSchnorr {
    
    pub fn new(
        commit_value: GtElement, 
        base_value: GtElement,
    ) -> Self {
        Self {
            commit: commit_value,
            base: base_value,
        }
    }

    pub fn prove(&self, trans_seq: &mut TranSeq, witness: ZpElement) {
        
        trans_seq.push(TranElem::Gt(self.commit));
        trans_seq.push(TranElem::Gt(self.base));

        let r = ZpElement::rand();
        let r_com = self.base * r;

        trans_seq.push(TranElem::Gt(r_com));

        let challenge = trans_seq.gen_challenge();

        let z = r + challenge * witness;

        trans_seq.push(TranElem::Zp(z));
    }

    pub fn verify_as_subprotocol(&self, trans_seq: &mut TranSeq
    ) -> bool {

        let pointer_old = trans_seq.pointer;

        trans_seq.pointer = pointer_old + 5;

        if (
            TranElem::Gt(self.commit),
            TranElem::Gt(self.base),
        ) != (
            trans_seq.data[pointer_old],
            trans_seq.data[pointer_old + 1],
        ) {
            println!("!! Invalid public input when verifying ZkSchnorr");
            return false;
        } 

        
        if let (
            TranElem::Gt(r_com),
            TranElem::Coin(challenge),
            TranElem::Zp(z),
        ) = (
            trans_seq.data[pointer_old + 2],
            trans_seq.data[pointer_old + 3],
            trans_seq.data[pointer_old + 4],
        ) {

           if r_com + challenge * self.commit == self.base * z {
               return true;
           } else {
               println!("!! ZkSchnorr equation check failed when verifying ZkSchnorr");
           } 
        } else {
            println!("!! Type check for transcript elements failed when verifying ZkSchnorr");
        }

        return false;
        
    } 

}

pub trait ToTranElem {
    fn to_gt(&self, srs: &SRS) -> GtElement;
    fn rand_general() -> Self;
    fn to_tran_elem(&self) -> TranElem;
    fn from_tran_elem(elem: TranElem) -> Self;
}

impl ToTranElem for ZpElement {
    fn to_gt(&self, srs: &SRS) -> GtElement {
        *self * srs.g_hat * srs.h_hat
    }

    fn rand_general() -> Self {
        ZpElement::rand()
    }

    fn to_tran_elem(&self) -> TranElem {
        TranElem::Zp(*self)
    }

    fn from_tran_elem(elem: TranElem) -> Self {
        if let TranElem::Zp(value) = elem {
            value
        } else {
            println!("Type mismatch when converting TranElem to ZpElement");
            ZpElement::zero()
        }
    }
}

impl ToTranElem for G1Element {
    fn to_gt(&self, srs: &SRS) -> GtElement {
        *self * srs.h_hat
    }

    fn rand_general() -> Self {
        ZpElement::rand() * G1Element::from(1 as u64)
    }

    fn to_tran_elem(&self) -> TranElem {
        TranElem::G1(*self)
    }

    fn from_tran_elem(elem: TranElem) -> Self {
        if let TranElem::G1(value) = elem {
            value
        } else {
            println!("Type mismatch when converting TranElem to G1Element");
            G1Element::zero()
        }
    }
}

impl ToTranElem for G2Element {
    fn to_gt(&self, srs: &SRS) -> GtElement {
        *self * srs.g_hat
    }

    fn rand_general() -> Self {
        ZpElement::rand() * G2Element::from(1 as u64)
    }

    fn to_tran_elem(&self) -> TranElem {
        TranElem::G2(*self)
    }


    fn from_tran_elem(elem: TranElem) -> Self {
        if let TranElem::G2(value) = elem {
            value
        } else {
            println!("Type mismatch when converting TranElem to G2Element");
            G2Element::zero()
        }
    }
}

impl ToTranElem for GtElement{
    fn to_gt(&self, _srs: &SRS) -> GtElement {
        *self
    }

    fn rand_general() -> Self {
        ZpElement::rand() * GtElement::from(1 as u64)
    }

    fn to_tran_elem(&self) -> TranElem {
        TranElem::Gt(*self)
    
    }


    fn from_tran_elem(elem: TranElem) -> Self {
        if let TranElem::Gt(value) = elem {
            value
        } else {
            println!("Type mismatch when converting TranElem to GtElement");
            GtElement::zero()
        }
    }
}

/// This is the simplified case of Algorithm 5 in the zkMatrix paper
impl ZkMulScalar
{
    
    pub fn new(
        c_com_value: GtElement, 
        a_com_value: GtElement,
        b_com_value: GtElement,
    ) -> Self {
        Self {
            c_com: c_com_value,
            a_com: a_com_value,
            b_com: b_com_value,
        }
    }

    pub fn prove<T, U, V>(
        &self, 
        srs: &SRS, 
        trans_seq: &mut TranSeq, 
        a: U, b: V,
        c_tilde: ZpElement, a_tilde: ZpElement, b_tilde: ZpElement,
    ) where
        T: 'static + Copy + Clone + ToTranElem 
            + Add<T, Output = T> + Mul<ZpElement, Output = T>,
        U: 'static + Copy + Clone + ToTranElem + Mul<V, Output = T> 
            + Add<U, Output=U> + Mul<ZpElement, Output = U>,
        V: 'static + Copy + Clone + ToTranElem 
            + Add<V, Output=V> + Mul<ZpElement, Output = V>,
     {
        
        trans_seq.push(TranElem::Gt(self.c_com));
        trans_seq.push(TranElem::Gt(self.a_com));
        trans_seq.push(TranElem::Gt(self.b_com));

        let alpha: U = U::rand_general();
        let beta: V = V::rand_general();

        let alpha_tilde = ZpElement::rand();
        let beta_tilde = ZpElement::rand();
        let gamma_tilde = ZpElement::rand();
        let delta_tilde = ZpElement::rand();

        let blind_base = srs.blind_base;
        let alpha_com 
            = alpha.to_gt(srs) + alpha_tilde * blind_base;
        let beta_com
            = beta.to_gt(srs) + beta_tilde * blind_base;
        let gamma_com =
            (alpha * beta).to_gt(srs) + gamma_tilde * blind_base;
        let delta_com =
            (alpha * b + a * beta).to_gt(srs) + delta_tilde * blind_base;
         
        trans_seq.push(TranElem::Gt(alpha_com));
        trans_seq.push(TranElem::Gt(beta_com));
        trans_seq.push(TranElem::Gt(gamma_com));
        trans_seq.push(TranElem::Gt(delta_com));

        let x = trans_seq.gen_challenge();

        let z_a = 
            - alpha_tilde - x * a_tilde;
        let z_b = 
            - beta_tilde - x * b_tilde;
        let z_c = 
            - gamma_tilde - x * delta_tilde - x * x * c_tilde;

        trans_seq.push(TranElem::Zp(z_a));
        trans_seq.push(TranElem::Zp(z_b));
        trans_seq.push(TranElem::Zp(z_c));

        let a_blind: U = 
            (alpha + a * x) * x;
        let b_blind: V = 
            (beta + b * x) * x;
        
        trans_seq.push(a_blind.to_tran_elem());
        trans_seq.push(b_blind.to_tran_elem());

    }


    pub fn verify_as_subprotocol<T, U, V>(
        &self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool 
    where
    T: 'static + Copy + Clone + ToTranElem 
        + Add<T, Output = T> + Mul<ZpElement, Output = T>,
    U: 'static + Copy + Clone + ToTranElem + Mul<V, Output = T> 
        + Add<U, Output=U> + Mul<ZpElement, Output = U>,
    V: 'static + Copy + Clone + ToTranElem 
        + Add<V, Output=V> + Mul<ZpElement, Output = V>,     
    {

        let pointer_old = trans_seq.pointer;

        trans_seq.pointer = pointer_old + 13;

        if (
            TranElem::Gt(self.c_com),
            TranElem::Gt(self.a_com),
            TranElem::Gt(self.b_com),
        ) != (
            trans_seq.data[pointer_old],
            trans_seq.data[pointer_old + 1],
            trans_seq.data[pointer_old + 2],
        ) {
            println!("!! Invalid public input when verifying ZkMulScalar");
            return false;
        } 

        
        if let (
            TranElem::Gt(alpha_com),
            TranElem::Gt(beta_com),
            TranElem::Gt(gamma_com),
            TranElem::Gt(delta_com),
            TranElem::Coin(x),
            TranElem::Zp(z_a),
            TranElem::Zp(z_b),
            TranElem::Zp(z_c),
        ) = (
            trans_seq.data[pointer_old + 3],
            trans_seq.data[pointer_old + 4],
            trans_seq.data[pointer_old + 5],
            trans_seq.data[pointer_old + 6],
            trans_seq.data[pointer_old + 7],
            trans_seq.data[pointer_old + 8],
            trans_seq.data[pointer_old + 9],
            trans_seq.data[pointer_old + 10],
        ) {

            let a_blind: U = U::from_tran_elem(
                trans_seq.data[pointer_old + 11]
            );
            let b_blind: V = V::from_tran_elem(
                trans_seq.data[pointer_old + 12]
            );

            let blind_base = srs.blind_base;
            let p_a = 
                (alpha_com + self.a_com * x + blind_base * z_a) * x;
            let p_b = 
                (beta_com + self.b_com * x + blind_base * z_b) * x;
            let p_c =
                (gamma_com + delta_com * x + self.c_com * (x * x) 
                + blind_base * z_c) * (x * x);

            let check1: bool = 
                p_a == a_blind.to_gt(srs);
            let check2: bool =
                p_b == b_blind.to_gt(srs);
            let check3: bool =
                p_c == (a_blind * b_blind).to_gt(srs);

           if check1 || check2 || check3 {
            //    println!("Check passed");
               return true;
               
           } else {
               println!(
                "!! ZkSchnorr equation check failed when verifying ZkSchnorr"
            );
           } 
        } else {
            println!(
                "!! Type check for transcript elements failed when verifying ZkSchnorr"
            );
        }

        return false;
        
    } 

}

/// This is a simplied case of ZkMulScalar
impl <V> ZkSemiMulScalar<V> 
where
    V: 'static + Copy + Clone + Add  + ToTranElem,
{
    pub fn new(
        c_com_value: GtElement, 
        a_com_value: GtElement,
        b_value: V,
    ) -> Self {
        Self {
            c_com: c_com_value,
            a_com: a_com_value,
            b: b_value,
        }
    }

    pub fn prove<T, U>(
        &self, 
        srs: &SRS, 
        trans_seq: &mut TranSeq, 
        a: U, 
        c_tilde: ZpElement, a_tilde: ZpElement, 
    ) where
        T: 'static + Copy + Clone + ToTranElem 
            + Add<T, Output = T> + Mul<ZpElement, Output = T>,
        U: 'static + Copy + Clone + ToTranElem + Mul<V, Output = T> 
            + Add<U, Output=U> + Mul<ZpElement, Output = U>,
     {
        
        trans_seq.push(TranElem::Gt(self.c_com));
        trans_seq.push(TranElem::Gt(self.a_com));
        trans_seq.push(self.b.to_tran_elem());

        let alpha: U = U::rand_general();
        let gamma: T = alpha * self.b;

        let alpha_tilde = ZpElement::rand();
        let gamma_tilde = ZpElement::rand();
     
        let blind_base = srs.blind_base;
        let alpha_com 
            = alpha.to_gt(srs) + alpha_tilde * blind_base;
        let gamma_com
            = gamma.to_gt(srs) + gamma_tilde * blind_base;
         
        trans_seq.push(TranElem::Gt(alpha_com));
        trans_seq.push(TranElem::Gt(gamma_com));
       
        let x = trans_seq.gen_challenge();

        let z_a = 
            - alpha_tilde - x * a_tilde;
        let z_c = 
            - gamma_tilde - x * c_tilde;


        trans_seq.push(TranElem::Zp(z_a));
        trans_seq.push(TranElem::Zp(z_c));

        let a_blind: U = 
            alpha + a * x;
        
        trans_seq.push(a_blind.to_tran_elem());
       
    }

    pub fn verify_as_subprotocol<T, U>(
        &self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool 
    where
    T: 'static + Copy + Clone + ToTranElem 
        + Add<T, Output = T> + Mul<ZpElement, Output = T>,
    U: 'static + Copy + Clone + ToTranElem + Mul<V, Output = T> 
        + Add<U, Output=U> + Mul<ZpElement, Output = U>, 
    {

        let pointer_old = trans_seq.pointer;

        trans_seq.pointer = pointer_old + 9;

        if (
            TranElem::Gt(self.c_com),
            TranElem::Gt(self.a_com),
            self.b.to_tran_elem(),
        ) != (
            trans_seq.data[pointer_old],
            trans_seq.data[pointer_old + 1],
            trans_seq.data[pointer_old + 2],
        ) {
            println!("!! Invalid public input when verifying ZkMulSemiScalar");
            return false;
        } 

        
        if let (
            TranElem::Gt(alpha_com),
            TranElem::Gt(gamma_com),
            TranElem::Coin(x),
            TranElem::Zp(z_a),
            TranElem::Zp(z_c),
        ) = (
            trans_seq.data[pointer_old + 3],
            trans_seq.data[pointer_old + 4],
            trans_seq.data[pointer_old + 5],
            trans_seq.data[pointer_old + 6],
            trans_seq.data[pointer_old + 7],
        ) {

            let a_blind: U = U::from_tran_elem(
                trans_seq.data[pointer_old + 8]
            );
        
            let blind_base = srs.blind_base;
            let p_a = 
                alpha_com + self.a_com * x + blind_base * z_a;
            let p_c = 
                gamma_com + self.c_com * x + blind_base * z_c;
      
            let check1: bool = 
                p_a == a_blind.to_gt(srs);
            let check2: bool =
                p_c == (a_blind * self.b).to_gt(srs);

           if check1 || check2 {
            //    println!("Check passed");
               return true;
               
           } else {
               println!(
                "!! ZkMulSemiScalar equation check failed"
            );
           } 
        } else {
            println!(
                "!! Type check for transcript elements failed when verifying ZkMulSemiScalar"
            );
        }

        return false;
        
    } 


}


#[cfg(test)]
mod tests {
    
    use super::*;
    use crate::setup::SRS;

    use crate::utils::curve::{
        ZpElement, G1Element, G2Element, GtElement,
    };
   
    #[test]
    fn test_zk_scalar() {

        let srs = SRS::new(2);

        let a = ZpElement::rand();
        let b = ZpElement::rand();
        let c = a * b;

        let b_g1 = ZpElement::rand() * G1Element::from(1 as u64);
        let c_g1 = a * b_g1;

        let a_g2 = ZpElement::rand() * G2Element::from(1 as u64);
        let c_gt = a_g2 * b_g1;

        let base = GtElement::from(5 as u64);

        let a_base = a * base;

        let a_tilde = ZpElement::rand();
        let b_tilde = ZpElement::rand();
        let c_tilde = ZpElement::rand();
        
        let blind_base = srs.blind_base;

        let a_com = a.to_gt(&srs) + a_tilde * blind_base;
        let b_com = b.to_gt(&srs) + b_tilde * blind_base;
        let c_com = c.to_gt(&srs) + c_tilde * blind_base;

        let b_g1_com = b_g1.to_gt(&srs) + b_tilde * blind_base;
        let c_g1_com = c_g1.to_gt(&srs) + c_tilde * blind_base;

        let a_g2_com = a_g2.to_gt(&srs) + a_tilde * blind_base;
        let c_gt_com = c_gt.to_gt(&srs) + c_tilde * blind_base;

        let mut trans_seq = TranSeq::new();

        let schnorr = ZkSchnorr::new(
            a_base, base
        );

        let mul_scalar = ZkMulScalar::new(
            c_com, a_com, b_com
        );

        let mul_semi = ZkSemiMulScalar::new(
            c_gt_com, a_g2_com, b_g1
        );

        let mul_scalar_g1 = ZkMulScalar::new(
            c_g1_com, a_com, b_g1_com
        );

        schnorr.prove(&mut trans_seq, a);
        mul_scalar.prove(
            &srs, &mut trans_seq, a, b, c_tilde, a_tilde, b_tilde
        );
        mul_semi.prove(&srs, &mut trans_seq, a_g2, c_tilde, a_tilde);

        mul_scalar_g1.prove(
            &srs, &mut trans_seq, a, b_g1, c_tilde, a_tilde, b_tilde
        );

        let check1 = schnorr.verify_as_subprotocol(&mut trans_seq);
        let check2 = 
            mul_scalar.verify_as_subprotocol::<ZpElement, ZpElement, ZpElement>(
                &srs, &mut trans_seq
            );
        let check3 =
            mul_semi.verify_as_subprotocol::<GtElement, G2Element>(
                &srs, &mut trans_seq
            );
        let check4 = 
            mul_scalar_g1.verify_as_subprotocol::<G1Element, ZpElement, G1Element>(
                &srs, &mut trans_seq
            );

        assert_eq!(check1, true);
        assert_eq!(check2, true);
        assert_eq!(check3, true);
        assert_eq!(check4, true);

    }

    
}
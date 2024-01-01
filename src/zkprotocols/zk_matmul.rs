//! Zero-Knowledge implementation of the zkMatMul protocol
//!
//! Details of this protocol can be found in the DualMatrix paper 
//!
//! To prove that holding three secret matrix \bm{c}, \bm{a}, \bm{b}
//! and three secret Zp elements c_tilde, a_tilde, b_tilde such that 
//! 
//! \bm{c} = \bm{a} \bm{b},
//! C_c = < \vec{G}, \bm{c}, \vec{H} > + c_tilde * blind_base, 
//! C_a = < \vec{G}, \bm{a}, \vec{H} > + a_tilde * blind_base,
//! C_b = < \vec{G}, \bm{b}, \vec{H} > + b_tilde * blind_base,
//! 
use crate::mat::Mat;
use crate::setup::SRS;

use crate::utils::curve::{
    ZpElement, GtElement, G1Element, G2Element,
    ConvertToZp,
};
use crate::utils::dirac::{self, BraKetZp};
use crate::utils::fiat_shamir::{TranElem, TranSeq};

use crate::zkprotocols::zk_left_proj::{ZkLeftProjPoly, ZkLeftProjProof};
use crate::zkprotocols::zk_right_proj::{ZkRightProjPoly, ZkRightProjProof};
use crate::zkprotocols::zk_scalar_proj::{ZkScalarProjPoly, ZkScalarProjCmProof};
use crate::zkprotocols::zk_ip_gt::ZkIpGt;

use crate::zkprotocols::zk_trans::ZkTranSeqProver;


/// Interface for the matMul protocol
pub struct ZkMatMul {
    pub c_blind: GtElement,
    pub a_blind: GtElement,
    pub b_blind: GtElement,
    pub m: usize,
    pub n: usize,
    pub l: usize,
}


impl ZkMatMul {
    pub fn new(
        c_com_value: GtElement, 
        a_com_value: GtElement, 
        b_com_value: GtElement,  
        m_value: usize,
        n_value: usize,
        l_value: usize,
    ) -> Self {
        Self {
            c_blind: c_com_value,
            a_blind: a_com_value,
            b_blind: b_com_value,
            m: m_value,
            n: n_value,
            l: l_value,
         }
    }

    pub fn prove<T, U, V>(
        &self, srs: &SRS, 
        zk_trans_seq: &mut ZkTranSeqProver, 
        mat_c: Mat<T>, 
        mat_a: Mat<U>, 
        mat_b: Mat<V>,
        cache_c: &Vec<G1Element>,
        cache_a: &Vec<G2Element>,
        cache_b: &Vec<G1Element>,
        c_tilde: ZpElement,
        a_tilde: ZpElement,
        b_tilde: ZpElement,
    ) where 
        T: 'static + ConvertToZp,
        Mat<T>: 'static + BraKetZp, 
        U: 'static + ConvertToZp,
        Mat<U>: 'static + BraKetZp, 
        V: 'static + ConvertToZp,
        Mat<V>: 'static + BraKetZp, 
    {
        
        zk_trans_seq.push_without_blinding(TranElem::Gt(self.c_blind));
        zk_trans_seq.push_without_blinding(TranElem::Gt(self.a_blind));
        zk_trans_seq.push_without_blinding(TranElem::Gt(self.b_blind));

        zk_trans_seq.push_without_blinding(TranElem::Size(self.m));
        zk_trans_seq.push_without_blinding(TranElem::Size(self.n));
        zk_trans_seq.push_without_blinding(TranElem::Size(self.l));
      
        let y = zk_trans_seq.gen_challenge();

        let m = self.m;
        let n = self.n;
        let l = self.l;


        if m & (m - 1) != 0 || n & (n - 1) != 0 || l & (l - 1) != 0
            || m != mat_c.shape.0 || m != mat_a.shape.0
            || n != mat_c.shape.1 || n != mat_b.shape.1
            || l != mat_a.shape.1 || l != mat_b.shape.0 
        {
            panic!("Invalid shape when proving MatMul");
        }

        
        let step = y.pow(n as u64);
        let y_l = std::iter::successors(
                Some(ZpElement::from(1 as u64)), 
                |&x| Some(x * step)
            ).take(m).collect::<Vec<ZpElement>>();
        let y_r = std::iter::successors(
                Some(ZpElement::from(1 as u64)),
                |&x| Some(x * y)
           ).take(n).collect::<Vec<ZpElement>>();

        let d = mat_c.braket_zp(&y_l, &y_r);
        let a_y = mat_a.bra_zp(&y_l);
        let b_y = mat_b.ket_zp(&y_r);

        let g_0 = srs.g_hat;
        let h_0 = srs.h_hat;
        let u_0 = g_0 * h_0;
        
        let d_com = d * u_0;
        let a_y_com = g_0 * dirac::inner_product(
            &a_y,
            &srs.h_hat_vec[0..l].to_vec(), 
        );
        let b_y_com = h_0 * dirac::inner_product(
            &b_y,
            &srs.g_hat_vec[0..l].to_vec(), 
        );


        let (d_blind, d_tilde) =
            zk_trans_seq.push_gen_blinding(TranElem::Gt(d_com));
        let (ay_blind, ay_tilde) =
            zk_trans_seq.push_gen_blinding(TranElem::Gt(a_y_com));
        let (by_blind, by_tilde) =
            zk_trans_seq.push_gen_blinding(TranElem::Gt(b_y_com));

        let ip1 = ZkScalarProjPoly::new(
            d_blind, 
            self.c_blind,
            (m,n), 
            y,
        );

        let ip2 = ZkLeftProjPoly::new(
            ay_blind, 
            self.a_blind,
            (m,l), 
            y,
            n,
        );

        let ip3 = ZkRightProjPoly::new(
            by_blind, 
            self.b_blind,
            (l, n), 
            y,
            1,
        );

        let ip4 = ZkIpGt::new(
            d_blind + ay_blind + by_blind,
            l, 
        );


        ip2.prove::<U>(
            srs, zk_trans_seq, 
            &mat_a, cache_a,
            ay_tilde, a_tilde,
        );

        std::mem::drop(mat_a);

        ip3.prove::<V>(
            srs, zk_trans_seq, 
            &mat_b, cache_b, 
            by_tilde, b_tilde
        );

        std::mem::drop(mat_b);

        ip1.prove_cm::<T>(
            srs, zk_trans_seq, 
            &mat_c, cache_c, 
            d_tilde, c_tilde,
        );

        std::mem::drop(mat_c);


        ip4.prove::<ZpElement>(
            srs, zk_trans_seq,
            &a_y, &b_y,
            d_tilde + ay_tilde + by_tilde,
        );
     
    }

    pub fn verify_as_subprotocol(
        &self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        let pointer_old = trans_seq.pointer;
        
        if (
            TranElem::Gt(self.c_blind),
            TranElem::Gt(self.a_blind),
            TranElem::Gt(self.b_blind),
            TranElem::Size(self.m),
            TranElem::Size(self.n),
            TranElem::Size(self.l),
        ) != (
            trans_seq.data[pointer_old],
            trans_seq.data[pointer_old + 1],
            trans_seq.data[pointer_old + 2], 
            trans_seq.data[pointer_old + 3],
            trans_seq.data[pointer_old + 4],
            trans_seq.data[pointer_old + 5],
        ) {
            println!("!! Invalid public input when verifying MatMul");
            return false;
        } 

        let y: ZpElement;
        
        if let TranElem::Coin(y_read) = 
            trans_seq.data[pointer_old + 6] 
        {
            y = y_read;
        } else {
            println!("!! Invalid transcript when verifying MatMul");
            return false;
        }

        let m = self.m;
        let n = self.n;
        let l = self.l;
        if m & (m - 1) != 0 || n & (n - 1) != 0 || l & (l - 1) != 0 {
            panic!("Invalid shape when verifying MatMul");
        }

        let d_blind: GtElement;
        let a_y_blind: GtElement;
        let b_y_blind: GtElement;

        if let (
            TranElem::Gt(d_blind_read),
            TranElem::Gt(a_y_blind_read),
            TranElem::Gt(b_y_blind_read),
        ) = (
            trans_seq.data[pointer_old + 7],
            trans_seq.data[pointer_old + 8],
            trans_seq.data[pointer_old + 9], 
        ) {
            d_blind = d_blind_read;
            a_y_blind = a_y_blind_read;
            b_y_blind = b_y_blind_read;
        } else {
            println!("!! Invalid Transcript type when verifying MatMul");
            return false;
        } 


        trans_seq.pointer = pointer_old + 10 ;

        let ip1 = ZkScalarProjPoly::new(
            d_blind, 
            self.c_blind,
            (m,n), 
            y,
        );

        let ip2 = ZkLeftProjPoly::new(
            a_y_blind, 
            self.a_blind,
            (m,l), 
            y,
            n,
        );

        let ip3 = ZkRightProjPoly::new(
            b_y_blind, 
            self.b_blind,
            (l, n), 
            y,
            1,
        );



        let ip4 = ZkIpGt::new(
            d_blind + a_y_blind + b_y_blind,
            l, 
        );

        let check2 = ip2.verify_as_subprotocol(srs, trans_seq);
        let check3 = ip3.verify_as_subprotocol(srs, trans_seq);
        let check1 = ip1.verify_as_subprotocol_cm(srs, trans_seq);
        let check4 = ip4.verify_as_subprotocol(srs, trans_seq);

        println!("Check of Ip1 in MatMul: {:?}", check1);
        println!("Check of Ip2 in MatMul: {:?}", check2);
        println!("Check of Ip3 in MatMul: {:?}", check3);
        println!("Check of Ip4 in MatMul: {:?}", check4);

        return  check1 && check2 && check3 && check4;
        
    }

    pub fn verify(&self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        if trans_seq.check_fiat_shamir() == false {
            println!("!! Fiat shamir check failed when verifying LeftProj");
            return false;
        }

        self.verify_as_subprotocol(srs, trans_seq)
    }

}


#[cfg(test)]
mod tests {
    
    use super::*;
    use crate::commit_mat::CommitMat;

    use crate::utils::curve::ZpElement;
    use crate::utils::test_data;
   
    #[test]
    fn test_zk_matmul() {

        let srs = SRS::new(64);

        let (c, a, b) = 
            test_data::gen_test_matrices_not_square();

        let (c_com, c_cache_cm) =
            c.commit_col_major(
                &srs.g_hat_vec, &srs.h_hat_vec
            );

        let (a_com, a_cache_rm) = 
            a.commit_row_major(
                &srs.g_hat_vec, &srs.h_hat_vec
            );
        
        let (b_com, b_cache_cm) = 
            b.commit_col_major(
                &srs.g_hat_vec, &srs.h_hat_vec
            );
        
        let c_tilde = ZpElement::rand();
        let a_tilde = ZpElement::rand();
        let b_tilde = ZpElement::rand();

        let blind_base = srs.blind_base;

        let matmul_protocol = ZkMatMul::new(
            c_com + c_tilde * blind_base,
            a_com + a_tilde * blind_base,
            b_com + b_tilde * blind_base, 
            c.shape.0,
            c.shape.1,
            a.shape.1,
        );


        let mut zk_trans_seq = ZkTranSeqProver::new(&srs);

        matmul_protocol.prove::<i128, i64, i64>(
            &srs,
            &mut zk_trans_seq,
            c, a, b,
            &c_cache_cm, &a_cache_rm, &b_cache_cm, 
            c_tilde, a_tilde, b_tilde,
        );

        let mut trans_seq = zk_trans_seq.publish_trans();

        let result = matmul_protocol.verify(
            &srs, &mut trans_seq
        );

        assert_eq!(trans_seq.data.len(), trans_seq.pointer);

        assert_eq!(result, true);

        println!(" * Verification of MatMul passed");

    }

    
}

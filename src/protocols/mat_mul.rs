//! Implementation of the MatMul protocol
//!
//! Details of this protocol can be found in the DualMatrix paper 
//!
//! To prove that holding three secret matrix \bm{c}, \bm{a}, \bm{b} such that
//! 
//! \bm{c} = \bm{a} \bm{b},
//! C_c = < \vec{G}, \bm{c}, \vec{H} >, 
//! C_a = < \vec{G}, \bm{a}, \vec{H} >,
//! C_b = < \vec{G}, \bm{b}, \vec{H} >,
//! 
use crate::mat::Mat;
use crate::setup::SRS;

use crate::utils::curve::{
    ZpElement, GtElement, G1Element, G2Element,
    ConvertToZp,
};
use crate::utils::dirac::{self, BraKetZp};
use crate::utils::fiat_shamir::{TranElem, TranSeq};

use crate::protocols::left_proj::LeftProjPoly;
use crate::protocols::right_proj::RightProjPoly;
use crate::protocols::scalar_proj::ScalarProjPoly;
use crate::protocols::ip_gt::IpGt;

use super::left_proj::LeftProjProof;
use super::right_proj::RightProjProof;
use super::scalar_proj::ScalarProjCmProof;

/// Interface for the matMul protocol
pub struct MatMul {
    pub c_com: GtElement,
    pub a_com: GtElement,
    pub b_com: GtElement,
    pub m: usize,
    pub n: usize,
    pub l: usize,
}


impl MatMul {
    pub fn new(
        c_com_value: GtElement, 
        a_com_value: GtElement, 
        b_com_value: GtElement,  
        m_value: usize,
        n_value: usize,
        l_value: usize,
    ) -> Self {
        Self {
            c_com: c_com_value,
            a_com: a_com_value,
            b_com: b_com_value,
            m: m_value,
            n: n_value,
            l: l_value,
         }
    }

    pub fn prove<T, U, V>(
        &self, srs: &SRS, 
        trans_seq: &mut TranSeq, 
        mat_c: &Mat<T>, 
        mat_a: &Mat<U>, 
        mat_b: &Mat<V>,
        cache_c: &Vec<G1Element>,
        cache_a: &Vec<G2Element>,
        cache_b: &Vec<G1Element>,
    ) where 
        T: 'static + ConvertToZp,
        Mat<T>: 'static + BraKetZp, 
        U: 'static + ConvertToZp,
        Mat<U>: 'static + BraKetZp, 
        V: 'static + ConvertToZp,
        Mat<V>: 'static + BraKetZp, 
    {
        
        trans_seq.push(TranElem::Gt(self.c_com));
        trans_seq.push(TranElem::Gt(self.a_com));
        trans_seq.push(TranElem::Gt(self.b_com));

        trans_seq.push(TranElem::Size(self.m));
        trans_seq.push(TranElem::Size(self.n));
        trans_seq.push(TranElem::Size(self.l));
      
        let y = trans_seq.gen_challenge();

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


        trans_seq.push(TranElem::Gt(d_com));
        trans_seq.push(TranElem::Gt(a_y_com));
        trans_seq.push(TranElem::Gt(b_y_com));

        let ip1 = ScalarProjPoly::new(
            d_com, 
            self.c_com,
            (m,n), 
            y,
        );

        let ip2 = LeftProjPoly::new(
            a_y_com, 
            self.a_com,
            (m,l), 
            y,
            n,
        );

        let ip3 = RightProjPoly::new(
            b_y_com, 
            self.b_com,
            (l, n), 
            y,
            1,
        );

        let ip4 = IpGt::new(
            d_com + a_y_com + b_y_com,
            l, 
        );

        ip1.prove_cm::<T>(srs, trans_seq, mat_c, cache_c);
        ip2.prove::<U>(srs, trans_seq, mat_a, cache_a);
        ip3.prove::<V>(srs, trans_seq, mat_b, cache_b);
        ip4.prove::<ZpElement>(srs, trans_seq,&a_y, &b_y);
     
    }

    pub fn verify_as_subprotocol(
        &self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        let pointer_old = trans_seq.pointer;
        
        if (
            TranElem::Gt(self.c_com),
            TranElem::Gt(self.a_com),
            TranElem::Gt(self.b_com),
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

        let d_com: GtElement;
        let a_y_com: GtElement;
        let b_y_com: GtElement;

        if let (
            TranElem::Gt(d_com_read),
            TranElem::Gt(a_y_com_read),
            TranElem::Gt(b_y_com_read),
        ) = (
            trans_seq.data[pointer_old + 7],
            trans_seq.data[pointer_old + 8],
            trans_seq.data[pointer_old + 9], 
        ) {
            d_com = d_com_read;
            a_y_com = a_y_com_read;
            b_y_com = b_y_com_read;
        } else {
            println!("!! Invalid Transcript type when verifying MatMul");
            return false;
        } 


        trans_seq.pointer = pointer_old + 10 ;

        let ip1 = ScalarProjPoly::new(
            d_com, 
            self.c_com,
            (m,n), 
            y,
        );

        let ip2 = LeftProjPoly::new(
            a_y_com, 
            self.a_com,
            (m,l), 
            y,
            n,
        );

        let ip3 = RightProjPoly::new(
            b_y_com, 
            self.b_com,
            (l, n), 
            y,
            1,
        );



        let ip4 = IpGt::new(
            d_com + a_y_com + b_y_com,
            l, 
        );

        let check1 = ip1.verify_as_subprotocol_cm(srs, trans_seq);
        let check2 = ip2.verify_as_subprotocol(srs, trans_seq);
        let check3 = ip3.verify_as_subprotocol(srs, trans_seq);
        let check4 = ip4.verify_as_subprotocol(srs, trans_seq);

        println!("check 1: {:?}", check1);
        println!("check 2: {:?}", check2);
        println!("check 3: {:?}", check3);
        println!("check 4: {:?}", check4);

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
    use crate::experiment_data;
   
    #[test]
    fn test_matmul() {

        let srs = SRS::new(64);

        let (c, a, b) = 
            experiment_data::gen_matrices_dense(32);

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


        let matmul_protocol = MatMul::new(
            c_com,
            a_com,
            b_com, 
            32,
            32,
            32,
        );


        let mut trans_seq = TranSeq::new();

        matmul_protocol.prove::<i128, i64, i64>(
            &srs,
            &mut trans_seq,
            &c, &a, &b,
            &c_cache_cm, &a_cache_rm, &b_cache_cm, 
        );

        let result = matmul_protocol.verify(
            &srs, &mut trans_seq
        );

        assert_eq!(trans_seq.data.len(), trans_seq.pointer);

        assert_eq!(result, true);

        println!(" * Verification of MatMul passed");

    }

    
}

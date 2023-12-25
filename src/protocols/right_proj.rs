///! Implementation of the Scalar Projection protocol
///!
///! Details of this protocol can be found in the DualMatrix paper 
///!
///! To prove that holding a secret matrix \bm{a} such that
///! 
///! C_a = < \vec{G}, \bm{a}, \vec{H} >
///! C_c = (l^T \bm{a} r) e(\hat{G}, \hat{H})
/// 
use crate::mat::Mat;
use crate::setup::SRS;

use crate::utils::curve::{
    ZpElement, GtElement, G1Element,
    ConvertToZp,
};
use crate::utils::dirac::{self, BraKetZp};
use crate::utils::fiat_shamir::{TranElem, TranSeq};
use crate::utils::xi;

use crate::protocols::pip::{PipG1, PipG2};



pub struct RightProj {
    pub c_com: GtElement,
    pub a_com: GtElement,
    pub shape: (usize, usize),
    pub r_vec: Vec<ZpElement>,
}

pub struct RightProjPoly {
    pub c_com: GtElement,
    pub a_com: GtElement,
    pub shape: (usize, usize),
    pub y: ZpElement, 
}

pub trait RightProjInterface {
    fn get_c_com(&self) -> GtElement;
    fn get_a_com(&self) -> GtElement;
    fn get_shape(&self) -> (usize, usize);
    fn reduce_r(&self, challenges: &Vec<ZpElement>) -> ZpElement;
    fn get_r_vec(&self) -> Vec<ZpElement>;
}

impl RightProjInterface for RightProj {
    fn reduce_r(&self, challenges: &Vec<ZpElement>) -> ZpElement {
        xi::reduce_from_challenges(challenges, &self.get_r_vec())
    }


    fn get_r_vec(&self) -> Vec<ZpElement> {
        self.r_vec.clone()
    }

    fn get_c_com(&self) -> GtElement {
        self.c_com
    }

    fn get_a_com(&self) -> GtElement {
        self.a_com
    }

    fn get_shape(&self) -> (usize, usize) {
        self.shape
    }

}

impl RightProjInterface for RightProjPoly {
    fn reduce_r(&self, challenges: &Vec<ZpElement>) -> ZpElement {
        xi::phi_s(
            self.y, challenges,
            0 as usize, 
            1 as usize,
        )
    }

    fn get_r_vec(&self) -> Vec<ZpElement> {
        let y = self.y;
        let n = self.get_shape().1;
        
        std::iter::successors(
            Some(ZpElement::from(1 as u64)),
            |&x| Some(y * x),
        ).take(n).collect::<Vec<ZpElement>>()
    }
    fn get_c_com(&self) -> GtElement {
        self.c_com
    }

    fn get_a_com(&self) -> GtElement {
        self.a_com
    }

    fn get_shape(&self) -> (usize, usize) {
        self.shape
    }

}

impl RightProj {
    pub fn new(
        c_com_value: GtElement, 
        a_com_value: GtElement,  
        shape_value: (usize, usize),
        r_vec: &Vec<ZpElement>,
    ) -> Self {
        Self {
            c_com: c_com_value,
            a_com: a_com_value,
            shape: shape_value,
            r_vec: r_vec.clone(),
        }
    }
}

impl RightProjPoly {
    pub fn new(
        c_com_value: GtElement, 
        a_com_value: GtElement,  
        shape_value: (usize, usize),
        y_value: ZpElement,
    ) -> Self {
        Self {
            c_com: c_com_value,
            a_com: a_com_value,
            shape: shape_value,
            y: y_value,
        }
    }
}

impl RightProjProof for RightProj {}
impl RightProjProof for RightProjPoly {}


/// Use the column major caches for Right projection
/// 
pub trait RightProjProof: RightProjInterface {

    fn prove<T>(
        &self, srs: &SRS, 
        trans_seq: &mut TranSeq, 
        mat_a: &Mat<T>, 
        a_cache: &Vec<G1Element>,
    ) where 
        T: 'static + ConvertToZp,
        Mat<T>: 'static + BraKetZp, 
    {

        trans_seq.push(TranElem::Gt(self.get_c_com()));
        trans_seq.push(TranElem::Gt(self.get_a_com()));
        trans_seq.push(TranElem::Size(self.get_shape().0));
        trans_seq.push(TranElem::Size(self.get_shape().1));

        let x = trans_seq.gen_challenge();

        let m = self.get_shape().0;
        let n = self.get_shape().1;

        if m & (m - 1) != 0 || n & (n - 1) != 0 
            || m != mat_a.shape.0 || n != mat_a.shape.1 
        {
            panic!("Invalid shape when proving ScalarProjRm");
        }

        let log_m = (m as u64).ilog2() as usize;
        let log_n = (n as u64).ilog2() as usize;

        let h_0 = srs.h_hat;
        let r_vec = self.get_r_vec();

       
        let mut capital_a_current = a_cache[0..n].to_vec();
        let mut h_vec_current = srs.h_hat_vec[0..n].to_vec();
        let mut r_current = r_vec[0..n].to_vec();
        
        let mut challenges_n: Vec<ZpElement> = Vec::new();
        let mut challenges_inv_n: Vec<ZpElement> = Vec::new();
 

        for j in 0..log_n {
            let current_len = n / 2usize.pow(j as u32);
            

            let capital_a_left = 
                capital_a_current[0..current_len/2].to_vec();
            let capital_a_right = 
                capital_a_current[current_len/2..current_len].to_vec();
            
            let r_left = 
                r_current[0..current_len/2].to_vec();
            let r_right = 
                r_current[current_len/2..current_len].to_vec();
            

            let h_left = 
                h_vec_current[0..current_len/2].to_vec();
            let h_right = 
                h_vec_current[current_len/2..current_len].to_vec();

            let l_tr = 
                dirac::inner_product(&capital_a_left, &h_right)
                + h_0 * x * dirac::inner_product(&capital_a_left, &r_right);
            let r_tr = 
                dirac::inner_product(&capital_a_right, &h_left)
                + h_0 * x * dirac::inner_product(&capital_a_right, &r_left);

            trans_seq.push(TranElem::Gt(l_tr));
            trans_seq.push(TranElem::Gt(r_tr));
            
            let x_j = trans_seq.gen_challenge();
            let x_j_inv = x_j.inv();

            challenges_n.push(x_j);
            challenges_inv_n.push(x_j_inv);

            capital_a_current = dirac::vec_addition(
                &capital_a_left,
                &dirac::vec_scalar_mul(
                    &capital_a_right, x_j_inv),
            );

            h_vec_current = dirac::vec_addition(
                &h_left,
                &dirac::vec_scalar_mul(
                    &h_right, x_j),
            );

            r_current = dirac::vec_addition(
                &r_left,
                &dirac::vec_scalar_mul(
                    &r_right, x_j),
            );

        }

        let xi_n_inv = xi::xi_from_challenges(&challenges_inv_n);
        let a_xi_inv = mat_a.bra_zp(&xi_n_inv);


        let h_reduce = h_vec_current[0];
        let r_reduce = r_current[0];
        let w = h_reduce + x * r_reduce * h_0;

        let mut a_current = dirac::vec_convert_to_zp_vec(
            &a_xi_inv);
        let mut g_vec_current = srs.g_hat_vec[0..m].to_vec();
        
        let mut challenges_m: Vec<ZpElement> = Vec::new();
        let mut challenges_inv_m: Vec<ZpElement> = Vec::new();
        

        for j in 0..log_m {
            let current_len = m / 2usize.pow(j as u32);
            
            let a_left = 
                a_current[0..current_len/2].to_vec();
            let a_right = 
                a_current[current_len/2..current_len].to_vec();
            

            let g_left = 
                g_vec_current[0..current_len/2].to_vec();
            let g_right = 
                g_vec_current[current_len/2..current_len].to_vec();

            let l_tr = 
                w * dirac::inner_product(&a_left, &g_right);
            let r_tr = 
                w * dirac::inner_product(&a_right, &g_left);
                
            trans_seq.push(TranElem::Gt(l_tr));
            trans_seq.push(TranElem::Gt(r_tr));
            
            let x_j = trans_seq.gen_challenge();
            let x_j_inv = x_j.inv();

            challenges_m.push(x_j);
            challenges_inv_m.push(x_j_inv);

            a_current = dirac::vec_addition(
                &a_left,
                &dirac::vec_scalar_mul(
                    &a_right, x_j_inv),
            );

            g_vec_current = dirac::vec_addition(
                &g_left,
                &dirac::vec_scalar_mul(
                    &g_right, x_j),
            );

        }

        let a_reduce = a_current[0];

        let g_reduce = g_vec_current[0];

        trans_seq.push(TranElem::Zp(a_reduce));
        trans_seq.push(TranElem::G1(g_reduce));
        trans_seq.push(TranElem::G2(h_reduce));

        let pip_g1 = PipG1::new(
            g_reduce, &challenges_m
        );

        let pip_g2 = PipG2::new(
            h_reduce, &challenges_n
        );

        pip_g1.prove(srs, trans_seq);
        pip_g2.prove(srs, trans_seq);

    }

    fn verify_as_subprotocol(
        &self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        let pointer_old = trans_seq.pointer;
        
        if (
            TranElem::Gt(self.get_c_com()),
            TranElem::Gt(self.get_a_com()),
            TranElem::Size(self.get_shape().0),
            TranElem::Size(self.get_shape().1)
        ) != (
            trans_seq.data[pointer_old],
            trans_seq.data[pointer_old + 1],
            trans_seq.data[pointer_old + 2], 
            trans_seq.data[pointer_old + 3],
        ) {
            println!("!! Invalid public input when verifying RightProj");
            return false;
        } 

        let x: ZpElement;
        
        if let TranElem::Coin(x_read) = 
            trans_seq.data[pointer_old + 4] 
        {
            x = x_read;
        } else {
            println!("!! Invalid transcript when verifying RightProj");
            return false;
        }

        let m = self.get_shape().0;
        let n = self.get_shape().1;
        if m & (m - 1) != 0 || n & (n - 1) != 0 {
            panic!("Invalid shape when verifying RightProj");
        }
        let log_m = (m as u64).ilog2() as usize;
        let log_n = (n as u64).ilog2() as usize;

        trans_seq.pointer = 
            pointer_old + 5 + 3* log_m + 3* log_n + 3 ;

        let mut current_pointer = pointer_old + 5;
        let mut lhs: GtElement = 
            x * self.get_c_com().clone() + self.get_a_com().clone();
        
        let mut challenges_n: Vec<ZpElement> = Vec::new();
        let mut challenges_inv_n: Vec<ZpElement> = Vec::new();

        for _ in 0..log_n {

            if let (
                TranElem::Gt(l_tr),
                TranElem::Gt(r_tr),
                TranElem::Coin(x_j),
            ) = (
                trans_seq.data[current_pointer],
                trans_seq.data[current_pointer + 1],
                trans_seq.data[current_pointer + 2],
            ) {
                
                let x_j_inv = x_j.inv();
                lhs = lhs + l_tr * x_j + r_tr * x_j_inv;
                challenges_n.push(x_j);
                challenges_inv_n.push(x_j_inv);

            } else {
                println!("!! Invalid transcript when verifying RightProj");
                return false;
            }

            current_pointer += 3;
        }

        let mut challenges_m: Vec<ZpElement> = Vec::new();
        let mut challenges_inv_m: Vec<ZpElement> = Vec::new();

        for _ in 0..log_m {

            if let (
                TranElem::Gt(l_tr),
                TranElem::Gt(r_tr),
                TranElem::Coin(x_j),
            ) = (
                trans_seq.data[current_pointer],
                trans_seq.data[current_pointer + 1],
                trans_seq.data[current_pointer + 2],
            ) {
                
                let x_j_inv = x_j.inv();
                lhs = lhs + l_tr * x_j + r_tr * x_j_inv;
                challenges_m.push(x_j);
                challenges_inv_m.push(x_j_inv);

            } else {
                println!("!! Invalid transcript when verifying RightProj");
                return false;
            }

            current_pointer += 3;
        }


        let xi_r = self.reduce_r(&challenges_n);
        

        if let (
            TranElem::Zp(a_reduce),
            TranElem::G1(g_reduce),
            TranElem::G2(h_reduce),
        ) = (
            trans_seq.data[current_pointer],
            trans_seq.data[current_pointer+1],
            trans_seq.data[current_pointer+2],
        ) {

            let rhs = 
                a_reduce * x * xi_r * srs.h_hat * g_reduce
                + a_reduce * g_reduce * h_reduce;
            if lhs == rhs {
                
                let pip_g1 = PipG1::new(
                    g_reduce, &challenges_m
                );
        
                let pip_g2 = PipG2::new(
                    h_reduce, &challenges_n
                );

                let check_1 = pip_g1.verify_as_subprotocol(srs, trans_seq);
                let check_2 = pip_g2.verify_as_subprotocol(srs, trans_seq);

                return check_1 && check_2;

            }
        } else {
            println!("!! Invalid transcript when verifying RightProj");
            return false;
        }
    
        println!("!! Verification of RightProj failed");
        return false;
        
    }

    fn verify(&self, srs: &SRS, trans_seq: &mut TranSeq
    ) -> bool {

        if trans_seq.check_fiat_shamir() == false {
            println!("!! Fiat shamir check failed when verifying RightProj");
            return false;
        }

        self.verify_as_subprotocol(srs, trans_seq)
    }

}


#[cfg(test)]
mod tests {
    
    use super::*;

    use rand::Rng;

    use crate::commit_mat::CommitMat;
    use crate::utils::curve::Zero;
   
    #[test]
    fn test_right_proj() {

        let srs = SRS::new(64);

        let m = 16 as usize;
        let n = 32 as usize;

        let a_data = (0..m).into_iter()
            .zip((0..n).into_iter())
            .map(|(i, j)| 
                (i, j, rand::thread_rng().gen::<u32>() as i64 )
            ).collect::<Vec<(usize, usize, i64)>>();

        
        let a = Mat::new_from_data_vec(
            "a_test", 
            (m, n), 
            a_data
        );

        let (a_com, a_cache_cm) 
            = a.commit_col_major(
                &srs.g_hat_vec, &srs.h_hat_vec
            );

        let r_vec = (0..n).map(|_| 
            ZpElement::rand()
        ).collect::<Vec<ZpElement>>();

        let c = a.ket_zp(&r_vec);

        let c_com = 
            srs.h_hat 
            * dirac::inner_product(&c, &srs.g_hat_vec);


        let right_proj_protocol = RightProj::new(
            c_com,
            a_com, 
            (m, n), 
            &r_vec, 
        );


        let mut trans_seq_rm = TranSeq::new();

        right_proj_protocol.prove::<i64>(
            &srs, 
            &mut trans_seq_rm, 
            &a,
            &a_cache_cm
        );

        let result_cm = right_proj_protocol.verify_as_subprotocol(
            &srs, &mut trans_seq_rm
        );

        assert_eq!(trans_seq_rm.data.len(), trans_seq_rm.pointer);

        assert_eq!(result_cm, true);

        println!(" * Verification of RightProj passed");

        let y = ZpElement::rand();

        let mut right_proj_poly = RightProjPoly::new(
            GtElement::zero(),
            a_com, 
            (m, n), 
            y,
        );

        let y_l = right_proj_poly.get_r_vec();


        let ay = a.ket_zp(&y_l);
        let ay_com =
            srs.h_hat 
            * dirac::inner_product(&ay, &srs.g_hat_vec);

        right_proj_poly.c_com = ay_com;



        let mut trans_seq_cm = TranSeq::new();

        right_proj_poly.prove::<i64>(
            &srs, 
            &mut trans_seq_cm, 
            &a,
            &a_cache_cm);

        let result_cm = right_proj_poly.verify(
            &srs, &mut trans_seq_cm
        );

        assert_eq!(result_cm, true);

        println!(" * Verification of RightProjPoly passed");


    }

    
}

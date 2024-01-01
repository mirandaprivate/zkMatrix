//! Utility functions for the braket operation
//! 
//! Parallelized whenever possible
//! 
use core::fmt::Debug;
use std::marker::{Send, Sync};
use std::sync::{Arc, Mutex};
use std::thread;

use rayon::{ThreadPoolBuilder, prelude::*};


use crate::utils::curve::{
    Add, Mul, AddAssign, Zero, 
    ZpElement, G1Element, G2Element, GtElement,
    ConvertToZp,
};

use crate::mat::Mat;

use crate::config::NUM_THREADS;

use super::curve::Double;

pub trait BraKet{
    fn bra(&self, v_base: &Vec<G1Element>) -> Vec<G1Element>;
    
    fn ket(&self, v_base: &Vec<G2Element>) -> Vec<G2Element>;
    
    fn braket(
        &self, g_base: &Vec<G1Element>, h_base: &Vec<G2Element>
    ) -> GtElement {
        let left_proj = self.bra(g_base);
        inner_product(&left_proj, h_base)
    }
}

pub trait BraKetZp {
    fn bra_zp(&self, v_base: &Vec<ZpElement>) -> Vec<ZpElement>;
    fn ket_zp(&self, v_base: &Vec<ZpElement>) -> Vec<ZpElement>;

    fn braket_zp(
        &self, g_base: &Vec<ZpElement>, h_base: &Vec<ZpElement>
    ) -> ZpElement {
        let left_proj = self.bra_zp(g_base);
        inner_product(&left_proj, h_base)
    }
}

impl BraKet for Mat<u64> {
    fn bra(&self, v_base: &Vec<G1Element>) -> Vec<G1Element> {
        bra_opt_u64(&self, v_base)
    }

    fn ket(&self, v_base: &Vec<G2Element>) -> Vec<G2Element> {
        ket_opt_u64(&self, v_base)
    }

}

impl BraKet for Mat<i64> {
    fn bra(&self, v_base: &Vec<G1Element>) -> Vec<G1Element> {
        bra_opt_i64(&self, v_base)
    }

    fn ket(&self, v_base: &Vec<G2Element>) -> Vec<G2Element> {
        ket_opt_i64(&self, v_base)
    }
}

impl BraKet for Mat<i128> {
    fn bra(&self, v_base: &Vec<G1Element>) -> Vec<G1Element> {
        bra_opt_i128(&self, v_base)
    }

    fn ket(&self, v_base: &Vec<G2Element>) -> Vec<G2Element> {
        ket_opt_i128(&self, v_base)
    }
}

impl BraKet for Mat<ZpElement> {
    fn bra(&self, v_base: &Vec<G1Element>) -> Vec<G1Element> {
        proj_left(&self, v_base)
    }

    fn ket(&self, v_base: &Vec<G2Element>) -> Vec<G2Element> {
        proj_right(&self, v_base)
    }
}


impl BraKetZp for Mat<i64> {
    fn bra_zp(&self, v_base: &Vec<ZpElement>) -> Vec<ZpElement> {
        // proj_left_opt_i64_to_zp(&self, v_base)
        proj_left(&self, v_base)
    }

    fn ket_zp(&self, v_base: &Vec<ZpElement>) -> Vec<ZpElement> {
        // proj_right_opt_i64_to_zp(&self, v_base)
        proj_right(&self, v_base)

    }
}

impl BraKetZp for Mat<i128> {
    fn bra_zp(&self, v_base: &Vec<ZpElement>) -> Vec<ZpElement> {
        // proj_left_opt_i128_to_zp(&self, v_base)
        proj_left(&self, v_base)
    }

    fn ket_zp(&self, v_base: &Vec<ZpElement>) -> Vec<ZpElement> {
        // proj_right_opt_i128_to_zp(&self, v_base)
        proj_right(&self, v_base)
    }
}

impl BraKetZp for Mat<ZpElement> {
    fn bra_zp(&self, v_base: &Vec<ZpElement>) -> Vec<ZpElement> {
        proj_left(&self, v_base)
    }

    fn ket_zp(&self, v_base: &Vec<ZpElement>) -> Vec<ZpElement> {
        proj_right(&self, v_base)
    }
}

pub fn inner_product<'a, T, U, V>(vec_a: &Vec<T>, vec_b: &Vec<U>) -> V
where
    T: 'a + Clone + Copy + Send + Sync + Mul<U, Output =V>,
    U: 'a + Clone + Copy + Send + Sync,
    V: 'a + Clone + Copy + Send + Sync + Add + Zero,
{
    let length = std::cmp::min(vec_a.len(), vec_b.len());
    let vec_a: Vec<T> = vec_a[..length].into();
    let vec_b: Vec<U> = vec_b[..length].into();


    let pool = ThreadPoolBuilder::new()
        .num_threads(NUM_THREADS)
        .build()
        .unwrap();

    let inner_product = pool.install(|| 
        {
            vec_a.par_iter()
            .zip(vec_b.par_iter())
            .fold(
                || V::zero(),
                |acc, (&a_val, &b_val)| acc + a_val * b_val,
            )
            .reduce(
                || V::zero(),
                |acc, val| acc + val,
            )
        }
    );

    return inner_product;

}

pub fn vec_scalar_mul<T> (vector: &Vec<T>, scalar: ZpElement) -> Vec<T>
where
    T: 'static + Clone + Copy + Send + Sync + Mul<ZpElement, Output = T>,
{

    let pool = ThreadPoolBuilder::new()
        .num_threads(NUM_THREADS)
        .build()
        .unwrap();

    let result = pool.install(|| 
        {
            vector.par_iter()
            .map(|&val| val * scalar)
            .collect()
        }
    );

    return result;

}



pub fn vec_addition<T> (vec_a: &Vec<T>, vec_b: &Vec<T>) -> Vec<T>
where
    T: 'static + Clone + Copy + Send + Sync + Add<Output = T>, 
{

    let pool = ThreadPoolBuilder::new()
        .num_threads(NUM_THREADS)
        .build()
        .unwrap();

    let result = pool.install(|| 
        {
            vec_a.par_iter()
            .zip(vec_b.par_iter())
            .map(|(&a, &b)| a + b)
            .collect()
        }
    );

    return result;

}


pub fn vec_convert_to_zp_vec<T> (vec_a: &Vec<T>) -> Vec<ZpElement> 
where 
    T: 'static + ConvertToZp, 
{

    let pool = ThreadPoolBuilder::new()
        .num_threads(NUM_THREADS)
        .build()
        .unwrap();

    let result = pool.install(|| 
        {
            vec_a.par_iter()
            .map(|&a| a.to_zp())
            .collect()
        }
    );

    return result;

}

pub fn proj_left<T, U, V>(mat_a: &Mat<T>, v_base: &Vec<U>) -> Vec<V>
where
    T: 'static + Clone + Copy + Send + Sync + Mul<U, Output = V>,
    U: 'static + Clone + Copy + Send + Sync + Mul<T, Output = V>,
    V: 'static + Clone + Copy + Send + Sync + Debug + Add + AddAssign + Zero,
{
    let a_data  = &mat_a.data;
    let n_row = mat_a.shape.0;
    let n_col = mat_a.shape.1;
    let v_base_arc = Arc::new(v_base[..n_row].to_vec());

    let result = Arc::new(
        Mutex::new(vec![V::zero(); n_col]));

    let chunks: Vec<_> = a_data.chunks(a_data.len() / NUM_THREADS)
        .collect();

    let mut handles = Vec::new();
    for chunk in chunks {
        let chunk = chunk.to_owned();
        let v_base_clone = Arc::clone(&v_base_arc);
        let handle = thread::spawn(
            move || {
                let mut local_result = vec![V::zero(); n_col];
                for &(row, col, val) in &chunk {
                    local_result[col] += v_base_clone[row] * val;
                }
                local_result
            });
        handles.push(handle);
    }

    for handle in handles {
        let thread_result = handle.join().unwrap();
        let mut result = result.lock().unwrap();
        for i in 0..n_col {
            result[i] += thread_result[i];
        }
    }

    Arc::try_unwrap(result).unwrap().into_inner().unwrap()

    // // The following is slower than the current version
    // // Both version works.
    // let a_data  = &mat_a.data;
    // let v_base_arc = Arc::new(v_base.clone());
    
    // let num_nonzero_entries = a_data.len();
    // let v_len = v_base.len();

    // let result = a_data.chunks(
    //     num_nonzero_entries / NUM_THREADS)
    //     .map(|chunk| {
    //         let v_arc_clone = Arc::clone(&v_base_arc);
    //         let chunk = chunk.to_owned();
    //         thread::spawn(
    //             move || 
    //             chunk.iter().fold(
    //                 vec![V::zero(); v_len],
    //                 |mut acc, &(row, col, val)| {
    //                     acc[col] += v_arc_clone[row] * val;
    //                     acc
    //                 }))
    //     })
    //     .fold(
    //         vec![V::zero(); v_len],
    //         |acc, handle| 
    //         handle.join().unwrap().iter().enumerate().map(
    //             |(i, &l)| l + acc[i]
    //         ).collect()
    //     );
    // result
}

pub fn proj_right<T, U, V>(mat_a: &Mat<T>, v_base: &Vec<U>) -> Vec<V>
where
    T: 'static + Clone + Copy + Send + Sync + Mul<U, Output = V>,
    U: 'static + Clone + Copy + Send + Sync + Mul<T, Output = V>,
    V: 'static + Clone + Copy + Send + Sync + Debug + Add + AddAssign + Zero,
{
    let a_data  = &mat_a.data;
    let n_row = mat_a.shape.0;
    let n_col = mat_a.shape.1;
    let v_base_arc = Arc::new(v_base[..n_col].to_vec());

    let result = Arc::new(Mutex::new(vec![V::zero(); n_row]));

    let chunks: Vec<_> = a_data.chunks(a_data.len() / NUM_THREADS)
        .collect();

    let mut handles = Vec::new();
    for chunk in chunks {
        let chunk = chunk.to_owned();
        let v_base_clone = Arc::clone(&v_base_arc);
        let handle = thread::spawn(
            move || {
                let mut local_result = vec![V::zero(); n_row];
                for &(row, col, val) in &chunk {
                    local_result[row] += v_base_clone[col] * val;
                }
                local_result
            });
        handles.push(handle);
    }

    for handle in handles {
        let thread_result = handle.join().unwrap();
        let mut result = result.lock().unwrap();
        for i in 0..n_row {
            result[i] += thread_result[i];
        }
    }

    Arc::try_unwrap(result).unwrap().into_inner().unwrap()
}

pub fn dirac<T, U, V, R, W>(
    g_base: &Vec<T>,
    mat_a: &Mat<U>, 
    h_base: &Vec<R>,
) -> W
where
    T: 'static + Clone + Copy + Send + Sync + Mul<U, Output = V>,
    U: 'static + Clone + Copy + Send + Sync + Mul<T, Output = V>,
    V: 'static + Clone + Copy + Send + Sync + Debug + Add + AddAssign + Zero
        + Mul<R, Output = W>,
    R: 'static + Clone + Copy + Send + Sync,
    W: 'static + Clone + Copy + Send + Sync + Add + Zero,
{
    let proj_left = proj_left(mat_a, g_base);
    let n_col = mat_a.shape.1;
    let result = inner_product(&proj_left, &h_base[..n_col].to_vec());
    
    result
}


pub fn bra_opt_u64(mat_a: &Mat<u64>, v_base: &Vec<G1Element>) -> Vec<G1Element> {
    
    let a_data  = &mat_a.data;
    let n_row = mat_a.shape.0;
    let n_col = mat_a.shape.1;
    let mut v_base_mut = v_base[..n_row].to_vec();

    let result_mutexes: Vec<Mutex<G1Element>> = (0..n_col).map(
        |_| Mutex::new(G1Element::zero())
        ).collect();
    
    let pool = ThreadPoolBuilder::new()
        .num_threads(NUM_THREADS)
        .build()
        .unwrap();
    
    for j_bit in 0..64{
        let v_base_arc = Arc::new(v_base_mut.clone());
        
        pool.install(|| {
            a_data.par_iter().for_each(|&(row, col, val)| {
                let v_base_clone = Arc::clone(&v_base_arc);
                let result_mutex = &result_mutexes[col];
                let mut result = result_mutex.lock().unwrap();
                if (val >> j_bit) & 1 == 1 {
                    *result += v_base_clone[row];
                }
            });
        });

        pool.install(|| {
            v_base_mut.par_iter_mut().for_each(|v| {
                *v = v.double();
            });
        });

    }

    // let v_base_arc = Arc::new(v_base_mut.clone());
    // pool.install(|| {
    //     a_data.par_iter().for_each(|&(row, col, val)| {
    //         let v_base_clone = Arc::clone(&v_base_arc);
    //         let result_mutex = &result_mutexes[col];
    //         let mut result = result_mutex.lock().unwrap();
    //         *result += v_base_clone[row] * val;
    //     });
    // });

    result_mutexes.into_iter().map(
        |mutex| mutex.into_inner().unwrap()
    ).collect()

}



pub fn ket_opt_u64(mat_a: &Mat<u64>, v_base: &Vec<G2Element>) -> Vec<G2Element> {
    
    let a_data  = &mat_a.data;
    let n_row = mat_a.shape.0;
    let n_col = mat_a.shape.1;
    let mut v_base_mut = v_base[..n_col].to_vec();

    let result_mutexes: Vec<Mutex<G2Element>> = (0..n_row).map(
        |_| Mutex::new(G2Element::zero())
    ).collect();
    
    let pool = ThreadPoolBuilder::new()
        .num_threads(NUM_THREADS)
        .build()
        .unwrap();
    
    for j_bit in 0..64{
        let v_base_arc = Arc::new(v_base_mut.clone());
        
        pool.install(|| {
            a_data.par_iter().for_each(|&(row, col, val)| {
                let v_base_clone = Arc::clone(&v_base_arc);
                let result_mutex = &result_mutexes[row];
                let mut result = result_mutex.lock().unwrap();
                if (val >> j_bit) & 1 == 1 {
                    *result += v_base_clone[col];
                }
            });
        });

        pool.install(|| {
            v_base_mut.par_iter_mut().for_each(|v| {
                *v = v.double();
            });
        });

    }

    // let v_base_arc = Arc::new(v_base_mut.clone());
    // pool.install(|| {
    //     a_data.par_iter().for_each(|&(row, col, val)| {
    //         let v_base_clone = Arc::clone(&v_base_arc);
    //         let result_mutex = &result_mutexes[col];
    //         let mut result = result_mutex.lock().unwrap();
    //         *result += v_base_clone[row] * val;
    //     });
    // });

    result_mutexes.into_iter().map(
        |mutex| mutex.into_inner().unwrap()
    ).collect()

}

pub fn bra_opt_i64(mat_a: &Mat<i64>, v_base: &Vec<G1Element>) -> Vec<G1Element> {

    let n_row = mat_a.shape.0;
    let n_col = mat_a.shape.1;
    let mut v_base_mut = v_base[..n_row].to_vec();
       
    let pool = ThreadPoolBuilder::new()
        .num_threads(NUM_THREADS)
        .build()
        .unwrap();

    let a_abs_sign: Vec<(usize, usize, u64, bool)> =  pool.install(
            || {
                mat_a.data.par_iter().map(|&(row, col, val)| (
                        row, col, val.abs() as u64, val < 0,
                )).collect()
            }
    );

    let result_mutexes: Vec<Mutex<G1Element>> = (0..n_col).map(
        |_| Mutex::new(G1Element::zero())
    ).collect();
    
    
    for j_bit in 0..64{
        let v_base_arc = Arc::new(v_base_mut.clone());
        
        pool.install(|| {
            a_abs_sign.par_iter().for_each(|&(row, col, val, sign)| {
                let v_base_clone = Arc::clone(&v_base_arc);
                let result_mutex = &result_mutexes[col];
                let mut result = result_mutex.lock().unwrap();
                if (val >> j_bit) & 1 == 1 {
                    if sign {
                        *result += - v_base_clone[row];
                    } else {
                        *result += v_base_clone[row];
                    }
                }
            });
        });

        pool.install(|| {
            v_base_mut.par_iter_mut().for_each(|v| {
                *v = v.double();
            });
        });

    }


    result_mutexes.into_iter().map(
        |mutex| mutex.into_inner().unwrap()
    ).collect()
    
}


pub fn ket_opt_i64(mat_a: &Mat<i64>, v_base: &Vec<G2Element>) -> Vec<G2Element> {

    let n_row = mat_a.shape.0;
    let n_col = mat_a.shape.1;
    let mut v_base_mut = v_base[..n_col].to_vec();
       
    let pool = ThreadPoolBuilder::new()
        .num_threads(NUM_THREADS)
        .build()
        .unwrap();

    let a_abs_sign: Vec<(usize, usize, u64, bool)> =  pool.install(
        || {
            mat_a.data.par_iter().map(|&(row, col, val)| (
                    row, col, val.abs() as u64, val < 0,
                )
            ).collect()
        }
    );

    let result_mutexes: Vec<Mutex<G2Element>> = (0..n_row).map(
        |_| Mutex::new(G2Element::zero())
    ).collect();
    
    
    for j_bit in 0..64{
        let v_base_arc = Arc::new(v_base_mut.clone());
        
        pool.install(|| {
            a_abs_sign.par_iter().for_each(|&(row, col, val, sign)| {
                let v_base_clone = Arc::clone(&v_base_arc);
                let result_mutex = &result_mutexes[row];
                let mut result = result_mutex.lock().unwrap();
                if (val >> j_bit) & 1 == 1 {
                    if sign {
                        *result += - v_base_clone[col];
                    } else {
                        *result += v_base_clone[col];
                    }
                }
            });
        });

        pool.install(|| {
            v_base_mut.par_iter_mut().for_each(|v| {
                *v = v.double();
            });
        });

    }


    result_mutexes.into_iter().map(
        |mutex| mutex.into_inner().unwrap()
    ).collect()
    
}



pub fn bra_opt_i128(mat_a: &Mat<i128>, v_base: &Vec<G1Element>) -> Vec<G1Element> {

    let n_row = mat_a.shape.0;
    let n_col = mat_a.shape.1;
    let mut v_base_mut = v_base[..n_row].to_vec();
       
    let pool = ThreadPoolBuilder::new()
        .num_threads(NUM_THREADS)
        .build()
        .unwrap();

    let a_abs_sign: Vec<(usize, usize, u128, bool)> =  pool.install(
            || {
                mat_a.data.par_iter().map(|&(row, col, val)| 
                    (
                        row, col, val.abs() as u128, val < 0,
                    )
                ).collect()
            }
        );

    let result_mutexes: Vec<Mutex<G1Element>> = (0..n_col).map(
        |_| Mutex::new(G1Element::zero())
        ).collect();
    
    
    for j_bit in 0..128{
        let v_base_arc = Arc::new(v_base_mut.clone());
        
        pool.install(|| {
            a_abs_sign.par_iter().for_each(|&(row, col, val, sign)| {
                let v_base_clone = Arc::clone(&v_base_arc);
                let result_mutex = &result_mutexes[col];
                let mut result = result_mutex.lock().unwrap();
                if (val >> j_bit) & 1 == 1 {
                    if sign {
                        *result += - v_base_clone[row];
                    } else {
                        *result += v_base_clone[row];
                    }
                }
            });
        });

        pool.install(|| {
            v_base_mut.par_iter_mut().for_each(|v| {
                *v = v.double();
            });
        });

    }


    result_mutexes.into_iter().map(
        |mutex| mutex.into_inner().unwrap()
    ).collect()
    
}


pub fn ket_opt_i128(mat_a: &Mat<i128>, v_base: &Vec<G2Element>) -> Vec<G2Element> {

    let n_row = mat_a.shape.0;
    let n_col = mat_a.shape.1;
    let mut v_base_mut = v_base[..n_col].to_vec();
       
    let pool = ThreadPoolBuilder::new()
        .num_threads(NUM_THREADS)
        .build()
        .unwrap();

    let a_abs_sign: Vec<(usize, usize, u128, bool)> =  pool.install(
            || {
                mat_a.data.par_iter().map(|&(row, col, val)| (
                        row, col, val.abs() as u128, val < 0,
                )).collect()
            }
    );

    let result_mutexes: Vec<Mutex<G2Element>> = (0..n_row).map(
        |_| Mutex::new(G2Element::zero())
    ).collect();
    
    
    for j_bit in 0..128{
        let v_base_arc = Arc::new(v_base_mut.clone());
        
        pool.install(|| {
            a_abs_sign.par_iter().for_each(|&(row, col, val, sign)| {
                let v_base_clone = Arc::clone(&v_base_arc);
                let result_mutex = &result_mutexes[row];
                let mut result = result_mutex.lock().unwrap();
                if (val >> j_bit) & 1 == 1 {
                    if sign {
                        *result += - v_base_clone[col];
                    } else {
                        *result += v_base_clone[col];
                    }
                }
            });
        });

        pool.install(|| {
            v_base_mut.par_iter_mut().for_each(|v| {
                *v = v.double();
            });
        });

    }


    result_mutexes.into_iter().map(
        |mutex| mutex.into_inner().unwrap()
    ).collect()
    
}





#[cfg(test)]
mod tests {
    use super::*;

    use std::time::Instant;
    use rand::Rng;

    use crate::setup::SRS;
    use crate::utils::curve::{G1Element, G2Element, GtElement};


    use crate::experiment_data;
    use crate::utils::test_data::*;

    use crate::config::SQRT_MATRIX_DIM_TEST;

    #[test]
    fn test_dirac() {
        // ga = a left-muliplied by g 
        let mat_a: Mat<u64> = gen_mat_a_u64_direct_test();
        let vec_g: Vec<G1Element> = gen_vec_v_g1_direct_test();
        let vec_ga: Vec<G1Element> = gen_vec_va_g1_from_kronecker_test();
        let vec_h: Vec<G2Element> = gen_vec_v_g2_direct_test();
        let gah: GtElement = gen_vav_gt_direct_test();

        let proj_test_left: Vec<G1Element> = proj_left(&mat_a, &vec_g);
        assert_eq!(proj_test_left, vec_ga);
        let gah_test_1 = inner_product(&proj_test_left, &vec_h);
        let bra_test = mat_a.bra(&vec_g);
        assert_eq!(bra_test, vec_ga);

        let bra_opt_test = bra_opt_u64(&mat_a, &vec_g);
        assert_eq!(bra_opt_test, vec_ga);

        let proj_test_right: Vec<G2Element> = proj_right(&mat_a, &vec_h);
        let gah_test_2 = inner_product(&vec_g, &proj_test_right);
        let ket_test = mat_a.ket(&vec_h);
        assert_eq!(ket_test, proj_test_right);
        
        let gah_test_3 = dirac(&vec_g, &mat_a, &vec_h);
        let gah_test_4 = mat_a.braket(&vec_g, &vec_h);

        assert_eq!(gah_test_1, gah);
        assert_eq!(gah_test_2, gah);
        assert_eq!(gah_test_3, gah);
        assert_eq!(gah_test_4, gah);

    }

    #[test]
    fn test_opt(){

        let srs = SRS::new(2* SQRT_MATRIX_DIM_TEST * SQRT_MATRIX_DIM_TEST);
    
        let length = SQRT_MATRIX_DIM_TEST * SQRT_MATRIX_DIM_TEST;
    
        let y = ZpElement::from(rand::thread_rng().gen::<u64>());
    
        let y_vec: Vec<ZpElement> = std::iter::successors(
            Some(y), 
            |&x| Some(x * y)
        ).take(length).collect();
    
        let (c, a, b) 
            = experiment_data::gen_matrices_sparse(SQRT_MATRIX_DIM_TEST);
    
        let c_zp = experiment_data::sprs_i128_to_zp(&c);
        let a_zp = experiment_data::sprs_i64_to_zp(&a);
        let b_zp = experiment_data::sprs_i64_to_zp(&b);
    
        let timer_bra_opt = Instant::now();
        
        let result_bra_opt = bra_opt_i64(
            &b, &srs.g_hat_vec);
        
        println!(" ** Bra i64 opt time: {:?}", timer_bra_opt.elapsed());
    
        let timer_bra_no_opt = Instant::now();
    
        let result_bra_no_opt = proj_left(
            &b_zp, &srs.g_hat_vec);
    
        println!(" ** Bra i64 no opt time: {:?}", timer_bra_no_opt.elapsed());
    
        assert_eq!(result_bra_opt, result_bra_no_opt);
        println!(" * Assert Equal between opt Bra and no opt Bra");
        
        let result_ket_no_opt = proj_right(
            &a, &srs.h_hat_vec);

        let timer_ket_opt = Instant::now();
        
        let result_ket_opt = ket_opt_i64(
            &a, &srs.h_hat_vec);

        assert_eq!(result_ket_opt, result_ket_no_opt);
        
        println!(" ** Ket i64 opt time: {:?}", timer_ket_opt.elapsed());
    
        let timer_ket_no_opt = Instant::now();
    
        let result_ket_no_opt = proj_right(
            &a_zp, &srs.h_hat_vec);

        println!(" ** Ket i64 no opt time: {:?}", timer_ket_no_opt.elapsed());
    
        assert_eq!(result_ket_opt, result_ket_no_opt);
        println!(" * Assert Equal between opt Ket and no opt Ket");
    
        let timer_bra_opt_128 = Instant::now();
        
        let result_bra_opt_128 = bra_opt_i128(
            &c, &srs.g_hat_vec);
        
        println!(" ** Bra i128 opt time: {:?}", timer_bra_opt_128.elapsed());
    
        let timer_bra_no_opt_128 = Instant::now();
    
        let result_bra_no_opt_128 = proj_left(
            &c_zp, &srs.g_hat_vec);
    
        println!(" ** Bra i128 no opt time: {:?}", timer_bra_no_opt_128.elapsed());
    
        assert_eq!(result_bra_opt_128, result_bra_no_opt_128);
        println!(" * Assert Equal of opt Bra and no opt Bra");
        
    }
}
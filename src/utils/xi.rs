//! Compute the n-vector xi by using log_2 n random challenges.
//! 
//! Compute the n-vector psi by using the n-vector xi and a random s.
//!  
//! 
use crate::utils::curve::ZpElement;
use crate::utils::dirac;
use crate::utils::curve::{Add, Mul, Zero};

/// Compute the n-vector xi by using log_2 n random challenges.
/// 
pub fn xi_from_challenges(challenges: &Vec<ZpElement>) -> Vec<ZpElement> {
    
    let log_n = challenges.len() as usize;

    let mut xi = vec![ZpElement::from(1 as u64)];

    for j in 0..log_n {
        let mut xi_left = xi.clone();
        let mut xi_right = xi.iter().map(
            |a| 
            *a * challenges[log_n - j - 1 as usize]
        ).collect::<Vec<ZpElement>>();
        xi_left.append(&mut xi_right);
        xi = xi_left;
    }

    xi 
}

/// Reduce a n-vector to a scalar by using log_2 n random challenges.
/// 
pub fn reduce_from_challenges<T>(
    challenges: &Vec<ZpElement>, vector: &Vec<T>
) -> T 
    where
    T: 'static + Clone + Copy + Send + Sync 
        + Mul<ZpElement, Output =T> + Add + Zero,
{    

    let xi = xi_from_challenges(challenges);

    if xi.len() > vector.len() {
        panic!("Vector length is too short when reducing from challenges.");
    }

    let vector_slice = vector[0..xi.len()].to_vec();

    dirac::inner_product(&vector_slice, &xi) 
}

/// Compute phi(s) in logarithmic time.
/// 
pub fn phi_s(
    s: ZpElement,
    challenges: &Vec<ZpElement>,
    exp_start: usize,
    exp_step: usize,
) -> ZpElement {
    let log_n = challenges.len() as usize;
    
    let start = s.pow(exp_start as u64);
    let step = s.pow(exp_step as u64);

    let s_pow_vec = std::iter::successors(
        Some(step),
        |&x| Some(x * x),
    ).take(log_n).collect::<Vec<ZpElement>>();


    let product = challenges.iter().rev()
        .zip(s_pow_vec.iter())
        .map(
            |(&a, &b)| 
            ZpElement::from(1 as u64)+ a * b
        ).fold(start, |acc, x| acc * x);

    product
}

/// Compute the n-vector psi by using the n-vector xi and a random s.
/// 
/// Psi is used for the verification accelleration.
/// 
pub fn psi_from_xi(xi: &Vec<ZpElement>, s: ZpElement) -> Vec<ZpElement> {
    
    let length = xi.len();

    let mut psi_rev: Vec<ZpElement> = Vec::with_capacity(length);

    psi_rev.push(xi[length-1]);

    for i in 1..length {
        psi_rev.push(psi_rev[i-1] * s + xi[length-1-i]);
    }

    psi_rev.reverse();
    
    psi_rev
}

/// Compute phi(s) directly. For unit test purpose only.
/// 
#[cfg(test)]
fn phi_s_direct(
    s: ZpElement, 
    challenges: &Vec<ZpElement>, 
    exp_start: usize, 
    exp_step: usize,
) -> ZpElement {
    
    let log_n = challenges.len() as usize;
    let n = 2usize.pow(log_n as u32);

    let start = s.pow(exp_start as u64);
    let step = s.pow(exp_step as u64);

    let s_vec: Vec<ZpElement> = std::iter::successors(
        Some(start), 
        |&x| Some(x * step)
    ).take(n).collect();

    reduce_from_challenges(&challenges, &s_vec)
}

/// Compute phi(s) by using the reduce method. For unit test purpose only.
/// 
#[cfg(test)]
fn phi_s_reduce(
    s: ZpElement, 
    challenges: &Vec<ZpElement>, 
    exp_start: usize, 
    exp_step: usize,
) -> ZpElement {
    
    let log_n = challenges.len() as usize;
    let n = 2usize.pow(log_n as u32);
 
    let start = s.pow(exp_start as u64);
    let step = s.pow(exp_step as u64);

    let s_vec: Vec<ZpElement> = std::iter::successors(
        Some(start), 
        |&x| Some(x * step)
    ).take(n).collect();

    let mut s_vec_current = s_vec.clone();

    for i in 0..log_n {
        let s_vec_left 
            = s_vec_current[0..n/2usize.pow(i as u32 + 1)].to_vec();
        let s_vec_right 
            = s_vec_current[n/2usize.pow(i as u32 + 1)..].to_vec();
        s_vec_current = s_vec_left.iter().zip(
            s_vec_right.iter()
            ).map(|(&a, &b)| 
                a + b * challenges[i]
            ).collect::<Vec<ZpElement>>();
    }

    s_vec_current[0]
}



#[cfg(test)]
mod tests{

    use super::*;

    use crate::utils::curve::G1Element;

    use rand::Rng;
    use std::time::Instant;

    fn reduce_from_challenges_g1_direct(
        challenges: &Vec<ZpElement>,
        vector: &Vec<G1Element>,
    ) -> G1Element 
    {
        let log_n = challenges.len() as usize;
        let n = 2usize.pow(log_n as u32);
    
        let mut vec_current: Vec<G1Element> = vector.clone();

        for i in 0..log_n {
            let vec_left 
                = vec_current[0..n/2usize.pow(i as u32 + 1)].to_vec();
            let vec_right 
                = vec_current[n/2usize.pow(i as u32 + 1)..].to_vec();
            vec_current = vec_left.iter().zip(
                vec_right.iter()
                ).map(|(&a, &b)| 
                    a + b * challenges[i]
                ).collect::<Vec<G1Element>>();
        }

        vec_current[0]
    }

    #[test]
    fn test_xi(){
        let log_n = 10;

        let challenges = (0..log_n).map(|_| 
            ZpElement::rand()
        ).collect::<Vec<ZpElement>>();

        let s = ZpElement::rand();

        let timer_direct = Instant::now();
        let phi_s_direct = phi_s_direct(
            s, &challenges, 3, 2);
        println!(
            " ** Compute phi_s direct time: {:?}", timer_direct.elapsed()
        );
        
        let timer_opt = Instant::now();
        let phi_s_opt = phi_s(
            s, &challenges, 3, 2);
        println!(
            " ** Compute phi_s opt time: {:?}", timer_opt.elapsed()
        );
        
        let timer_reduce = Instant::now();
        let phi_s_reduce = phi_s_reduce(
            s, &challenges, 3, 2);
        println!(
            " ** Compute phi_s reduce time: {:?}", timer_reduce.elapsed()
        );

        assert_eq!(phi_s_direct, phi_s_reduce);
        assert_eq!(phi_s_reduce, phi_s_opt);

        println!(" * Assert Equal across three methods ");
        
        let n = 2usize.pow(log_n as u32);
        let vector_g1 = (0..n).map(|_| 
            G1Element::from(
                rand::thread_rng().gen::<u64>()
            )
        ).collect::<Vec<G1Element>>();

        let timer_g1_reduce = Instant::now();
        let g1_reduce = reduce_from_challenges(
            &challenges, &vector_g1);
        println!(" ** Compute g1 reduce opt time: {:?}", 
                timer_g1_reduce.elapsed());

        let timer_g1_direct = Instant::now();
        let g1_reduce_direct = reduce_from_challenges_g1_direct(
            &challenges, &vector_g1);
        println!(" ** Compute g1 reduce direct time: {:?}", 
                timer_g1_direct.elapsed());


        assert_eq!(g1_reduce_direct, g1_reduce);

        println!(" * Assert Equal between two methods ");

        let s_hat = ZpElement::rand();

        let s_hat_vec: Vec<ZpElement> = std::iter::successors(
            Some(ZpElement::from(1 as u64)), 
            |&x| Some(x * s_hat)
        ).take(n).collect();


        let xi = xi_from_challenges(&challenges);
        let psi = psi_from_xi(&xi, s);

        let xi_s = phi_s(
            s, &challenges, 1, 1);
        let xi_s_hat = phi_s(
            s_hat, &challenges, 1, 1);
        let psi_s_hat = dirac::inner_product(
            &psi, &s_hat_vec);
        assert_eq!((s-s_hat) * psi_s_hat, xi_s - xi_s_hat );
  
        println!(" * Assert correct computation of psi ");

    }

}
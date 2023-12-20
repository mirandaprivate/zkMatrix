//! Compute the n-vector xi by using log_2 n random challenges.
//! 
//! 
//!  
//! 
use crate::utils::curve::ZpElement;

pub fn xi_from_challenges_recursive(challenges: &Vec<ZpElement>) -> Vec<ZpElement> {
    
    let log_n = challenges.len() as usize;

    let mut xi = vec![ZpElement::from(1 as u64)];

    for j in 0..log_n {
        let mut xi_left = xi.clone();
        let mut xi_right = xi.iter().map(
            |a| 
            *a * challenges[log_n -1 - j as usize]
        ).collect::<Vec<ZpElement>>();
        xi_left.append(&mut xi_right);
        xi = xi_left;
    }

    xi 
}
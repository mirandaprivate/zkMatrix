use bls12_381::{Scalar, G1Projective, G2Projective, pairing, Gt, G1Affine, G2Affine};
use rand::Rng;
use std::thread;

fn random_scalar<R: Rng + ?Sized>(rng: &mut R) -> Scalar {
    Scalar::from_raw([
        rng.gen::<u64>(),
        rng.gen::<u64>(),
        rng.gen::<u64>(),
        rng.gen::<u64>()
    ])
}

fn random_g1<R: Rng + ?Sized>(rng: &mut R) -> G1Projective {
    G1Projective::generator() * random_scalar(rng)
}

fn random_g2<R: Rng + ?Sized>(rng: &mut R) -> G2Projective {
    G2Projective::generator() * random_scalar(rng)
}

const MATRIX_SIZE: usize = 64;
const NUM_NONZERO_ENTRIES: usize = 256;
const NUM_THREADS: usize = 8;

fn main() {
    let mut rng = rand::thread_rng();
    let mut m: Vec<(usize, usize, Scalar)> = Vec::new();

    for _ in 0..NUM_NONZERO_ENTRIES {
        let i = rng.gen_range(0..MATRIX_SIZE);
        let j = rng.gen_range(0..MATRIX_SIZE);
        m.push((i, j, random_scalar(&mut rng)));
    }

    let v1: Vec<G1Projective> = (0..MATRIX_SIZE).map(|_| random_g1(&mut rng)).collect();
    let v2: Vec<G2Projective> = (0..MATRIX_SIZE).map(|_| random_g2(&mut rng)).collect();

    let chunk_size = NUM_NONZERO_ENTRIES / NUM_THREADS;
    let m_chunks: Vec<Vec<(usize, usize, Scalar)>> = m.chunks(chunk_size).map(|chunk| chunk.to_vec()).collect();

    let mut aggregate_intermediate = vec![G1Projective::identity(); MATRIX_SIZE];
    let mut handles: Vec<_> = m_chunks.into_iter().map(|chunk| {
        let v1_clone = v1.clone();
        thread::spawn(move || {
            let mut local_intermediate = vec![G1Projective::identity(); MATRIX_SIZE];
            for &(row, col, ref val) in &chunk {
                local_intermediate[row] += v1_clone[col] * val;
            }
            local_intermediate
        })
    }).collect();

    for handle in handles.drain(..) {
        let local_intermediate = handle.join().unwrap();
        for i in 0..MATRIX_SIZE {
            aggregate_intermediate[i] += local_intermediate[i];
        }
    }

    let mut result = Gt::identity();
    for i in 0..MATRIX_SIZE {
        result += pairing(&G1Affine::from(aggregate_intermediate[i]), &G2Affine::from(v2[i]));
    }

    println!("Result: {:?}", result);
}

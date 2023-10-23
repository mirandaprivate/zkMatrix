# zkMatrix: Conquering Zero-Knowledge Proof for Large Matrix Multiplication

This folder contains codes for our paper
*Pairing-Based zkMatrix: Conquering Zero-Knowledge Proof for Large Matrix Multiplication*.

## Overview

Given a pairing
$e: \mathbb{G}_1 \times \mathbb{G}_2 \mapsto \mathbb{G}_T$, 
and two vectors 
$\vec{\mathbf{G}} \in \mathbb{G}_1^{\star}$ 
and 
$\vec{\mathbf{H}} \in \mathbb{G}_2^{\star}$ ,

then the two-tier commitment of a $m \times n$ matrix 
$C_a = \mathbf{a} = \{a_{ij}\} \in \mathbb{Z}_p^{m\times n}$ is defined by:

$$
\langle \vec{\mathbf{G}}  |  \mathbf{a}   |  \vec{\mathbf{H}} \rangle
: = \bigoplus_{i=1}^m \bigoplus_{j=1}^n a_{ij} e(G_i, H_j)
$$ 

Suppose the prover has made commitments to three $m \times n$ matrix 
$\mathbf{a}$, $\mathbf{b}$, and $\mathbf{c}$ as follows:
$$ 
C_a = \langle \vec{\mathbf{G}}  |  \mathbf{a}   |  \vec{\mathbf{H}} \rangle 
\in \mathbb{G}_T, 
\\\\
C_b =  \langle \vec{\mathbf{G}}  |  \mathbf{b}   |  \vec{\mathbf{H}} \rangle 
\in \mathbb{G}_T,
\\\\
C_c =  \langle \vec{\mathbf{G}}  |  \mathbf{c}  |  \vec{\mathbf{H}} \rangle
\in \mathbb{G}_T .
$$

Then, the prover can generate a zero-knowledge proof with $O(n)$ time complexity
for the relation:
$$
\mathcal{R} 
= \\{ 
     C_c \in \mathbb{G}_T, C_a \in \mathbb{G}_T, C_b \in \mathbb{G}_T;
    \vec{\mathbf{G}} \in \mathbb{G}_1^{\star} , \vec{\mathbf{H}} \in \mathbb{G}_2^{\star} 
\\\\
: \mathbf{a} \in \mathbb{Z}_p^{m\times l},
    \mathbf{b} \in \mathbb{Z}_p^{l \times n},
    \mathbf{c} \in \mathbb{Z}_p^{m \times n}
    \\\\
| \mathbf{c} = \mathbf{a} \mathbf{b} 
    \wedge C_c =
     \langle \vec{\mathbf{G}}  |  \mathbf{c}   |  \vec{\mathbf{H}} \rangle
    \wedge C_a =
     \langle \vec{\mathbf{G}}  |  \mathbf{a}   |  \vec{\mathbf{H}} \rangle
    \wedge C_b =
     \langle \vec{\mathbf{G}}  |  \mathbf{b}   |  \vec{\mathbf{H}} \rangle     
\\}.
$$

We employ the random oracle approach.

![alg](assets/alg7.png)

---

## Introduction 



## Getting Started

### Running the Code

To execute our program, run the following command:
```bash
cd /path/to/zkmatrix
cargo run
```

## Compatibility Note

- **Rust Toolchain:** 3.10.12
- **Environment:** Ubuntu 22.04

---

## Algorithms

### Matrix Commitment

### Inner-Product Argument

### Semi-Inner-Product Argument

### High-Dimensional Semi-Inner-Product Argument

### Protocol MatMul

---

## Directory Contents

- **lib.rs:** The primary entry point for the program.
- **util/** Utility functions.
- **op/** Zero-knowledge proofs for various matrix oprators.

For more

--- 

## Citing

If our work benefits to your research, please cite our paper as follows:

This is math in README file:

$$ \sqrt{\mathbf{x}}$$
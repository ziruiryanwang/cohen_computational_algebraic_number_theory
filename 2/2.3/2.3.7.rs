use num_rational::BigRational;
use num_traits::Zero;

use crate::{
    inverse_image_matrix,
    matrix::{mat_mul_rect},
    supplement_basis,
    Matrix,
};

pub fn supplement_subspace(v: Matrix, m: Matrix) -> Option<Matrix> {
    if v.len() != m.len() {
        return None;
    }
    let n = m.first().map(|row| row.len()).unwrap_or(0);
    let r = v.first().map(|row| row.len()).unwrap_or(0);
    if r > n {
        return None;
    }
    if m.iter().any(|row| row.len() != n) || v.iter().any(|row| row.len() != r) {
        return None;
    }

    let x = inverse_image_matrix(m.clone(), v)?;
    let b = supplement_basis(x)?;

    let extra = n - r;
    let mut c = vec![vec![BigRational::zero(); extra]; n];
    for i in 0..n {
        for j in 0..extra {
            c[i][j] = b[i][r + j].clone();
        }
    }

    mat_mul_rect(&m, &c)
}

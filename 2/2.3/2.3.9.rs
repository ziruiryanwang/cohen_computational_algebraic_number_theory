use num_rational::BigRational;
use num_traits::Zero;

use crate::{
    image_basis, kernel_basis,
    matrix::{columns_to_matrix, mat_mul_rect},
    Matrix,
};

pub fn intersection_subspaces(m: Matrix, m_prime: Matrix) -> Option<Matrix> {
    if m.len() != m_prime.len() {
        return None;
    }
    let rows = m.len();
    let n = m.first().map(|row| row.len()).unwrap_or(0);
    let n_prime = m_prime.first().map(|row| row.len()).unwrap_or(0);
    if m.iter().any(|row| row.len() != n) || m_prime.iter().any(|row| row.len() != n_prime) {
        return None;
    }

    let mut m1 = Vec::with_capacity(rows);
    for i in 0..rows {
        let mut row = m[i].clone();
        row.extend_from_slice(&m_prime[i]);
        m1.push(row);
    }

    let ker = kernel_basis(m1);
    let p = ker.len();
    let total_cols = n + n_prime;
    let mut n_mat = vec![vec![BigRational::zero(); p]; total_cols];
    for (col_idx, v) in ker.iter().enumerate() {
        for row_idx in 0..total_cols {
            n_mat[row_idx][col_idx] = v[row_idx].clone();
        }
    }

    let n1: Matrix = n_mat.into_iter().take(n).collect();
    let m2 = mat_mul_rect(&m, &n1)?;
    let basis = image_basis(m2);
    if basis.is_empty() {
        return Some(vec![vec![]; rows]);
    }
    Some(columns_to_matrix(&basis))
}

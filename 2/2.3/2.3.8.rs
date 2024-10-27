use crate::{image_basis, matrix::columns_to_matrix, Matrix};

pub fn sum_subspaces(m: Matrix, m_prime: Matrix) -> Option<Matrix> {
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

    let basis = image_basis(m1);
    if basis.is_empty() {
        return Some(vec![vec![]; rows]);
    }
    Some(columns_to_matrix(&basis))
}

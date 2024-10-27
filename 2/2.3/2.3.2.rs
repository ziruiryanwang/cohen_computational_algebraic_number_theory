use num_traits::Zero;

use crate::matrix::{Matrix, Vector};

pub fn image_basis(m: Matrix) -> Vec<Vector> {
    let (_rref, pivots) = rref_with_pivots(m.clone());
    let mut basis = Vec::new();
    for &pivot_col in pivots.iter() {
        let col_vec: Vector = m.iter().map(|row| row[pivot_col].clone()).collect();
        basis.push(col_vec);
    }
    basis
}

fn rref_with_pivots(mut a: Matrix) -> (Matrix, Vec<usize>) {
    let m_rows = a.len();
    if m_rows == 0 {
        return (a, Vec::new());
    }
    let n_cols = a[0].len();
    let mut row = 0usize;
    let mut pivots = Vec::new();
    for col in 0..n_cols {
        if row >= m_rows {
            break;
        }
        let mut pivot = None;
        for r in row..m_rows {
            if !a[r][col].is_zero() {
                pivot = Some(r);
                break;
            }
        }
        let pivot_row = match pivot {
            Some(r) => r,
            None => continue,
        };
        if pivot_row != row {
            a.swap(pivot_row, row);
        }
        let inv = a[row][col].clone().recip();
        for c in col..n_cols {
            a[row][c] *= inv.clone();
        }
        for r in 0..m_rows {
            if r == row {
                continue;
            }
            if a[r][col].is_zero() {
                continue;
            }
            let factor = a[r][col].clone();
            for c in col..n_cols {
                a[r][c] = a[r][c].clone() - &factor * &a[row][c];
            }
        }
        pivots.push(col);
        row += 1;
    }
    (a, pivots)
}

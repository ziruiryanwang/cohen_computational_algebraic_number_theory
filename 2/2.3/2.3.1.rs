use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::matrix::Matrix;

pub fn kernel_basis(m: Matrix) -> Vec<Vec<BigRational>> {
    let (rref, pivots) = rref_with_pivots(m);
    let n_cols = if let Some(row) = rref.first() { row.len() } else { 0 };
    let pivot_set: std::collections::HashSet<usize> = pivots.iter().cloned().collect();
    let mut basis = Vec::new();
    for free_col in 0..n_cols {
        if pivot_set.contains(&free_col) {
            continue;
        }
        let mut v = vec![BigRational::zero(); n_cols];
        v[free_col] = BigRational::one();
        for (row, &pivot_col) in pivots.iter().enumerate() {
            v[pivot_col] = -rref[row][free_col].clone();
        }
        basis.push(v);
    }
    basis
}

fn rref_with_pivots(mut a: Matrix) -> (Matrix, Vec<usize>) {
    let m = a.len();
    if m == 0 {
        return (a, Vec::new());
    }
    let n = a[0].len();
    let mut row = 0;
    let mut pivots = Vec::new();
    for col in 0..n {
        if row >= m {
            break;
        }
        let mut pivot = None;
        for r in row..m {
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
        for c in col..n {
            a[row][c] *= inv.clone();
        }
        for r in 0..m {
            if r == row {
                continue;
            }
            if a[r][col].is_zero() {
                continue;
            }
            let factor = a[r][col].clone();
            for c in col..n {
                a[r][c] = a[r][c].clone() - &factor * &a[row][c];
            }
        }
        pivots.push(col);
        row += 1;
    }
    (a, pivots)
}

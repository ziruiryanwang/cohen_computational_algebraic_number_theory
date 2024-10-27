use num_bigint::BigInt;

use crate::{lll_dependent, lll_integral};

fn to_columns(a: &[Vec<BigInt>]) -> Vec<Vec<BigInt>> {
    if a.is_empty() {
        return Vec::new();
    }
    let m = a.len();
    let n = a[0].len();
    let mut cols = vec![vec![BigInt::from(0); m]; n];
    for i in 0..m {
        for j in 0..n {
            cols[j][i] = a[i][j].clone();
        }
    }
    cols
}

fn columns_to_rows(cols: &[Vec<BigInt>]) -> Vec<Vec<BigInt>> {
    if cols.is_empty() {
        return Vec::new();
    }
    let n = cols.len();
    let m = cols[0].len();
    let mut rows = vec![vec![BigInt::from(0); n]; m];
    for j in 0..n {
        for i in 0..m {
            rows[i][j] = cols[j][i].clone();
        }
    }
    rows
}

pub fn kernel_lll(a: Vec<Vec<BigInt>>) -> Option<Vec<Vec<BigInt>>> {
    let m = a.len();
    if m == 0 {
        return Some(Vec::new());
    }
    let n = a[0].len();
    if a.iter().any(|row| row.len() != n) {
        return None;
    }

    let cols = to_columns(&a);
    let (_basis, h, p) = lll_dependent(cols)?;
    let r = n - p;
    if r == 0 {
        return Some(Vec::new());
    }
    let h_cols = to_columns(&h);
    let kernel_cols: Vec<Vec<BigInt>> = h_cols.into_iter().take(r).collect();
    let (reduced_kernel, _) = lll_integral(kernel_cols)?;
    Some(columns_to_rows(&reduced_kernel))
}

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::Zero;

use crate::{
    inverse,
    lll_dependent,
    lll_integral,
    Matrix,
};

fn to_columns(a: Vec<Vec<BigInt>>) -> Vec<Vec<BigInt>> {
    if a.is_empty() {
        return Vec::new();
    }
    let m = a.len();
    let n = a[0].len();
    let mut cols = vec![vec![BigInt::zero(); m]; n];
    for i in 0..m {
        for j in 0..n {
            cols[j][i] = a[i][j].clone();
        }
    }
    cols
}

fn columns_to_rows(cols: Vec<Vec<BigInt>>) -> Vec<Vec<BigInt>> {
    if cols.is_empty() {
        return Vec::new();
    }
    let n = cols.len();
    let m = cols[0].len();
    let mut rows = vec![vec![BigInt::zero(); n]; m];
    for j in 0..n {
        for i in 0..m {
            rows[i][j] = cols[j][i].clone();
        }
    }
    rows
}

fn dot(a: &[BigInt], b: &[BigInt]) -> BigInt {
    let mut s = BigInt::zero();
    for (x, y) in a.iter().zip(b.iter()) {
        s += x * y;
    }
    s
}

pub fn kernel_image_lll(a: Vec<Vec<BigInt>>) -> Option<(Vec<Vec<BigInt>>, usize)> {
    let m = a.len();
    if m == 0 {
        return None;
    }
    let n = a[0].len();
    if a.iter().any(|row| row.len() != n) {
        return None;
    }

    let cols = to_columns(a);
    let (_b_lll, mut h, p) = lll_dependent(cols)?;
    let r = n - p;

    if r > 0 {
        let h_cols = to_columns(h.clone());
        let kernel_basis: Vec<Vec<BigInt>> = h_cols.iter().take(r).cloned().collect();
        if let Some((reduced_kernel, _)) = lll_integral(kernel_basis) {
            let kernel_rows = columns_to_rows(reduced_kernel);
            for i in 0..h.len() {
                for j in 0..r {
                    h[i][j] = kernel_rows[i][j].clone();
                }
            }
        }
    }

    if r == 0 {
        return Some((h, p));
    }

    let mut gram: Matrix = vec![vec![BigRational::zero(); r]; r];
    let h_cols = to_columns(h.clone());
    for j in 0..r {
        for k in 0..r {
            gram[j][k] = BigRational::from_integer(dot(&h_cols[j], &h_cols[k]));
        }
    }
    let d_inv = inverse(gram)?;

    let h_cols_mut = h_cols;
    let mut h_rows = columns_to_rows(h_cols_mut.clone());

    for i_idx in r..n {
        let hi = &h_cols_mut[i_idx];
        let mut v = vec![BigRational::zero(); r];
        for j in 0..r {
            v[j] = BigRational::from_integer(dot(hi, &h_cols_mut[j]));
        }
        let mut dv = vec![BigRational::zero(); r];
        for row in 0..r {
            for col in 0..r {
                dv[row] += &d_inv[row][col] * &v[col];
            }
        }

        let mut mjs = vec![BigInt::zero(); r];
        for j in 0..r {
            let num = dv[j].numer().clone();
            let den = dv[j].denom().clone();
            mjs[j] = num.div_floor(&den);
        }

        for row in 0..m {
            let mut adjustment = BigInt::zero();
            for j in 0..r {
                adjustment += &mjs[j] * &h_rows[row][j];
            }
            h_rows[row][i_idx] = h_rows[row][i_idx].clone() - adjustment;
        }
    }

    let h_final = h_rows;
    Some((h_final, p))
}

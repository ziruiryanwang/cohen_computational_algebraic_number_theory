use num_bigint::BigInt;
use num_traits::{ToPrimitive, Zero};

use crate::{
    cholesky_decomposition, lll_reduction, short_vectors,
};

fn upper_tri_inverse(r: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = r.len();
    let mut inv_cols = vec![vec![0f64; n]; n];
    for j in 0..n {
        let mut x = vec![0f64; n];
        for i_rev in 0..n {
            let i = n - 1 - i_rev;
            let mut sum = 0.0;
            for k in (i + 1)..n {
                sum += r[i][k] * x[k];
            }
            let rhs = if i == j { 1.0 } else { 0.0 };
            if r[i][i] == 0.0 {
                return None;
            }
            x[i] = (rhs - sum) / r[i][i];
        }
        for i in 0..n {
            inv_cols[j][i] = x[i];
        }
    }
    let mut inv = vec![vec![0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            inv[i][j] = inv_cols[j][i];
        }
    }
    Some(inv)
}

fn mat_mul_f64(a: &[Vec<f64>], b: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    if a.is_empty() || b.is_empty() {
        return Some(Vec::new());
    }
    let a_rows = a.len();
    let a_cols = a[0].len();
    if b.len() != a_cols {
        return None;
    }
    let b_cols = b[0].len();
    let mut res = vec![vec![0f64; b_cols]; a_rows];
    for i in 0..a_rows {
        for k in 0..a_cols {
            for j in 0..b_cols {
                res[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    Some(res)
}

fn transpose_f64(m: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if m.is_empty() {
        return Vec::new();
    }
    let r = m.len();
    let c = m[0].len();
    let mut t = vec![vec![0f64; r]; c];
    for i in 0..r {
        for j in 0..c {
            t[j][i] = m[i][j];
        }
    }
    t
}

fn col_norm2(m: &[Vec<f64>], j: usize) -> f64 {
    let mut s = 0.0;
    for i in 0..m.len() {
        s += m[i][j] * m[i][j];
    }
    s
}

fn apply_permutation_cols(m: &[Vec<f64>], perm: &[usize]) -> Vec<Vec<f64>> {
    let rows = m.len();
    let cols = perm.len();
    let mut res = vec![vec![0f64; cols]; rows];
    for i in 0..rows {
        for (new_j, &old_j) in perm.iter().enumerate() {
            res[i][new_j] = m[i][old_j];
        }
    }
    res
}

fn big_mat_vec_mul(h: &[Vec<BigInt>], v: &[BigInt]) -> Vec<BigInt> {
    let rows = h.len();
    let cols = if rows == 0 { 0 } else { h[0].len() };
    let mut res = vec![BigInt::zero(); rows];
    for i in 0..rows {
        let mut s = BigInt::zero();
        for j in 0..cols {
            s += &h[i][j] * &v[j];
        }
        res[i] = s;
    }
    res
}

pub fn fincke_pohst(
    a: Vec<Vec<f64>>,
    c: f64,
    scale: i64,
) -> Option<Vec<(Vec<BigInt>, f64)>> {
    let n = a.len();
    if n == 0 || a.iter().any(|row| row.len() != n) || scale == 0 {
        return None;
    }
    let (_q_tmp, r) = cholesky_decomposition(a)?;
    let r_inv = upper_tri_inverse(&r)?;
    let mut rows_int = Vec::with_capacity(n);
    for row in r_inv.iter() {
        let mut v = Vec::with_capacity(n);
        for &val in row.iter() {
            let scaled = (val * scale as f64).round() as i64;
            v.push(BigInt::from(scaled));
        }
        rows_int.push(v);
    }

    let (_b_red, h) = lll_reduction(rows_int)?;

    let mut h_f64 = vec![vec![0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            h_f64[i][j] = h[i][j].to_f64().unwrap_or(0.0) / scale as f64;
        }
    }
    let s = mat_mul_f64(&r, &h_f64)?;

    let mut norms: Vec<(usize, f64)> = (0..n).map(|j| (j, col_norm2(&s, j))).collect();
    norms.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    let perm: Vec<usize> = norms.iter().map(|&(idx, _)| idx).collect();
    let s_perm = apply_permutation_cols(&s, &perm);

    let s_t = transpose_f64(&s_perm);
    let a1 = mat_mul_f64(&s_t, &s_perm)?;
    let (q1, _r1) = cholesky_decomposition(a1)?;

    let sv = short_vectors(q1, c);
    let mut res = Vec::new();
    for (y, qval) in sv {
        let x = big_mat_vec_mul(&h, &y);
        res.push((x, qval));
    }
    Some(res)
}

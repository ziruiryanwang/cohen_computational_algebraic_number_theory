use num_rational::BigRational;
use num_traits::Zero;

use crate::{Matrix, Vector};

pub fn inverse_image_matrix(m: Matrix, v: Matrix) -> Option<Matrix> {
    let m_rows = m.len();
    if m_rows == 0 || v.len() != m_rows || m[0].len() < v[0].len() {
        return None;
    }
    let n = m[0].len();
    let r = v[0].len();

    // copy M to work on, and B <- V
    let mut m_work = m.clone();
    let mut b_work = v.clone();

    let mut j: isize = -1;
    while {
        j += 1;
        (j as usize) < n
    } {
        let col = j as usize;
        let mut pivot = None;
        for i in col..m_rows {
            if !m_work[i][col].is_zero() {
                pivot = Some(i);
                break;
            }
        }
        let i = match pivot {
            Some(v) => v,
            None => return None,
        };

        if i > col {
            m_work.swap(i, col);
            b_work.swap(i, col);
        }

        let mjj_inv = m_work[col][col].clone().recip();
        let mut c: Vec<BigRational> = vec![BigRational::zero(); m_rows];
        for k in (col + 1)..m_rows {
            c[k] = &mjj_inv * &m_work[k][col];
        }

        for k in (col + 1)..m_rows {
            for l in (col + 1)..n {
                m_work[k][l] = &m_work[k][l] - &c[k] * &m_work[col][l];
            }
            for l in 0..r {
                b_work[k][l] = &b_work[k][l] - &c[k] * &b_work[col][l];
            }
        }
    }

    // Back substitution on first n rows (upper triangular)
    let mut x = vec![vec![BigRational::zero(); r]; n];
    for i_rev in 0..n {
        let i = n - 1 - i_rev;
        for col_r in 0..r {
            let mut sum = BigRational::zero();
            for j in (i + 1)..n {
                sum += &m_work[i][j] * &x[j][col_r];
            }
            if m_work[i][i].is_zero() {
                return None;
            }
            x[i][col_r] = (&b_work[i][col_r] - sum) / m_work[i][i].clone();
        }
    }

    // Check rest of rows
    for k in n..m_rows {
        for col_r in 0..r {
            let mut sum = BigRational::zero();
            for j in 0..n {
                sum += &m[k][j] * &x[j][col_r];
            }
            if sum != b_work[k][col_r] {
                return None;
            }
        }
    }

    Some(x)
}

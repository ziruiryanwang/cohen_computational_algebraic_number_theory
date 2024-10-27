use num_traits::Zero;

use crate::Matrix;

pub fn column_echelon_form(mut m: Matrix) -> Matrix {
    let rows = m.len();
    let n_cols = m.first().map(|row| row.len()).unwrap_or(0);
    if m.iter().any(|row| row.len() != n_cols) {
        return m;
    }
    if rows == 0 {
        return m;
    }

    let mut i: isize = rows as isize - 1;
    let mut k: isize = n_cols as isize - 1;

    while i >= 0 {
        if k < 0 {
            if i == 0 {
                break;
            }
            i -= 1;
            continue;
        }

        let mut pivot_j: Option<isize> = None;
        let mut j = k;
        while j >= 0 {
            if !m[i as usize][j as usize].is_zero() {
                pivot_j = Some(j);
                break;
            }
            j -= 1;
        }

        let pivot_col = match pivot_j {
            Some(v) => v as usize,
            None => {
                if i == 0 {
                    break;
                }
                i -= 1;
                continue;
            }
        };
        let k_idx = k as usize;
        let d = m[i as usize][pivot_col].clone().recip();

        for l in 0..=i as usize {
            let t = &d * &m[l][pivot_col];
            if pivot_col != k_idx {
                m[l][pivot_col] = m[l][k_idx].clone();
            }
            m[l][k_idx] = t;
        }

        let pivot_row = m[i as usize].clone();
        for j2 in 0..n_cols {
            if j2 == k_idx {
                continue;
            }
            let factor = pivot_row[j2].clone();
            if factor.is_zero() {
                continue;
            }
            for l in 0..=i as usize {
                let val = m[l][j2].clone() - m[l][k_idx].clone() * factor.clone();
                m[l][j2] = val;
            }
        }

        k -= 1;
        if i == 0 {
            break;
        }
        i -= 1;
    }

    m
}

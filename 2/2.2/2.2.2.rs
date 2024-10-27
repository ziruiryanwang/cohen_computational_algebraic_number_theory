use num_rational::BigRational;
use num_traits::Zero;

use crate::matrix::{identity, Matrix};

pub fn inverse(mut m: Matrix) -> Option<Matrix> {
    let n = m.len();
    if m.iter().any(|row| row.len() != n) {
        return None;
    }
    let mut b = identity(n);

    let mut j: isize = -1;
    while {
        j += 1;
        (j as usize) < n
    } {
        let col = j as usize;
        let mut pivot = None;
        for i in col..n {
            if !m[i][col].is_zero() {
                pivot = Some(i);
                break;
            }
        }
        let i = pivot?;

        if i > col {
            m.swap(i, col);
            b.swap(i, col);
        }

        let mjj = m[col][col].clone();
        let d = mjj.recip();
        let mut c: Vec<BigRational> = vec![BigRational::zero(); n];
        for k in (col + 1)..n {
            c[k] = &d * &m[k][col];
        }

        for k in (col + 1)..n {
            for l in (col + 1)..n {
                m[k][l] = &m[k][l] - &c[k] * &m[col][l];
            }
            for l in 0..n {
                b[k][l] = &b[k][l] - &c[k] * &b[col][l];
            }
        }
    }

    let mut x = vec![vec![BigRational::zero(); n]; n];
    for i_rev in 0..n {
        let i = n - 1 - i_rev;
        for r in 0..n {
            let mut sum = BigRational::zero();
            for j in (i + 1)..n {
                sum += &m[i][j] * &x[j][r];
            }
            if m[i][i].is_zero() {
                return None;
            }
            x[i][r] = (&b[i][r] - sum) / m[i][i].clone();
        }
    }
    Some(x)
}

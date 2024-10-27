use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::matrix::Matrix;

pub fn determinant(mut m: Matrix) -> Option<BigRational> {
    let n = m.len();
    if m.iter().any(|row| row.len() != n) {
        return None;
    }
    let mut det = BigRational::one();
    let mut sign = BigRational::one();
    for j in 0..n {
        let mut pivot = None;
        for i in j..n {
            if !m[i][j].is_zero() {
                pivot = Some(i);
                break;
            }
        }
        let i = match pivot {
            Some(v) => v,
            None => return Some(BigRational::zero()),
        };
        if i > j {
            m.swap(i, j);
            sign = -sign;
        }
        let mjj = m[j][j].clone();
        det *= &mjj;
        let d = mjj.recip();
        for k in (j + 1)..n {
            let c = &d * &m[k][j];
            for l in (j + 1)..n {
                m[k][l] = &m[k][l] - &c * &m[j][l];
            }
        }
    }
    Some(det * sign)
}

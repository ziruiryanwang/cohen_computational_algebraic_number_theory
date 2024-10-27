use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

pub fn gauss_bareiss_det(mut m: Vec<Vec<BigInt>>) -> Option<BigInt> {
    let n = m.len();
    if n == 0 || m.iter().any(|row| row.len() != n) {
        return None;
    }
    let mut c = BigInt::one();
    let mut sign = BigInt::one();
    for k in 0..n {
        if k == n - 1 {
            return Some(sign * m[k][k].clone());
        }
        if m[k][k].is_zero() {
            let mut pivot = None;
            for i in (k + 1)..n {
                if !m[i][k].is_zero() {
                    pivot = Some(i);
                    break;
                }
            }
            if let Some(i) = pivot {
                m.swap(k, i);
                sign = -sign;
            } else {
                return Some(BigInt::zero());
            }
        }
        let p = m[k][k].clone();
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let t = &p * &m[i][j] - &m[i][k] * &m[k][j];
                if c.is_one() {
                    m[i][j] = t;
                } else {
                    let (q, r) = t.div_rem(&c);
                    if !r.is_zero() {
                        return None;
                    }
                    m[i][j] = q;
                }
            }
        }
        c = p;
    }
    None
}

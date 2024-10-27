use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

use crate::extended_euclid;

fn identity_bigint(n: usize) -> Vec<Vec<BigInt>> {
    let mut m = vec![vec![BigInt::zero(); n]; n];
    for i in 0..n {
        m[i][i] = BigInt::one();
    }
    m
}

pub fn kernel_integer(mut a: Vec<Vec<BigInt>>) -> Option<Vec<Vec<BigInt>>> {
    let m = a.len();
    if m == 0 {
        return Some(Vec::new());
    }
    let n = a[0].len();
    if a.iter().any(|row| row.len() != n) {
        return None;
    }

    let mut u = identity_bigint(n);
    let l = if m <= n { 0isize } else { (m - n) as isize };
    let mut i = m as isize - 1;
    let mut k = n as isize - 1;
    let mut j = k;

    loop {
        if j <= 0 {
            let mut b = a[i as usize][k as usize].clone();
            if b.is_zero() {
                return None;
            }
            if b.is_negative() {
                for row in 0..m {
                    a[row][k as usize] = -a[row][k as usize].clone();
                    u[row][k as usize] = -u[row][k as usize].clone();
                }
                b = -b;
            }

            for col in (k as usize + 1)..n {
                let q = a[i as usize][col].div_floor(&b);
                if q.is_zero() {
                    continue;
                }
                for row in 0..m {
                    a[row][col] = a[row][col].clone() - &q * &a[row][k as usize];
                }
                for row in 0..n {
                    u[row][col] = u[row][col].clone() - &q * &u[row][k as usize];
                }
            }

            if i as usize == l as usize {
                let basis_cols = k as usize;
                let mut m_out = vec![vec![BigInt::zero(); basis_cols]; n];
                for row in 0..n {
                    for col in 0..basis_cols {
                        m_out[row][col] = u[row][col].clone();
                    }
                }
                return Some(m_out);
            } else {
                i -= 1;
                k -= 1;
                j = k;
                continue;
            }
        }

        j -= 1;
        if a[i as usize][j as usize].is_zero() {
            continue;
        }

        let aik = a[i as usize][k as usize].clone();
        let aij = a[i as usize][j as usize].clone();
        let (u_c, v_c, d) = extended_euclid(aik.clone(), aij.clone());
        if d.is_zero() {
            return None;
        }
        let factor_k = aik.div_floor(&d);
        let factor_j = aij.div_floor(&d);

        let mut b_col_a = vec![BigInt::zero(); m];
        let mut new_a_j = vec![BigInt::zero(); m];
        for row in 0..m {
            b_col_a[row] = &u_c * a[row][k as usize].clone() + &v_c * a[row][j as usize].clone();
            new_a_j[row] = factor_k.clone() * a[row][j as usize].clone()
                - factor_j.clone() * a[row][k as usize].clone();
        }
        for row in 0..m {
            a[row][k as usize] = b_col_a[row].clone();
            a[row][j as usize] = new_a_j[row].clone();
        }

        let mut b_col_u = vec![BigInt::zero(); n];
        let mut new_u_j = vec![BigInt::zero(); n];
        for row in 0..n {
            b_col_u[row] = &u_c * u[row][k as usize].clone() + &v_c * u[row][j as usize].clone();
            new_u_j[row] = factor_k.clone() * u[row][j as usize].clone()
                - factor_j.clone() * u[row][k as usize].clone();
        }
        for row in 0..n {
            u[row][k as usize] = b_col_u[row].clone();
            u[row][j as usize] = new_u_j[row].clone();
        }

        continue;
    }
}

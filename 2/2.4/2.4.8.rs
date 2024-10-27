use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::Zero;

use crate::extended_euclid;

fn mod_pos(a: BigInt, m: &BigInt) -> BigInt {
    if m.is_zero() {
        return a;
    }
    a.mod_floor(m)
}

pub fn hermite_modulo_d_general(mut a: Vec<Vec<BigInt>>, d: BigInt) -> Option<Vec<Vec<BigInt>>> {
    let m = a.len();
    if m == 0 || d.is_zero() {
        return None;
    }
    let n = a[0].len();
    if a.iter().any(|row| row.len() != n) || n < m {
        return None;
    }

    let mut r_mod = d;
    let mut i: isize = m as isize - 1;
    let mut k: isize = n as isize - 1;
    let mut j: isize = k;

    while i >= 0 && k >= 0 {
        if j <= 0 {
            let aik = a[i as usize][k as usize].clone();
            let (u, _v, g) = extended_euclid(aik.clone(), r_mod.clone());
            let mut wi = vec![BigInt::zero(); m];
            for row in 0..m {
                wi[row] = mod_pos(u.clone() * a[row][k as usize].clone(), &r_mod);
            }
            if wi[i as usize].is_zero() {
                wi[i as usize] = g.clone();
            }
            for row in 0..m {
                a[row][i as usize] = wi[row].clone();
            }

            let pivot = wi[i as usize].clone();
            if pivot.is_zero() {
                return None;
            }
            for col in (i as usize + 1)..m {
                let q = a[col][i as usize].div_floor(&pivot);
                if q.is_zero() {
                    continue;
                }
                for row in 0..m {
                    a[row][col] = a[row][col].clone() - &q * &wi[row];
                }
            }

            if i == 0 {
                let mut w = vec![vec![BigInt::zero(); m]; m];
                for row in 0..m {
                    for col in 0..m {
                        w[row][col] = a[row][col].clone();
                    }
                }
                return Some(w);
            }

            if g.is_zero() {
                return None;
            }
            r_mod = r_mod / g;
            i -= 1;
            k -= 1;
            j = k;
            continue;
        } else {
            j -= 1;
            if a[i as usize][j as usize].is_zero() {
                continue;
            }
            let aik = a[i as usize][k as usize].clone();
            let aij = a[i as usize][j as usize].clone();
            let (u, v, g) = extended_euclid(aik.clone(), aij.clone());

            let mut b_col = vec![BigInt::zero(); m];
            for row in 0..m {
                b_col[row] = &u * a[row][k as usize].clone() + &v * a[row][j as usize].clone();
            }

            let fac_k = aik.div_floor(&g);
            let fac_j = aij.div_floor(&g);
            for row in 0..m {
                let val = fac_k.clone() * a[row][j as usize].clone()
                    - fac_j.clone() * a[row][k as usize].clone();
                a[row][j as usize] = mod_pos(val, &r_mod);
            }
            for row in 0..m {
                a[row][k as usize] = mod_pos(b_col[row].clone(), &r_mod);
            }
            continue;
        }
    }

    None
}

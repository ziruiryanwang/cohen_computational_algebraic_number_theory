use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, Zero};

use crate::{extended_euclid, gauss_bareiss_det};

fn mod_pos(a: BigInt, m: &BigInt) -> BigInt {
    if m.is_zero() {
        return a;
    }
    a.mod_floor(m)
}

pub fn smith_normal_form(mut a: Vec<Vec<BigInt>>) -> Option<Vec<BigInt>> {
    let n = a.len();
    if n == 0 || a.iter().any(|row| row.len() != n) {
        return None;
    }

    if n == 1 {
        let det = a[0][0].abs();
        return Some(vec![det]);
    }

    let det = gauss_bareiss_det(a.clone())?.abs();
    if det.is_zero() {
        return None;
    }
    let mut r = det;
    let mut i: isize = n as isize - 1;

    while i >= 0 {
        let mut j = i;
        let mut c_count = 0;
        // Row reduction stage
        loop {
            if j <= 0 {
                break;
            }
            j -= 1;
            if a[i as usize][j as usize].is_zero() {
                continue;
            }
            let aii = a[i as usize][i as usize].clone();
            let aij = a[i as usize][j as usize].clone();
            let (u, v, d) = extended_euclid(aii.clone(), aij.clone());
            if d.is_zero() {
                return None;
            }
            let factor_i = aii.div_floor(&d);
            let factor_j = aij.div_floor(&d);

            let mut b_col = vec![BigInt::zero(); n];
            let mut new_col_j = vec![BigInt::zero(); n];
            for row in 0..n {
                b_col[row] = &u * a[row][i as usize].clone() + &v * a[row][j as usize].clone();
                new_col_j[row] = factor_i.clone() * a[row][j as usize].clone()
                    - factor_j.clone() * a[row][i as usize].clone();
                new_col_j[row] = mod_pos(new_col_j[row].clone(), &r);
            }
            for row in 0..n {
                a[row][i as usize] = mod_pos(b_col[row].clone(), &r);
                a[row][j as usize] = new_col_j[row].clone();
            }
            j = i; // restart row sweep
        }

        // Column reduction stage
        j = i;
        loop {
            if j <= 0 {
                break;
            }
            j -= 1;
            if a[j as usize][i as usize].is_zero() {
                continue;
            }
            let aii = a[i as usize][i as usize].clone();
            let aji = a[j as usize][i as usize].clone();
            let (u, v, d) = extended_euclid(aii.clone(), aji.clone());
            if d.is_zero() {
                return None;
            }
            let factor_i = aii.div_floor(&d);
            let factor_j = aji.div_floor(&d);

            let mut b_row = vec![BigInt::zero(); n];
            let mut new_row_j = vec![BigInt::zero(); n];
            for col in 0..n {
                b_row[col] = &u * a[i as usize][col].clone() + &v * a[j as usize][col].clone();
                new_row_j[col] = factor_i.clone() * a[j as usize][col].clone()
                    - factor_j.clone() * a[i as usize][col].clone();
                new_row_j[col] = mod_pos(new_row_j[col].clone(), &r);
            }
            for col in 0..n {
                a[i as usize][col] = mod_pos(b_row[col].clone(), &r);
                a[j as usize][col] = new_row_j[col].clone();
            }

            c_count += 1;
            if c_count > 0 {
                j = i;
            }
        }

        if c_count > 0 {
            i += 0; // remain same i, restart
            continue;
        }

        let b = a[i as usize][i as usize].clone();
        for row in 0..(i as usize) {
            for col in 0..(i as usize) {
                if !b.is_zero() && (&a[row][col] % &b) != BigInt::zero() {
                    for col_i in 0..n {
                        a[i as usize][col_i] =
                            mod_pos(a[i as usize][col_i].clone() + a[row][col_i].clone(), &r);
                    }
                    for row_i in 0..n {
                        a[row_i][i as usize] =
                            mod_pos(a[row_i][i as usize].clone() + a[row_i][col].clone(), &r);
                    }
                    j = i;
                    continue;
                }
            }
        }

        let di = {
            let g = num_integer::Integer::gcd(&b, &r);
            if g.is_zero() { b.abs() } else { g.abs() }
        };
        let mut diag = Vec::new();
        diag.push(di.clone());
        r = r / di;
        if i == 1 {
            let d1 = num_integer::Integer::gcd(&a[0][0], &r).abs();
            diag.push(d1);
            diag.reverse();
            return Some(diag);
        } else {
            i -= 1;
        }
    }

    None
}

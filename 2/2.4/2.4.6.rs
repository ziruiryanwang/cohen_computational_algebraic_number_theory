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

pub fn hermite_modulo_d(l: Vec<Vec<BigInt>>, d: BigInt) -> Option<Vec<Vec<BigInt>>> {
    let m = l.len();
    if m == 0 || l.iter().any(|row| row.len() != m) || d.is_zero() {
        return None;
    }
    let mut w = l.clone();
    let mut b = d;
    let mut i: isize = m as isize - 1;

    loop {
        let lii = w[i as usize][i as usize].clone();
        let (u, _v, g) = extended_euclid(lii.clone(), b.clone());
        for r in 0..m {
            let val = u.clone() * w[r][i as usize].clone();
            w[r][i as usize] = mod_pos(val, &b);
        }
        if g == b {
            w[i as usize][i as usize] = g.clone();
        }

        if i > 0 {
            if g.is_zero() {
                return None;
            }
            b = b / g;
            i -= 1;
            continue;
        } else {
            break;
        }
    }

    for ii in (0..m.saturating_sub(1)).rev() {
        let diag = w[ii][ii].clone();
        if diag.is_zero() {
            return None;
        }
        for j in (ii + 1)..m {
            let q = w[ii][j].div_floor(&diag);
            if q.is_zero() {
                continue;
            }
            for r in 0..m {
                w[r][j] = w[r][j].clone() - &q * &w[r][ii];
            }
        }
    }

    Some(w)
}

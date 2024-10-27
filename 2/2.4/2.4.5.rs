use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, Zero};

fn extended_gcd(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    let mut old_r = a;
    let mut r = b;
    let mut old_s = BigInt::from(1);
    let mut s = BigInt::from(0);
    let mut old_t = BigInt::from(0);
    let mut t = BigInt::from(1);

    while !r.is_zero() {
        let q = &old_r / &r;
        let tmp_r = old_r - &q * &r;
        old_r = r;
        r = tmp_r;

        let tmp_s = old_s - &q * &s;
        old_s = s;
        s = tmp_s;

        let tmp_t = old_t - &q * &t;
        old_t = t;
        t = tmp_t;
    }

    if old_r.is_negative() {
        (-old_s, -old_t, -old_r)
    } else {
        (old_s, old_t, old_r)
    }
}

pub fn hermite_normal_form_euclid(mut a: Vec<Vec<BigInt>>) -> Option<Vec<Vec<BigInt>>> {
    let m = a.len();
    if m == 0 {
        return Some(a);
    }
    let n = a[0].len();
    if a.iter().any(|row| row.len() != n) {
        return None;
    }

    let l = if m <= n { 0isize } else { (m - n) as isize };
    let mut i = m as isize - 1;
    let mut k = n as isize - 1;
    let mut j = k;

    loop {
        if k < 0 {
            return Some(a);
        }

        loop {
            if j == 0 {
                break;
            }
            j -= 1;
            if a[i as usize][j as usize].is_zero() {
                continue;
            }

            let a_ik = a[i as usize][k as usize].clone();
            let a_ij = a[i as usize][j as usize].clone();
            let (u, v, d) = extended_gcd(a_ik.clone(), a_ij.clone());

            let mut b_col = vec![BigInt::zero(); m];
            for row in 0..m {
                b_col[row] = &u * &a[row][k as usize] + &v * &a[row][j as usize];
            }

            let factor_j = a_ik.div_floor(&d);
            let factor_k = a_ij.div_floor(&d);
            for row in 0..m {
                a[row][j as usize] =
                    &factor_j * &a[row][j as usize] - &factor_k * &a[row][k as usize];
            }

            for row in 0..m {
                a[row][k as usize] = b_col[row].clone();
            }

            continue;
        }

        let mut b = a[i as usize][k as usize].clone();
        if b.is_negative() {
            for row in 0..m {
                a[row][k as usize] = -a[row][k as usize].clone();
            }
            b = -b;
        }

        if b.is_zero() {
            k += 1;
        } else {
            for col in ((k + 1) as usize)..n {
                let q = a[i as usize][col].div_floor(&b);
                if q.is_zero() {
                    continue;
                }
                for row in 0..m {
                    a[row][col] = a[row][col].clone() - &q * &a[row][k as usize];
                }
            }
        }

        if i == l {
            let start = k as usize;
            let mut w = vec![vec![BigInt::zero(); n - start]; m];
            for row in 0..m {
                for col in start..n {
                    w[row][col - start] = a[row][col].clone();
                }
            }
            return Some(w);
        } else {
            i -= 1;
            k -= 1;
            if k < 0 {
                return Some(a);
            }
            j = k;
        }
    }
}

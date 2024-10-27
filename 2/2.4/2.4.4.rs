use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, Zero};

pub fn hermite_normal_form(mut a: Vec<Vec<BigInt>>) -> Option<Vec<Vec<BigInt>>> {
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

    loop {
        if k < 0 {
            return Some(a);
        }

        let mut row_finished = true;
        for j in 0..k {
            if !a[i as usize][j as usize].is_zero() {
                row_finished = false;
                break;
            }
        }
        if row_finished {
            if a[i as usize][k as usize].is_negative() {
                for row in 0..m {
                    a[row][k as usize] = -a[row][k as usize].clone();
                }
            }
            let b = a[i as usize][k as usize].clone();
            if b.is_zero() {
                k += 1;
            } else {
                for j in ((k + 1) as usize)..n {
                    let q = a[i as usize][j].div_floor(&b);
                    if q.is_zero() {
                        continue;
                    }
                    for row in 0..m {
                        a[row][j] = a[row][j].clone() - &q * &a[row][k as usize];
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
            }
            continue;
        }

        let mut j0 = None;
        let mut min_abs: Option<BigInt> = None;
        for j in 0..=k as usize {
            let val = a[i as usize][j].clone();
            if val.is_zero() {
                continue;
            }
            let abs_val = val.abs();
            match &min_abs {
                None => {
                    min_abs = Some(abs_val);
                    j0 = Some(j);
                }
                Some(curr) => {
                    if abs_val < *curr {
                        min_abs = Some(abs_val);
                        j0 = Some(j);
                    }
                }
            }
        }
        let j0 = j0?;

        if j0 < k as usize {
            for row in 0..m {
                a[row].swap(j0, k as usize);
            }
        }
        if a[i as usize][k as usize].is_negative() {
            for row in 0..m {
                a[row][k as usize] = -a[row][k as usize].clone();
            }
        }
        let b = a[i as usize][k as usize].clone();
        for j in 0..k as usize {
            let q = a[i as usize][j].div_floor(&b);
            if q.is_zero() {
                continue;
            }
            for row in 0..m {
                a[row][j] = a[row][j].clone() - &q * &a[row][k as usize];
            }
        }
    }
}

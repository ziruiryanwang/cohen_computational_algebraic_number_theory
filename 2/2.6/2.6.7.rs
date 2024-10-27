use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

fn identity_bigint(n: usize) -> Vec<Vec<BigInt>> {
    let mut m = vec![vec![BigInt::zero(); n]; n];
    for i in 0..n {
        m[i][i] = BigInt::one();
    }
    m
}

fn dot(a: &[BigInt], b: &[BigInt]) -> BigInt {
    let mut s = BigInt::zero();
    for (x, y) in a.iter().zip(b.iter()) {
        s += x * y;
    }
    s
}

fn nearest_integer(r: &BigRational) -> BigInt {
    let num = r.numer().clone();
    let den = r.denom().clone();
    let (q, rem) = num.div_rem(&den);
    let double = &rem * 2u32;
    if double.abs() < den {
        q
    } else if double.abs() > den {
        if num.is_negative() {
            q - BigInt::one()
        } else {
            q + BigInt::one()
        }
    } else if q.is_even() {
        q
    } else if num.is_negative() {
        q - BigInt::one()
    } else {
        q + BigInt::one()
    }
}

fn recompute(
    b: &Vec<Vec<BigInt>>,
    upto: usize,
    lambda: &mut Vec<Vec<BigInt>>,
    d: &mut Vec<BigInt>,
) -> Option<()> {
    let n = upto + 1;
    if b.is_empty() {
        return Some(());
    }
    d[0] = dot(&b[0], &b[0]);
    if d[0].is_zero() {
        return None;
    }
    for k in 1..n {
        for j in 0..=k {
            let mut u = dot(&b[k], &b[j]);
            for i in 0..j {
                let di = d[i].clone();
                let dim1 = if i == 0 { BigInt::one() } else { d[i - 1].clone() };
                let lki = if k > i { lambda[k][i].clone() } else { BigInt::zero() };
                let lji = if j > i { lambda[j][i].clone() } else { BigInt::zero() };
                let numerator = di * u - lki * lji;
                u = numerator.div_floor(&dim1);
            }
            if j < k {
                lambda[k][j] = u.clone();
            } else {
                d[k] = u;
            }
        }
        if d[k].is_zero() {
            return None;
        }
    }
    Some(())
}

fn redi(
    b: &mut Vec<Vec<BigInt>>,
    h: &mut Vec<Vec<BigInt>>,
    lambda: &mut Vec<Vec<BigInt>>,
    d: &mut Vec<BigInt>,
    k: usize,
    l: usize,
) -> Option<()> {
    let denom = d[l].clone();
    if denom.is_zero() {
        return None;
    }
    let mu = BigRational::from_integer(lambda[k][l].clone()) / BigRational::from_integer(denom);
    let q = nearest_integer(&mu);
    if q.is_zero() {
        return Some(());
    }
    let bl = b[l].clone();
    let bk = &mut b[k];
    for idx in 0..bk.len() {
        bk[idx] -= &q * &bl[idx];
    }
    let hl: Vec<BigInt> = h.iter().map(|row| row[l].clone()).collect();
    for (idx, row) in h.iter_mut().enumerate() {
        row[k] -= &q * &hl[idx];
    }
    recompute(b, k, lambda, d)
}

fn swap(
    b: &mut Vec<Vec<BigInt>>,
    h: &mut Vec<Vec<BigInt>>,
    lambda: &mut Vec<Vec<BigInt>>,
    d: &mut Vec<BigInt>,
    k: usize,
) -> Option<()> {
    b.swap(k, k - 1);
    for row in 0..h.len() {
        h[row].swap(k, k - 1);
    }
    recompute(b, k, lambda, d)
}

pub fn lll_integral(
    basis: Vec<Vec<BigInt>>,
) -> Option<(Vec<Vec<BigInt>>, Vec<Vec<BigInt>>)> {
    let n = basis.len();
    if n == 0 {
        return Some((basis, Vec::new()));
    }
    let m = basis[0].len();
    if basis.iter().any(|v| v.len() != m) {
        return None;
    }

    let mut b = basis;
    let mut h = identity_bigint(n);
    let mut lambda = vec![vec![BigInt::zero(); n]; n];
    let mut d = vec![BigInt::zero(); n];

    let mut k: usize = 1;
    if recompute(&b, 1.min(n - 1), &mut lambda, &mut d).is_none() {
        return None;
    }
    let mut k_max = 1usize;

    while k < n {
        if k > k_max {
            if recompute(&b, k, &mut lambda, &mut d).is_none() {
                return None;
            }
            k_max = k;
        }

        if redi(&mut b, &mut h, &mut lambda, &mut d, k, k - 1).is_none() {
            return None;
        }

        let condition = if k >= 2 {
            let left = BigRational::from_integer(d[k].clone() * d[k - 2].clone());
            let three = BigInt::from(3);
            let four = BigInt::from(4);
            let right_num =
                three * &d[k - 1] * &d[k - 1] - four.clone() * &lambda[k][k - 1] * &lambda[k][k - 1];
            let right = BigRational::from_integer(right_num) / BigRational::from_integer(four);
            left < right
        } else {
            false
        };

        if condition {
            if swap(&mut b, &mut h, &mut lambda, &mut d, k).is_none() {
                return None;
            }
            k_max = k;
            if k > 1 {
                k -= 1;
            }
            continue;
        }

        for l in (0..k - 1).rev() {
            if redi(&mut b, &mut h, &mut lambda, &mut d, k, l).is_none() {
                return None;
            }
        }

        k += 1;
    }

    Some((b, h))
}

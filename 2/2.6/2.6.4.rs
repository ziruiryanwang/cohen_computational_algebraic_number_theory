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

fn dot_int(a: &[BigInt], b: &[BigInt]) -> BigInt {
    let mut s = BigInt::zero();
    for (x, y) in a.iter().zip(b.iter()) {
        s += x * y;
    }
    s
}

fn dot_rat(a: &[BigRational], b: &[BigRational]) -> BigRational {
    let mut s = BigRational::zero();
    for (x, y) in a.iter().zip(b.iter()) {
        s += x * y;
    }
    s
}

fn to_rat(v: &[BigInt]) -> Vec<BigRational> {
    v.iter()
        .map(|x| BigRational::from_integer(x.clone()))
        .collect()
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
    b_star: &mut Vec<Vec<BigRational>>,
    b_norm: &mut Vec<BigRational>,
    mu: &mut Vec<Vec<BigRational>>,
) -> Option<()> {
    let m = b.first().map(|v| v.len()).unwrap_or(0);
    for i in 0..=upto {
        let mut bi_star = to_rat(&b[i]);
        for j in 0..i {
            let numerator = BigRational::from_integer(dot_int(&b[i], &b[j]));
            let mu_ij = numerator / b_norm[j].clone();
            mu[i][j] = mu_ij.clone();
            for t in 0..m {
                bi_star[t] -= &mu_ij * &b_star[j][t];
            }
        }
        let norm = dot_rat(&bi_star, &bi_star);
        if norm.is_zero() {
            return None;
        }
        b_star[i] = bi_star;
        b_norm[i] = norm;
    }
    Some(())
}

fn red(
    b: &mut Vec<Vec<BigInt>>,
    h: &mut Vec<Vec<BigInt>>,
    mu: &mut Vec<Vec<BigRational>>,
    k: usize,
    l: usize,
) {
    let q = nearest_integer(&mu[k][l]);
    if q.is_zero() {
        return;
    }
    let bl = b[l].clone();
    let bk = &mut b[k];
    for i in 0..bk.len() {
        bk[i] -= &q * &bl[i];
    }
    let hl: Vec<BigInt> = h.iter().map(|row| row[l].clone()).collect();
    for (idx, row) in h.iter_mut().enumerate() {
        row[k] -= &q * &hl[idx];
    }
}

fn insert(
    b: &mut Vec<Vec<BigInt>>,
    h: &mut Vec<Vec<BigInt>>,
    k: usize,
    i: usize,
) {
    if k <= i {
        return;
    }
    let bk = b[k].clone();
    let mut hk = vec![BigInt::zero(); h.len()];
    for row in 0..h.len() {
        hk[row] = h[row][k].clone();
    }
    for idx in (i + 1..=k).rev() {
        b[idx] = b[idx - 1].clone();
        for row in 0..h.len() {
            h[row][idx] = h[row][idx - 1].clone();
        }
    }
    b[i] = bk;
    for row in 0..h.len() {
        h[row][i] = hk[row].clone();
    }
}

pub fn lll_deep(
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

    let delta = BigRational::new(BigInt::from(3), BigInt::from(4));
    let mut b = basis;
    let mut h = identity_bigint(n);
    let mut b_star = vec![vec![BigRational::zero(); m]; n];
    let mut b_norm = vec![BigRational::zero(); n];
    let mut mu = vec![vec![BigRational::zero(); n]; n];

    let mut k: usize = 1;
    while k < n {
        if recompute(&b, k, &mut b_star, &mut b_norm, &mut mu).is_none() {
            return None;
        }
        if b_norm[k].is_zero() {
            return None;
        }

        for l in (0..k).rev() {
            red(&mut b, &mut h, &mut mu, k, l);
            if recompute(&b, k, &mut b_star, &mut b_norm, &mut mu).is_none() {
                return None;
            }
        }

        let mut b_proj = BigRational::from_integer(dot_int(&b[k], &b[k]));
        let mut i: usize = 0;
        loop {
            if i == k {
                k += 1;
                break;
            }
            let lhs = (&delta * b_norm[i].clone()) * BigRational::new(BigInt::one(), BigInt::from(1));
            if lhs <= b_proj {
                let mu_sq = &mu[k][i] * &mu[k][i];
                b_proj -= mu_sq * b_norm[i].clone();
                i += 1;
                continue;
            } else {
                insert(&mut b, &mut h, k, i);
                if i >= 1 {
                    k = i - 1;
                    b_proj = BigRational::from_integer(dot_int(&b[k], &b[k]));
                    i = 0;
                    if recompute(&b, k, &mut b_star, &mut b_norm, &mut mu).is_none() {
                        return None;
                    }
                    continue;
                } else {
                    k = 0;
                    break;
                }
            }
        }
    }

    Some((b, h))
}

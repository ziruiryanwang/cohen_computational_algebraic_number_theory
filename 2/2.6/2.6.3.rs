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
    let (num, den) = (r.numer().clone(), r.denom().clone());
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
    } else {
        if q.is_even() {
            q
        } else if num.is_negative() {
            q - BigInt::one()
        } else {
            q + BigInt::one()
        }
    }
}

fn recompute(
    b: &Vec<Vec<BigInt>>,
    upto: usize,
    b_star: &mut Vec<Vec<BigRational>>,
    b_norm: &mut Vec<BigRational>,
    mu: &mut Vec<Vec<BigRational>>,
) -> Option<()> {
    for i in 0..=upto {
        let mut bi_star = to_rat(&b[i]);
        for j in 0..i {
            let numerator = BigRational::from_integer(dot_int(&b[i], &b[j]));
            let mu_ij = numerator / b_norm[j].clone();
            mu[i][j] = mu_ij.clone();
            for t in 0..bi_star.len() {
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

pub fn lll_reduction(
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

        if !mu[k][k - 1].is_zero() {
            let q = nearest_integer(&mu[k][k - 1]);
            if !q.is_zero() {
                for i in 0..m {
                    let t = &q * &b[k - 1][i];
                    b[k][i] -= t;
                }
                for i in 0..n {
                    let t = &q * &h[i][k - 1];
                    h[i][k] -= t;
                }
                if recompute(&b, k, &mut b_star, &mut b_norm, &mut mu).is_none() {
                    return None;
                }
            }
        }

        let lhs = b_norm[k].clone();
        let mu_sq = &mu[k][k - 1] * &mu[k][k - 1];
        let rhs = (delta.clone() - mu_sq) * b_norm[k - 1].clone();

        if lhs < rhs {
            b.swap(k, k - 1);
            h.swap(k, k - 1);
            if k > 1 {
                k -= 1;
            }
            continue;
        }

        for l in (0..(k - 1)).rev() {
            if mu[k][l].is_zero() {
                continue;
            }
            let q = nearest_integer(&mu[k][l]);
            if q.is_zero() {
                continue;
            }
            for i in 0..m {
                let t = &q * &b[l][i];
                b[k][i] -= t;
            }
            for i in 0..n {
                let t = &q * &h[i][l];
                h[i][k] -= t;
            }
            if recompute(&b, k, &mut b_star, &mut b_norm, &mut mu).is_none() {
                return None;
            }
        }

        k += 1;
    }

    Some((b, h))
}

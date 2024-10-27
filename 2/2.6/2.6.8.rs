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

fn to_rat(v: &[BigInt]) -> Vec<BigRational> {
    v.iter()
        .map(|x| BigRational::from_integer(x.clone()))
        .collect()
}

fn dot_int_rat(a: &[BigInt], b: &[BigRational]) -> BigRational {
    let mut s = BigRational::zero();
    for (x, y) in a.iter().zip(b.iter()) {
        s += BigRational::from_integer(x.clone()) * y;
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
    norms: &mut Vec<BigRational>,
    mu: &mut Vec<Vec<BigRational>>,
) -> Option<()> {
    let dim = b.first().map(|v| v.len()).unwrap_or(0);
    for i in 0..=upto {
        let mut bi_star = to_rat(&b[i]);
        for j in 0..i {
            if norms[j].is_zero() {
                mu[i][j] = BigRational::zero();
                continue;
            }
            let numerator = dot_int_rat(&b[i], &b_star[j]);
            let mu_ij = numerator / norms[j].clone();
            mu[i][j] = mu_ij.clone();
            for t in 0..dim {
                bi_star[t] -= &mu_ij * &b_star[j][t];
            }
        }
        let norm = dot_rat(&bi_star, &bi_star);
        b_star[i] = bi_star;
        norms[i] = norm;
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

fn swapg(b: &mut Vec<Vec<BigInt>>, h: &mut Vec<Vec<BigInt>>, k: usize) {
    b.swap(k, k - 1);
    for row in h.iter_mut() {
        row.swap(k, k - 1);
    }
}

pub fn lll_dependent(
    basis: Vec<Vec<BigInt>>,
) -> Option<(Vec<Vec<BigInt>>, Vec<Vec<BigInt>>, usize)> {
    let n = basis.len();
    if n == 0 {
        return Some((basis, Vec::new(), 0));
    }
    let dim = basis[0].len();
    if basis.iter().any(|v| v.len() != dim) {
        return None;
    }

    let delta = BigRational::new(BigInt::from(3), BigInt::from(4));
    let mut b = basis;
    let mut h = identity_bigint(n);
    let mut b_star = vec![vec![BigRational::zero(); dim]; n];
    let mut norms = vec![BigRational::zero(); n];
    let mut mu = vec![vec![BigRational::zero(); n]; n];

    if recompute(&b, 0, &mut b_star, &mut norms, &mut mu).is_none() {
        return None;
    }
    let mut k: usize = 1;
    let mut k_max: usize = 0;

    while k < n {
        if k > k_max {
            if recompute(&b, k, &mut b_star, &mut norms, &mut mu).is_none() {
                return None;
            }
            k_max = k;
        }

        if k == 0 {
            k = 1;
            continue;
        }
        red(&mut b, &mut h, &mut mu, k, k - 1);
        if recompute(&b, k, &mut b_star, &mut norms, &mut mu).is_none() {
            return None;
        }

        let lovasz = if norms[k - 1].is_zero() {
            false
        } else {
            let lhs = norms[k].clone();
            let mu_sq = &mu[k][k - 1] * &mu[k][k - 1];
            let rhs = (delta.clone() - mu_sq) * norms[k - 1].clone();
            lhs < rhs
        };

        if lovasz {
            swapg(&mut b, &mut h, k);
            if k > 1 {
                k -= 1;
            }
            k_max = k_max.max(k);
            continue;
        }

        for l in (0..k - 1).rev() {
            red(&mut b, &mut h, &mut mu, k, l);
            if recompute(&b, k, &mut b_star, &mut norms, &mut mu).is_none() {
                return None;
            }
        }

        k += 1;
    }

    let r = norms.iter().filter(|x| x.is_zero()).count();
    let p = n - r;
    Some((b, h, p))
}

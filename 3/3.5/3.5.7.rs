use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

use crate::berlekamp_factorization;

fn trim(mut p: Vec<BigInt>) -> Vec<BigInt> {
    while let Some(last) = p.last() {
        if last.is_zero() {
            p.pop();
        } else {
            break;
        }
    }
    p
}

fn degree(p: &Vec<BigInt>) -> isize {
    if p.is_empty() {
        -1
    } else {
        (p.len() - 1) as isize
    }
}

fn content(p: &Vec<BigInt>) -> BigInt {
    let mut g = BigInt::zero();
    for c in p.iter() {
        g = g.gcd(c);
    }
    g
}

fn scalar_div(vec: Vec<BigInt>, d: &BigInt) -> Vec<BigInt> {
    vec.into_iter().map(|c| c / d).collect()
}

fn primitive_part(p: Vec<BigInt>) -> Vec<BigInt> {
    let g = content(&p);
    if g.is_zero() || g.is_one() {
        return trim(p);
    }
    trim(scalar_div(p, &g))
}

fn add_mod_p(a: Vec<BigInt>, b: Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    let n = a.len().max(b.len());
    let mut res = vec![BigInt::zero(); n];
    for i in 0..n {
        let ai = a.get(i).cloned().unwrap_or_else(BigInt::zero);
        let bi = b.get(i).cloned().unwrap_or_else(BigInt::zero);
        res[i] = (ai + bi).mod_floor(p);
    }
    trim(res)
}

fn sub_mod_p(a: Vec<BigInt>, b: Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    let n = a.len().max(b.len());
    let mut res = vec![BigInt::zero(); n];
    for i in 0..n {
        let ai = a.get(i).cloned().unwrap_or_else(BigInt::zero);
        let bi = b.get(i).cloned().unwrap_or_else(BigInt::zero);
        res[i] = (ai - bi).mod_floor(p);
    }
    trim(res)
}

fn mul_mod_p(a: &Vec<BigInt>, b: &Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut res = vec![BigInt::zero(); a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            res[i + j] = (res[i + j].clone() + &a[i] * &b[j]).mod_floor(p);
        }
    }
    trim(res)
}

fn div_rem_mod_p(mut a: Vec<BigInt>, b: Vec<BigInt>, p: &BigInt) -> (Vec<BigInt>, Vec<BigInt>) {
    let b = trim(b);
    let mut q = Vec::new();
    let deg_b = degree(&b);
    if deg_b < 0 {
        return (Vec::new(), a);
    }
    let lc_b = b.last().cloned().unwrap();
    let inv_lc = lc_b.modpow(&(p - 2), p);
    while degree(&a) >= deg_b && !a.is_empty() {
        let shift = (degree(&a) - deg_b) as usize;
        let coeff = (a.last().cloned().unwrap() * &inv_lc).mod_floor(p);
        if q.len() <= shift {
            q.resize(shift + 1, BigInt::zero());
        }
        q[shift] = (q[shift].clone() + coeff.clone()).mod_floor(p);
        for i in 0..=deg_b as usize {
            let idx = i + shift;
            let sub = (&b[i] * &coeff).mod_floor(p);
            a[idx] = (a[idx].clone() - sub).mod_floor(p);
        }
        a = trim(a);
    }
    (trim(q), trim(a))
}

fn gcd_mod_p(mut a: Vec<BigInt>, mut b: Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    a = trim(a);
    b = trim(b);
    while !b.is_empty() {
        let r = div_rem_mod_p(a.clone(), b.clone(), p).1;
        a = b;
        b = r;
    }
    monic_mod_p(a, p)
}

fn monic_mod_p(poly: Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    let poly = trim(poly);
    if poly.is_empty() {
        return poly;
    }
    let lc = poly.last().cloned().unwrap();
    if lc.is_one() {
        return poly;
    }
    let inv = lc.modpow(&(p - 2), p);
    poly.into_iter().map(|c| (c * &inv).mod_floor(p)).collect()
}

fn mod_poly_coeffs(a: Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    a.into_iter().map(|c| c.mod_floor(p)).collect()
}

fn squarefree_mod_p(a: Vec<BigInt>, p: &BigInt) -> Option<Vec<BigInt>> {
    let a = mod_poly_coeffs(a, p);
    let der = deriv_mod_p(&a, p);
    let g = gcd_mod_p(a.clone(), der, p);
    if degree(&g) > 0 {
        return None;
    }
    Some(a)
}

fn deriv_mod_p(poly: &Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    if poly.len() <= 1 {
        return Vec::new();
    }
    let mut res = Vec::with_capacity(poly.len() - 1);
    for i in 1..poly.len() {
        let coeff = (&poly[i] * BigInt::from(i)).mod_floor(p);
        res.push(coeff);
    }
    trim(res)
}

pub fn factor_over_z(a: Vec<BigInt>) -> Vec<Vec<BigInt>> {
    let mut factors = Vec::new();
    let cont = content(&a);
    if !cont.is_zero() && !cont.is_one() {
        factors.push(vec![cont.clone()]);
    }
    let mut u = primitive_part(a);
    if degree(&u) <= 0 {
        if !u.is_empty() {
            factors.push(u);
        }
        return factors;
    }

    let primes = [2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29];
    for &pp in primes.iter() {
        let p = BigInt::from(pp);
        if let Some(up) = squarefree_mod_p(u.clone(), &p) {
            let fac_mod_p = berlekamp_factorization(up, &p);
            if fac_mod_p.len() > 1 {
                for f in fac_mod_p {
                    let lifted = lift_from_mod_p(f, &p);
                    factors.push(lifted.clone());
                }
                return factors;
            }
        }
    }

    factors.push(u);
    factors
}

fn lift_from_mod_p(f: Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    f.into_iter()
        .map(|c| {
            let mut v = c.mod_floor(p);
            if v > p.clone() / BigInt::from(2) {
                v = v - p;
            }
            v
        })
        .collect()
}

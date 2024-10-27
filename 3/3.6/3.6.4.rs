use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

use crate::{factor_over_z, subresultant_gcd};

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

fn deriv(p: &Vec<BigInt>) -> Vec<BigInt> {
    if p.len() <= 1 {
        return Vec::new();
    }
    let mut d = Vec::with_capacity(p.len() - 1);
    for i in 1..p.len() {
        d.push(&p[i] * BigInt::from(i as i64));
    }
    trim(d)
}

fn is_squarefree(p: &Vec<BigInt>) -> bool {
    let p = trim(p.clone());
    if degree(&p) <= 0 {
        return true;
    }
    let d = deriv(&p);
    if d.is_empty() {
        return true;
    }
    let g = subresultant_gcd(p, d);
    degree(&g) <= 0
}

fn binom(n: usize, k: usize) -> BigInt {
    if k == 0 || k == n {
        return BigInt::one();
    }
    let mut res = BigInt::one();
    let mut kk = k.min(n - k);
    for i in 1..=kk {
        res *= BigInt::from((n - kk + i) as i64);
        res /= BigInt::from(i as i64);
    }
    res
}

fn shift_poly(p: &Vec<BigInt>, k: &BigInt) -> Vec<BigInt> {
    if k.is_zero() {
        return trim(p.clone());
    }
    let mut res: Vec<BigInt> = Vec::new();
    for (i, coeff) in p.iter().enumerate() {
        if coeff.is_zero() {
            continue;
        }
        for j in 0..=i {
            let comb = binom(i, j);
            let k_pow = k.pow((i - j) as u32);
            let term = coeff * &comb * k_pow;
            if res.len() <= j {
                res.resize(j + 1, BigInt::zero());
            }
            res[j] += term;
        }
    }
    trim(res)
}

fn div_exact(mut a: Vec<BigInt>, b: Vec<BigInt>) -> Option<Vec<BigInt>> {
    let b = trim(b);
    if b.is_empty() {
        return None;
    }
    let deg_b = degree(&b);
    let lc_b = b.last().cloned().unwrap();
    let mut q: Vec<BigInt> = Vec::new();
    a = trim(a);
    while degree(&a) >= deg_b && !a.is_empty() {
        let shift = (degree(&a) - deg_b) as usize;
        let lc_a = a.last().cloned().unwrap();
        if !lc_a.is_multiple_of(&lc_b) {
            return None;
        }
        let coeff = lc_a / &lc_b;
        if q.len() <= shift {
            q.resize(shift + 1, BigInt::zero());
        }
        q[shift] += coeff.clone();
        for i in 0..=deg_b as usize {
            let idx = i + shift;
            a[idx] -= &b[i] * &coeff;
        }
        a = trim(a);
    }
    if a.is_empty() {
        Some(trim(q))
    } else {
        None
    }
}

fn squarefree_part(a: Vec<BigInt>) -> Vec<BigInt> {
    let a = trim(a);
    if degree(&a) <= 0 {
        return a;
    }
    let d = deriv(&a);
    if d.is_empty() {
        return vec![BigInt::one()];
    }
    let g = subresultant_gcd(a.clone(), d);
    if degree(&g) <= 0 {
        return a;
    }
    div_exact(a, g).unwrap_or_else(|| Vec::new())
}

pub fn factor_over_number_field(_t: Vec<BigInt>, a: Vec<BigInt>) -> Vec<Vec<BigInt>> {
    let a = trim(a);
    if a.is_empty() {
        return Vec::new();
    }

    let u = squarefree_part(a.clone());
    let mut k = BigInt::zero();
    let shifted = loop {
        let n = shift_poly(&u, &k);
        if is_squarefree(&n) {
            break n;
        }
        k += 1;
    };

    let mut factors = Vec::new();
    for f in factor_over_z(shifted) {
        let ai = shift_poly(&f, &(-k.clone()));
        if ai.is_empty() {
            continue;
        }
        let mut multiplicity = 0usize;
        let mut cur = a.clone();
        while let Some(q) = div_exact(cur.clone(), ai.clone()) {
            multiplicity += 1;
            cur = q;
        }
        if multiplicity == 0 {
            multiplicity = 1;
        }
        for _ in 0..multiplicity {
            factors.push(ai.clone());
        }
    }
    factors
}

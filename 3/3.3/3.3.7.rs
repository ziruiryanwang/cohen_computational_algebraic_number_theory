use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

use crate::pseudo_division;

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

pub fn resultant_subresultant(mut a: Vec<BigInt>, mut b: Vec<BigInt>) -> BigInt {
    a = trim(a);
    b = trim(b);
    if a.is_empty() || b.is_empty() {
        return BigInt::zero();
    }

    let mut ca = content(&a);
    let mut cb = content(&b);
    if !ca.is_zero() {
        a = scalar_div(a, &ca);
    }
    if !cb.is_zero() {
        b = scalar_div(b, &cb);
    }

    let mut g = BigInt::one();
    let mut h = BigInt::one();
    let mut s = BigInt::one();
    let mut t = ca.pow(degree(&b).max(0) as u32) * cb.pow(degree(&a).max(0) as u32);

    if degree(&a) < degree(&b) {
        std::mem::swap(&mut a, &mut b);
        let deg_a = degree(&a);
        let deg_b = degree(&b);
        if deg_a % 2 != 0 && deg_b % 2 != 0 {
            s = -s;
        }
    }

    loop {
        let deg_a = degree(&a);
        let deg_b = degree(&b);
        if deg_b < 0 {
            break;
        }
        if deg_b == 0 {
            break;
        }

        if deg_a % 2 != 0 && deg_b % 2 != 0 {
            s = -s;
        }

        let delta = deg_a - deg_b;
        let (_q, mut r) = match pseudo_division(a.clone(), b.clone()) {
            Some(res) => res,
            None => return BigInt::zero(),
        };
        r = trim(r);
        if r.is_empty() {
            break;
        }

        a = b;
        let h_pow = h.pow(delta as u32);
        let gh_delta = g.clone() * h_pow;
        b = scalar_div(r, &gh_delta);
        g = a.last().cloned().unwrap_or_else(BigInt::one);
        let g_pow = g.pow(delta as u32);
        let h_div = if delta > 1 {
            h.pow((delta - 1) as u32)
        } else {
            BigInt::one()
        };
        h = g_pow / h_div;
    }

    let deg_b = degree(&b);
    if deg_b > 0 {
        return BigInt::zero();
    }

    let deg_a = degree(&a);
    if deg_a < 0 {
        return BigInt::zero();
    }

    let h_pow = h.pow(deg_a as u32);
    let l_b = b.last().cloned().unwrap_or_else(BigInt::one);
    let l_a = a.last().cloned().unwrap_or_else(BigInt::one);
    s * t * h_pow * l_b.pow(deg_a as u32) * l_a
}

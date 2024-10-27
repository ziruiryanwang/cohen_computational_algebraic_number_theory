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

fn scalar_mul(vec: Vec<BigInt>, d: &BigInt) -> Vec<BigInt> {
    vec.into_iter().map(|c| c * d).collect()
}

pub fn subresultant_gcd(mut a: Vec<BigInt>, mut b: Vec<BigInt>) -> Vec<BigInt> {
    a = trim(a);
    b = trim(b);
    if degree(&b) > degree(&a) {
        std::mem::swap(&mut a, &mut b);
    }
    if b.is_empty() {
        return a;
    }
    let mut ca = content(&a);
    let mut cb = content(&b);
    let d = ca.gcd(&cb);
    if !ca.is_zero() && !ca.is_one() {
        a = scalar_div(a, &ca);
    }
    if !cb.is_zero() && !cb.is_one() {
        b = scalar_div(b, &cb);
    }
    let mut g = BigInt::one();
    let mut h = BigInt::one();
    let mut c = 0usize;

    loop {
        let delta = degree(&a) - degree(&b);
        let (_q, mut r) = match pseudo_division(a.clone(), b.clone()) {
            Some(res) => res,
            None => return Vec::new(),
        };
        r = trim(r);
        if r.is_empty() {
            break;
        }
        if degree(&r) == 0 {
            b = vec![BigInt::one()];
            break;
        }

        a = b;
        let h_pow_delta = h.pow(delta as u32);
        let denom = g.clone() * h_pow_delta;
        b = scalar_div(r, &denom);
        c += 1;
        if c < 10 {
            g = a.last().cloned().unwrap_or_else(BigInt::one);
            let g_pow = g.pow(delta as u32);
            let h_div = if delta > 1 {
                h.pow((delta - 1) as u32)
            } else {
                BigInt::one()
            };
            h = g_pow / h_div;
        } else {
            ca = content(&a);
            cb = content(&b);
            if !ca.is_zero() && !ca.is_one() {
                a = scalar_div(a, &ca);
            }
            if !cb.is_zero() && !cb.is_one() {
                b = scalar_div(b, &cb);
            }
            g = BigInt::one();
            h = BigInt::one();
            c = 0;
        }
    }

    let cont_b = content(&b);
    if !cont_b.is_zero() && !cont_b.is_one() {
        b = scalar_div(b, &cont_b);
    }
    scalar_mul(b, &d)
}

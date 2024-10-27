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

fn primitive_part(p: Vec<BigInt>) -> Vec<BigInt> {
    let g = content(&p);
    if g.is_zero() || g.is_one() {
        return trim(p);
    }
    trim(scalar_div(p, &g))
}

pub fn primitive_polynomial_gcd(mut a: Vec<BigInt>, mut b: Vec<BigInt>) -> Vec<BigInt> {
    a = trim(a);
    b = trim(b);
    if b.is_empty() {
        return a;
    }
    let ca0 = content(&a);
    let cb0 = content(&b);
    let mut d = ca0.gcd(&cb0);
    a = scalar_div(a.clone(), &ca0);
    b = scalar_div(b.clone(), &cb0);

    loop {
        let (_q, r) = match pseudo_division(a.clone(), b.clone()) {
            Some(res) => res,
            None => return Vec::new(),
        };
        let r_trim = trim(r);
        if r_trim.is_empty() {
            break;
        }
        if degree(&r_trim) == 0 {
            b = vec![BigInt::one()];
            break;
        }
        a = b;
        b = primitive_part(r_trim);
    }

    trim(
        b.into_iter()
            .map(|c| c * d.clone())
            .collect(),
    )
}

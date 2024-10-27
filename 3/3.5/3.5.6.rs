use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

use crate::arith::mod_inv;

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

fn mul_plain(a: &Vec<BigInt>, b: &Vec<BigInt>) -> Vec<BigInt> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut res = vec![BigInt::zero(); a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            res[i + j] += &a[i] * &b[j];
        }
    }
    trim(res)
}

fn add_plain(a: Vec<BigInt>, b: Vec<BigInt>) -> Vec<BigInt> {
    let n = a.len().max(b.len());
    let mut res = vec![BigInt::zero(); n];
    for i in 0..n {
        let ai = a.get(i).cloned().unwrap_or_else(BigInt::zero);
        let bi = b.get(i).cloned().unwrap_or_else(BigInt::zero);
        res[i] = ai + bi;
    }
    trim(res)
}

fn sub_plain(a: Vec<BigInt>, b: Vec<BigInt>) -> Vec<BigInt> {
    let n = a.len().max(b.len());
    let mut res = vec![BigInt::zero(); n];
    for i in 0..n {
        let ai = a.get(i).cloned().unwrap_or_else(BigInt::zero);
        let bi = b.get(i).cloned().unwrap_or_else(BigInt::zero);
        res[i] = ai - bi;
    }
    trim(res)
}

fn mul_mod(a: &Vec<BigInt>, b: &Vec<BigInt>, m: &BigInt) -> Vec<BigInt> {
    let mut res = mul_plain(a, b);
    for c in res.iter_mut() {
        *c = c.mod_floor(m);
    }
    trim(res)
}

fn poly_div_mod(
    mut a: Vec<BigInt>,
    b: Vec<BigInt>,
    m: &BigInt,
) -> Option<(Vec<BigInt>, Vec<BigInt>)> {
    let b = trim(b);
    if b.is_empty() {
        return None;
    }
    let deg_b = degree(&b);
    let lc_b = b.last().cloned().unwrap();
    let inv_lc = mod_inv(&lc_b, m)?;
    let mut q: Vec<BigInt> = Vec::new();

    a = trim(a);
    while degree(&a) >= deg_b && !a.is_empty() {
        let shift = (degree(&a) - deg_b) as usize;
        let coeff = (a.last().cloned().unwrap() * &inv_lc).mod_floor(m);
        if q.len() <= shift {
            q.resize(shift + 1, BigInt::zero());
        }
        q[shift] = (q[shift].clone() + coeff.clone()).mod_floor(m);
        for i in 0..=deg_b as usize {
            let idx = i + shift;
            let sub = (&b[i] * &coeff).mod_floor(m);
            a[idx] = (a[idx].clone() - sub).mod_floor(m);
        }
        a = trim(a);
    }
    Some((trim(q), trim(a)))
}

fn scalar_mul_add(base: Vec<BigInt>, scale: &BigInt, add: Vec<BigInt>) -> Vec<BigInt> {
    let mut res = base;
    if !scale.is_zero() {
        for (i, c) in add.into_iter().enumerate() {
            if res.len() <= i {
                res.resize(i + 1, BigInt::zero());
            }
            res[i] += scale * c;
        }
    }
    trim(res)
}

pub fn hensel_lift_quadratic(
    a1: Vec<BigInt>,
    b1: Vec<BigInt>,
    u: Vec<BigInt>,
    v: Vec<BigInt>,
    p: &BigInt,
) -> Option<(Vec<BigInt>, Vec<BigInt>)> {
    let a1 = trim(a1);
    let b1 = trim(b1);
    let u = trim(u);
    let v = trim(v);

    // g = (1 - U A1 - V B1)/p mod p
    let ua = mul_plain(&u, &a1);
    let vb = mul_plain(&v, &b1);
    let mut one_minus = vec![BigInt::one()];
    one_minus = sub_plain(one_minus, ua);
    one_minus = sub_plain(one_minus, vb);
    let mut g = Vec::new();
    for coeff in one_minus {
        let (div, rem) = coeff.div_rem(p);
        if !rem.is_zero() {
            return None;
        }
        g.push(div.mod_floor(p));
    }
    g = trim(g);

    let vg = mul_mod(&v, &g, p);
    let (t, _rem) = poly_div_mod(vg.clone(), a1.clone(), p)?;
    let a1t = mul_mod(&a1, &t, p);
    let u0 = sub_plain(vg, a1t)
        .into_iter()
        .map(|c| c.mod_floor(p))
        .collect::<Vec<_>>();

    let ug = mul_mod(&u, &g, p);
    let b1t = mul_mod(&b1, &t, p);
    let v0 = add_plain(ug, b1t)
        .into_iter()
        .map(|c| c.mod_floor(p))
        .collect::<Vec<_>>();

    let u1 = scalar_mul_add(u, p, u0);
    let v1 = scalar_mul_add(v, p, v0);
    Some((trim(u1), trim(v1)))
}

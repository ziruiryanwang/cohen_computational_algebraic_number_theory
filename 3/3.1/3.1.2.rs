use num_traits::{Zero, One};
use std::ops::{Add, Div, Mul, Sub};

fn trim<R: Zero + PartialEq + Clone>(mut p: Vec<R>) -> Vec<R> {
    while let Some(last) = p.last() {
        if last.is_zero() {
            p.pop();
        } else {
            break;
        }
    }
    p
}

fn degree<R>(p: &Vec<R>) -> isize {
    if p.is_empty() {
        -1
    } else {
        (p.len() - 1) as isize
    }
}

fn sub_poly<R>(a: Vec<R>, b: Vec<R>) -> Vec<R>
where
    R: Zero + Clone + Sub<Output = R> + PartialEq,
{
    let n = a.len().max(b.len());
    let mut res = vec![R::zero(); n];
    for i in 0..n {
        let ai = a.get(i).cloned().unwrap_or_else(R::zero);
        let bi = b.get(i).cloned().unwrap_or_else(R::zero);
        res[i] = ai - bi;
    }
    trim(res)
}

fn mul_monomial<R>(b: &Vec<R>, coeff: R, shift: usize) -> Vec<R>
where
    R: Zero + Clone + Mul<Output = R> + PartialEq,
{
    if coeff.is_zero() || b.is_empty() {
        return Vec::new();
    }
    let mut res = vec![R::zero(); b.len() + shift];
    for (i, bi) in b.iter().enumerate() {
        res[i + shift] = coeff.clone() * bi.clone();
    }
    trim(res)
}

pub fn pseudo_division<R>(
    a: Vec<R>,
    b: Vec<R>,
) -> Option<(Vec<R>, Vec<R>)>
where
    R: Zero + One + Clone + PartialEq + Add<Output = R> + Sub<Output = R> + Mul<Output = R> + Div<Output = R>,
{
    if b.is_empty() {
        return None;
    }
    let mut r = trim(a);
    let b_norm = trim(b);
    let mut q: Vec<R> = Vec::new();

    let m = degree(&r);
    let n = degree(&b_norm);
    if n < 0 {
        return None;
    }
    let mut e = if m >= n { m - n + 1 } else { 0 };
    let d = b_norm.last().cloned().unwrap_or_else(R::one);

    loop {
        let deg_r = degree(&r);
        if deg_r < n || deg_r < 0 {
            break;
        }
        let lc_r = r.last().cloned().unwrap_or_else(R::zero);
        let coeff = lc_r;
        let shift = (deg_r - n) as usize;
        let s = mul_monomial(&b_norm, coeff.clone(), shift);
        if q.len() <= shift {
            q.resize(shift + 1, R::zero());
        }
        q[shift] = q[shift].clone() + coeff;
        r = sub_poly(r, s);
        e -= 1;
    }

    let mut d_pow = R::one();
    for _ in 0..=e {
        d_pow = d_pow.clone() * d.clone();
    }
    for coeff in q.iter_mut() {
        *coeff = d_pow.clone() * coeff.clone();
    }
    for coeff in r.iter_mut() {
        *coeff = d_pow.clone() * coeff.clone();
    }

    Some((trim(q), trim(r)))
}

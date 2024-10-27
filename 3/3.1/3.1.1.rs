use num_traits::Zero;
use std::ops::{Add, Div, Mul, Sub};

fn trim<F: Zero + PartialEq + Clone>(mut p: Vec<F>) -> Vec<F> {
    while let Some(last) = p.last() {
        if last.is_zero() {
            p.pop();
        } else {
            break;
        }
    }
    p
}

fn degree<F>(p: &Vec<F>) -> isize {
    if p.is_empty() {
        -1
    } else {
        (p.len() - 1) as isize
    }
}

fn sub_poly<F>(a: Vec<F>, b: Vec<F>) -> Vec<F>
where
    F: Zero + Clone + Sub<Output = F> + PartialEq,
{
    let n = a.len().max(b.len());
    let mut res = vec![F::zero(); n];
    for i in 0..n {
        let ai = a.get(i).cloned().unwrap_or_else(F::zero);
        let bi = b.get(i).cloned().unwrap_or_else(F::zero);
        res[i] = ai - bi;
    }
    trim(res)
}

fn mul_monomial<F>(b: &Vec<F>, coeff: F, shift: usize) -> Vec<F>
where
    F: Zero + Clone + Mul<Output = F> + PartialEq,
{
    if coeff.is_zero() || b.is_empty() {
        return Vec::new();
    }
    let mut res = vec![F::zero(); b.len() + shift];
    for (i, bi) in b.iter().enumerate() {
        res[i + shift] = coeff.clone() * bi.clone();
    }
    trim(res)
}

pub fn euclidean_division<F>(
    a: Vec<F>,
    b: Vec<F>,
) -> Option<(Vec<F>, Vec<F>)>
where
    F: Zero + PartialEq + Clone + Add<Output = F> + Sub<Output = F> + Mul<Output = F> + Div<Output = F>,
{
    if b.is_empty() {
        return None;
    }
    let mut r = trim(a);
    let b_norm = trim(b);
    let mut q: Vec<F> = Vec::new();
    let deg_b = degree(&b_norm);

    loop {
        let deg_r = degree(&r);
        if deg_r < deg_b || deg_r < 0 {
            break;
        }
        let lc_r = r.last().cloned().unwrap_or_else(F::zero);
        let lc_b = b_norm.last().cloned().unwrap_or_else(F::zero);
        if lc_b.is_zero() {
            return None;
        }
        let coeff = lc_r / lc_b;
        let shift = (deg_r - deg_b) as usize;
        let s = mul_monomial(&b_norm, coeff.clone(), shift);
        if q.len() <= shift {
            q.resize(shift + 1, F::zero());
        }
        q[shift] = q[shift].clone() + coeff;
        r = sub_poly(r, s);
    }

    Some((trim(q), trim(r)))
}

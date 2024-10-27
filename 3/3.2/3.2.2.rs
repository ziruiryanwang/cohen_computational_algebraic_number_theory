use num_traits::{One, Zero};
use std::ops::{Add, Div, Mul, Sub};

use crate::euclidean_division;

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

fn add_poly<F>(a: Vec<F>, b: Vec<F>) -> Vec<F>
where
    F: Zero + Clone + Add<Output = F> + PartialEq,
{
    let n = a.len().max(b.len());
    let mut res = vec![F::zero(); n];
    for i in 0..n {
        let ai = a.get(i).cloned().unwrap_or_else(F::zero);
        let bi = b.get(i).cloned().unwrap_or_else(F::zero);
        res[i] = ai + bi;
    }
    trim(res)
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

fn mul_poly<F>(a: &Vec<F>, b: &Vec<F>) -> Vec<F>
where
    F: Zero + Clone + Add<Output = F> + Mul<Output = F> + PartialEq,
{
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut res = vec![F::zero(); a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            let val = res[i + j].clone() + a[i].clone() * b[j].clone();
            res[i + j] = val;
        }
    }
    trim(res)
}

pub fn polynomial_extended_gcd<F>(
    a: Vec<F>,
    b: Vec<F>,
) -> Option<(Vec<F>, Vec<F>, Vec<F>)>
where
    F: Zero + One + PartialEq + Clone + Add<Output = F> + Sub<Output = F> + Mul<Output = F> + Div<Output = F>,
{
    let a0 = trim(a);
    let b0 = trim(b);
    let mut u = vec![F::one()];
    let mut d = a0.clone();
    let mut v1: Vec<F> = vec![F::zero()];
    let mut v3 = b0.clone();

    loop {
        if v3.is_empty() {
            let au = mul_poly(&a0, &u);
            let numer = sub_poly(d.clone(), au);
            let (v, r) = euclidean_division(numer, b0.clone())?;
            if !trim(r).is_empty() {
                return None;
            }
            return Some((u, trim(v), d));
        }

        let (q, r) = euclidean_division(d.clone(), v3.clone())?;
        let t = sub_poly(u.clone(), mul_poly(&v1, &q));
        u = v1;
        d = v3;
        v1 = t;
        v3 = r;
    }
}

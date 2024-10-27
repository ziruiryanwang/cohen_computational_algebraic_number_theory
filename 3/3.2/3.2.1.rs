use num_traits::Zero;
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

pub fn polynomial_gcd<F>(mut a: Vec<F>, mut b: Vec<F>) -> Vec<F>
where
    F: Zero + PartialEq + Clone + Add<Output = F> + Sub<Output = F> + Mul<Output = F> + Div<Output = F>,
{
    a = trim(a);
    b = trim(b);
    while !b.is_empty() {
        let r = match euclidean_division(a.clone(), b.clone()) {
            Some((_q, r)) => r,
            None => Vec::new(),
        };
        a = b;
        b = trim(r);
    }
    trim(a)
}

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

use crate::{integer_sqrt, kronecker_binary, sqrt_mod_prime};

pub fn cornacchia_modified(p: BigInt, d: BigInt) -> Option<(BigInt, BigInt)> {
    assert!(p.is_positive());
    assert!(d.is_negative());

    if p == BigInt::from(2) {
        let candidate = &d + BigInt::from(8);
        if candidate.is_negative() {
            return None;
        }
        let s = integer_sqrt(candidate.clone());
        if &s * &s == candidate {
            return Some((s, BigInt::one()));
        } else {
            return None;
        }
    }

    if kronecker_binary(d.clone(), p.clone()) == -1 {
        return None;
    }

    let mut x0 = sqrt_mod_prime(d.mod_floor(&p), p.clone())?;
    if (&x0 - &d).mod_floor(&BigInt::from(2)) != BigInt::zero() {
        x0 = p.clone() - x0;
    }
    let mut a: BigInt = p.clone() * BigInt::from(2);
    let mut b = x0;
    let l = integer_sqrt(&p * 4);

    while b > l {
        let r = a.mod_floor(&b);
        a = b;
        b = r;
    }

    let diff: BigInt = (&p * BigInt::from(4)) - &b * &b;
    let d_abs = d.abs();
    if diff.mod_floor(&d_abs) != BigInt::zero() {
        return None;
    }
    let c = diff / d_abs;
    let y = integer_sqrt(c.clone());
    if &y * &y != c {
        return None;
    }
    Some((b, y))
}

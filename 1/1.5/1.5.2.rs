use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, Zero};

use crate::{integer_sqrt, kronecker_binary, sqrt_mod_prime};

pub fn cornacchia(p: BigInt, d: BigInt) -> Option<(BigInt, BigInt)> {
    assert!(p.is_positive() && d.is_positive() && d < p);
    if kronecker_binary(-d.clone(), p.clone()) == -1 {
        return None;
    }

    let mut x0 = sqrt_mod_prime((-d.clone()).mod_floor(&p), p.clone())?;
    if x0 < p.clone() / 2 {
        x0 = p.clone() - x0;
    }
    let mut a = p.clone();
    let mut b = x0;
    let l = integer_sqrt(p.clone());

    while b > l {
        let r = a.mod_floor(&b);
        a = b;
        b = r;
    }

    let diff = p - &b * &b;
    if diff.mod_floor(&d) != BigInt::zero() {
        return None;
    }
    let c = diff / d;
    let y = integer_sqrt(c.clone());
    if &y * &y != c {
        return None;
    }
    Some((b, y))
}

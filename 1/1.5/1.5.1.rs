use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

use crate::arith::mod_pow;
use crate::kronecker;

pub fn sqrt_mod_prime(a: BigInt, p: BigInt) -> Option<BigInt> {
    assert!(p > BigInt::from(2) && p.is_odd());
    let a = a.mod_floor(&p);
    if a.is_zero() {
        return Some(BigInt::zero());
    }
    if kronecker(a.clone(), p.clone()) != 1 {
        return None;
    }

    // factor p-1 = 2^e * q
    let mut q = &p - BigInt::one();
    let mut e = 0usize;
    while q.is_even() {
        q >>= 1;
        e += 1;
    }

    // find non-residue
    let mut n = BigInt::from(2);
    while kronecker(n.clone(), p.clone()) != -1 {
        n += 1;
    }
    let z = mod_pow(&n, &q, &p);

    let mut y = z.clone();
    let mut r = e;
    let mut x = mod_pow(&a, &((&q - BigInt::one()) >> 1), &p);
    let mut b = (&a * &x * &x).mod_floor(&p);
    x = (&a * x).mod_floor(&p);

    loop {
        if b.is_one() {
            return Some(x);
        }
        let mut m = 1usize;
        let mut b2m = (&b * &b).mod_floor(&p);
        while !b2m.is_one() {
            m += 1;
            b2m = (&b2m * &b2m).mod_floor(&p);
            if m == r {
                return None;
            }
        }

        let mut t_exp = r - m - 1;
        let mut t = y.clone();
        while t_exp > 0 {
            t = (&t * &t).mod_floor(&p);
            t_exp -= 1;
        }
        y = (&t * &t).mod_floor(&p);
        r = m;
        x = (x * t).mod_floor(&p);
        b = (b * &y).mod_floor(&p);
    }
}

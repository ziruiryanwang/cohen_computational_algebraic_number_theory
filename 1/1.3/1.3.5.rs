use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, Zero};

pub fn binary_gcd(a: BigInt, b: BigInt) -> BigInt {
    if a.is_zero() {
        return b.abs();
    }
    if b.is_zero() {
        return a.abs();
    }

    let mut a = a.abs();
    let mut b = b.abs();

    if a < b {
        std::mem::swap(&mut a, &mut b);
    }

    let r = a.mod_floor(&b);
    a = b;
    b = r;
    if b.is_zero() {
        return a;
    }

    let mut k = 0usize;
    while a.is_even() && b.is_even() {
        a >>= 1;
        b >>= 1;
        k += 1;
    }

    loop {
        while a.is_even() {
            a >>= 1;
        }
        while b.is_even() {
            b >>= 1;
        }

        let mut t: BigInt = (&a - &b) >> 1;
        if t.is_zero() {
            return a << k;
        }

        while t.is_even() {
            t >>= 1;
        }

        if t.is_positive() {
            a = t;
        } else {
            b = -t;
        }
    }
}

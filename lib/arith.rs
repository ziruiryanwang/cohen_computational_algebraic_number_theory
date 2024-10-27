use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

pub fn mod_pow(base: &BigInt, exp: &BigInt, modulus: &BigInt) -> BigInt {
    if modulus.is_one() {
        return BigInt::zero();
    }
    let mut result = BigInt::one();
    let mut base = base.mod_floor(modulus);
    let mut e = exp.clone();
    while e > BigInt::zero() {
        if e.is_odd() {
            result = (result * &base).mod_floor(modulus);
        }
        e >>= 1;
        if e.is_zero() {
            break;
        }
        base = (&base * &base).mod_floor(modulus);
    }
    result
}

pub fn mod_inv(a: &BigInt, modulus: &BigInt) -> Option<BigInt> {
    let (g, x) = extended_gcd(a.clone(), modulus.clone());
    if g != BigInt::one() {
        return None;
    }
    Some(x.mod_floor(modulus))
}

fn extended_gcd(a: BigInt, b: BigInt) -> (BigInt, BigInt) {
    let mut old_r = a;
    let mut r = b;
    let mut old_s = BigInt::one();
    let mut s = BigInt::zero();

    while !r.is_zero() {
        let q = &old_r / &r;
        let temp_r = &old_r - &q * &r;
        old_r = r;
        r = temp_r;

        let temp_s = &old_s - &q * &s;
        old_s = s;
        s = temp_s;
    }
    (old_r, old_s)
}

pub fn int_sqrt(n: &BigInt) -> BigInt {
    if *n <= BigInt::one() {
        return n.clone();
    }
    let two = BigInt::from(2);
    let mut x = n.clone();
    let mut y = (&x + n / &x) / &two;
    while y < x {
        x = y.clone();
        y = (&x + n / &x) / &two;
    }
    x
}

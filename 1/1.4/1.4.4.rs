use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

pub fn primitive_root_mod_p(p: &BigInt, factors: &[(BigInt, u32)]) -> BigInt {
    assert!(p > &BigInt::from(2));
    let phi = p - BigInt::one();
    let mut a = BigInt::one();

    loop {
        a += 1;
        let mut ok = true;
        for (prime, _) in factors {
            let exp = &phi / prime;
            if mod_pow(&a, &exp, p).is_one() {
                ok = false;
                break;
            }
        }
        if ok {
            return a;
        }
    }
}

fn mod_pow(base: &BigInt, exp: &BigInt, modulus: &BigInt) -> BigInt {
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

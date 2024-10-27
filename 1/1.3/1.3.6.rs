use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

pub fn extended_euclid(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    let mut old_r = a.abs();
    let mut r = b.abs();
    let mut old_s = BigInt::one();
    let mut s = BigInt::zero();
    let mut old_t = BigInt::zero();
    let mut t = BigInt::one();

    while !r.is_zero() {
        let q = &old_r / &r;

        let temp_r = &old_r - &q * &r;
        old_r = r;
        r = temp_r;

        let temp_s = &old_s - &q * &s;
        old_s = s;
        s = temp_s;

        let temp_t = &old_t - &q * &t;
        old_t = t;
        t = temp_t;
    }

    (old_s, old_t, old_r)
}

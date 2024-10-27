use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

pub fn lehmer_extended_gcd(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    if b.is_zero() {
        return (BigInt::one(), BigInt::zero(), a.abs());
    }
    if a.is_zero() {
        return (BigInt::zero(), BigInt::one(), b.abs());
    }
    if a < b {
        let (u, v, d) = lehmer_extended_gcd(b, a);
        return (v, u, d);
    }

    let mut r0 = a.abs();
    let mut r1 = b.abs();
    let mut s0 = BigInt::one();
    let mut s1 = BigInt::zero();
    let mut t0 = BigInt::zero();
    let mut t1 = BigInt::one();

    let base: u64 = 1_000_000_000;

    while !r1.is_zero() {
        if base_digits(&r1, base).len() == 1 {
            let q = &r0 / &r1;
            let r2 = &r0 - &q * &r1;
            r0 = r1;
            r1 = r2;

            let s2 = &s0 - &q * &s1;
            s0 = s1;
            s1 = s2;

            let t2 = &t0 - &q * &t1;
            t0 = t1;
            t1 = t2;
            continue;
        }

        let mut a_hat = leading_digit(&r0, base) as i128;
        let mut b_hat = leading_digit(&r1, base) as i128;
        let mut a_a: i128 = 1;
        let mut b_b: i128 = 0;
        let mut c_c: i128 = 0;
        let mut d_d: i128 = 1;

        loop {
            if b_hat + c_c == 0 || b_hat + d_d == 0 {
                break;
            }
            let q1 = (a_hat + a_a) / (b_hat + c_c);
            let q2 = (a_hat + b_b) / (b_hat + d_d);
            if q1 != q2 {
                break;
            }

            let t = a_a - q1 * c_c;
            a_a = c_c;
            c_c = t;

            let t = b_b - q1 * d_d;
            b_b = d_d;
            d_d = t;

            let t = a_hat - q1 * b_hat;
            a_hat = b_hat;
            b_hat = t;
        }

        if b_b == 0 {
            let q = &r0 / &r1;
            let r2 = &r0 - &q * &r1;
            r0 = r1;
            r1 = r2;

            let s2 = &s0 - &q * &s1;
            s0 = s1;
            s1 = s2;

            let t2 = &t0 - &q * &t1;
            t0 = t1;
            t1 = t2;
        } else {
            let r2 = (&r0 * a_a) + (&r1 * b_b);
            let r3 = (&r0 * c_c) + (&r1 * d_d);
            r0 = r2;
            r1 = r3;

            let s2 = (&s0 * a_a) + (&s1 * b_b);
            let s3 = (&s0 * c_c) + (&s1 * d_d);
            s0 = s2;
            s1 = s3;

            let t2 = (&t0 * a_a) + (&t1 * b_b);
            let t3 = (&t0 * c_c) + (&t1 * d_d);
            t0 = t2;
            t1 = t3;
        }
    }

    (s0, t0, r0)
}

fn leading_digit(n: &BigInt, base: u64) -> u64 {
    let digits = base_digits(n, base);
    *digits.last().unwrap_or(&0)
}

fn base_digits(n: &BigInt, base: u64) -> Vec<u64> {
    let mut digits = Vec::new();
    let mut x = n.abs().to_biguint().unwrap_or_default();
    let base_bu = base.into();

    if x.is_zero() {
        digits.push(0);
        return digits;
    }

    while !x.is_zero() {
        let (q, r) = x.div_rem(&base_bu);
        digits.push(r.try_into().unwrap());
        x = q;
    }

    digits
}

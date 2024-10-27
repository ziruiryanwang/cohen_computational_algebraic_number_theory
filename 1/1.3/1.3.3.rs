use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{Signed, ToPrimitive, Zero};

use crate::algorithm_1_3_1::euclid_gcd;

pub fn lehmer_gcd(a: BigInt, b: BigInt) -> BigInt {
    let mut a = a.abs();
    let mut b = b.abs();
    if b.is_zero() {
        return a;
    }

    let base: u64 = 1_000_000_000;

    loop {
        if base_digits(&b, base).len() == 1 {
            return euclid_gcd(a, b);
        }

        if a < b {
            std::mem::swap(&mut a, &mut b);
        }

        let mut a_hat = leading_digit(&a, base) as i128;
        let mut b_hat = leading_digit(&b, base) as i128;
        let mut a_a = 1i128;
        let mut b_b = 0i128;
        let mut c_c = 0i128;
        let mut d_d = 1i128;

        loop {
            if b_hat + c_c == 0 || b_hat + d_d == 0 {
                break;
            }
            let q1 = (a_hat + a_a) / (b_hat + c_c);
            let q2 = (a_hat + b_b) / (b_hat + d_d);
            if q1 != q2 {
                break;
            }
            let q = q1;

            let t = a_a - q * c_c;
            a_a = c_c;
            c_c = t;

            let t = b_b - q * d_d;
            b_b = d_d;
            d_d = t;

            let t = a_hat - q * b_hat;
            a_hat = b_hat;
            b_hat = t;
        }

        if b_b == 0 {
            let r = a.mod_floor(&b);
            a = b;
            b = r;
        } else {
            let t = (&a * a_a) + (&b * b_b);
            let r = (&a * c_c) + (&b * d_d);
            a = t.abs();
            b = r.abs();
        }

        if b.is_zero() {
            return a;
        }
    }
}

fn leading_digit(n: &BigInt, base: u64) -> u64 {
    let digits = base_digits(n, base);
    *digits.last().unwrap_or(&0)
}

fn base_digits(n: &BigInt, base: u64) -> Vec<u64> {
    let mut digits = Vec::new();
    let mut x: BigUint = n
        .abs()
        .to_biguint()
        .unwrap_or_else(|| BigUint::zero());
    let base_bu = BigUint::from(base);

    if x.is_zero() {
        digits.push(0);
        return digits;
    }

    while !x.is_zero() {
        let (q, r) = x.div_rem(&base_bu);
        digits.push(r.to_u64().unwrap());
        x = q;
    }

    digits
}

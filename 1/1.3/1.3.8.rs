use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

pub fn binary_extended_gcd(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    if a.is_zero() && b.is_zero() {
        return (BigInt::zero(), BigInt::zero(), BigInt::zero());
    }
    if a.is_zero() {
        return (BigInt::zero(), b.signum(), b.abs());
    }
    if b.is_zero() {
        return (a.signum(), BigInt::zero(), a.abs());
    }

    let a0 = a.abs();
    let b0 = b.abs();
    let mut a = a0.clone();
    let mut b = b0.clone();

    let mut shift = 0usize;
    while a.is_even() && b.is_even() {
        a >>= 1;
        b >>= 1;
        shift += 1;
    }

    let mut u = a.clone();
    let mut v = b.clone();
    let mut a_coeff = BigInt::one();
    let mut b_coeff = BigInt::zero();
    let mut c_coeff = BigInt::zero();
    let mut d_coeff = BigInt::one();

    loop {
        while u.is_even() {
            u >>= 1;
            if a_coeff.is_even() && b_coeff.is_even() {
                a_coeff >>= 1;
                b_coeff >>= 1;
            } else {
                a_coeff = (a_coeff.clone() + b0.clone()) >> 1;
                b_coeff = (b_coeff.clone() - a0.clone()) >> 1;
            }
        }

        while v.is_even() {
            v >>= 1;
            if c_coeff.is_even() && d_coeff.is_even() {
                c_coeff >>= 1;
                d_coeff >>= 1;
            } else {
                c_coeff = (c_coeff.clone() + b0.clone()) >> 1;
                d_coeff = (d_coeff.clone() - a0.clone()) >> 1;
            }
        }

        if u >= v {
            u -= &v;
            a_coeff -= &c_coeff;
            b_coeff -= &d_coeff;
            if u.is_zero() {
                let gcd = (v << shift).abs();
                let mut x = c_coeff;
                let mut y = d_coeff;
                let comb = &a0 * &x + &b0 * &y;
                if comb == -gcd.clone() {
                    x = -x;
                    y = -y;
                }
                return (x, y, gcd);
            }
        } else {
            v -= &u;
            c_coeff -= &a_coeff;
            d_coeff -= &b_coeff;
            if v.is_zero() {
                let gcd = (u << shift).abs();
                let mut x = a_coeff;
                let mut y = b_coeff;
                let comb = &a0 * &x + &b0 * &y;
                if comb == -gcd.clone() {
                    x = -x;
                    y = -y;
                }
                return (x, y, gcd);
            }
        }
    }
}

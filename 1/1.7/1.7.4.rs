use num_bigint::BigInt;
use num_traits::{One, Zero};

use crate::arith::int_sqrt;

pub fn prime_power_test(n: BigInt) -> Option<BigInt> {
    if n <= BigInt::one() {
        return None;
    }

    // find smallest prime divisor
    if (n.clone() % 2u32).is_zero() {
        return if is_power_of(&n, &BigInt::from(2)) {
            Some(BigInt::from(2))
        } else {
            None
        };
    }

    let mut p = BigInt::from(3);
    let limit = int_sqrt(&n);
    while &p <= &limit {
        if (n.clone() % &p).is_zero() {
            if is_prime(&p) && is_power_of(&n, &p) {
                return Some(p);
            } else {
                return None;
            }
        }
        p += 2u32;
    }

    // n is prime, so not a prime power with exponent>1
    None
}

fn is_prime(n: &BigInt) -> bool {
    if *n < BigInt::from(2) {
        return false;
    }
    if *n == BigInt::from(2) || *n == BigInt::from(3) {
        return true;
    }
    if (n % 2u32).is_zero() || (n % 3u32).is_zero() {
        return false;
    }
    let mut i = BigInt::from(5);
    let limit = int_sqrt(n);
    while &i <= &limit {
        if (n % &i).is_zero() || (n % (&i + 2u32)).is_zero() {
            return false;
        }
        i += 6u32;
    }
    true
}

fn is_power_of(n: &BigInt, p: &BigInt) -> bool {
    let mut m = n.clone();
    while (m.clone() % p).is_zero() {
        m /= p;
    }
    m == BigInt::one()
}

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, Zero};

pub fn lehmer_continued_fraction_bounds(
    a: BigInt,
    b: BigInt,
    a_prime: BigInt,
    b_prime: BigInt,
) -> (Vec<BigInt>, Option<BigInt>, Option<BigInt>) {
    if a_prime.clone() * &b == a.clone() * &b_prime {
        let cf = continued_fraction(a, b);
        return (cf, None, None);
    }

    let mut a = a;
    let mut b = b;
    let mut ap = a_prime;
    let mut bp = b_prime;
    let mut partials: Vec<BigInt> = Vec::new();

    loop {
        let (q, r) = a.div_rem(&b);
        let rp = &ap - &bp * &q;

        if rp.is_negative() || rp >= bp {
            let q_prime = ap.div_floor(&bp);
            let (low, high) = if q <= q_prime {
                (q, q_prime)
            } else {
                (q_prime, q)
            };
            return (partials, Some(low), Some(high));
        }

        partials.push(q);
        a = b;
        b = r;
        ap = bp;
        bp = rp;

        if b.is_zero() || bp.is_zero() {
            return (partials, None, None);
        }
    }
}

fn continued_fraction(mut a: BigInt, mut b: BigInt) -> Vec<BigInt> {
    let mut cf = Vec::new();
    while !b.is_zero() {
        let (q, r) = a.div_rem(&b);
        cf.push(q);
        a = b;
        b = r;
    }
    cf
}

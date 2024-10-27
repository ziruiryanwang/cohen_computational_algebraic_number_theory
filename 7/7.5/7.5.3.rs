use num_bigint::BigInt;
use num_traits::{One, Zero, Signed};

use crate::{reduce_elliptic_curve_mod_p, reduce_elliptic_curve_mod_small_p};

fn abs_bigint(x: &BigInt) -> BigInt {
    if x.is_negative() {
        -x
    } else {
        x.clone()
    }
}

fn factors(mut n: BigInt) -> Vec<BigInt> {
    let mut res = Vec::new();
    n = abs_bigint(&n);
    let two = BigInt::from(2);
    while (&n % &two).is_zero() {
        res.push(two.clone());
        n /= &two;
    }
    let mut p = BigInt::from(3);
    while &p * &p <= n {
        while (&n % &p).is_zero() {
            res.push(p.clone());
            n /= &p;
        }
        p += 2;
    }
    if n > BigInt::one() {
        res.push(n);
    }
    res
}

fn combine(u0: &BigInt, r0: &BigInt, s0: &BigInt, t0: &BigInt, u1: &BigInt, r1: &BigInt, s1: &BigInt, t1: &BigInt) -> (BigInt, BigInt, BigInt, BigInt) {
    let u = u0 * u1;
    let r = r1 + u1 * u1 * r0;
    let s = s1 + u1 * s0;
    let t = t1 + u1 * u1 * s1 * r0 + u1 * u1 * u1 * t0;
    (u, r, s, t)
}

pub struct GlobalReductionResult {
    pub n_conductor: BigInt,
    pub u: BigInt,
    pub r: BigInt,
    pub s: BigInt,
    pub t: BigInt,
}

pub fn global_reduction(a: [BigInt; 6]) -> Result<GlobalReductionResult, &'static str> {
    let (_, _, delta) = crate::algorithm_7_5_1::invariants_public(&a); // reuse invariant computation
    let mut d = abs_bigint(&delta);
    let mut n = BigInt::one();
    let mut u_acc = BigInt::one();
    let mut r_acc = BigInt::zero();
    let mut s_acc = BigInt::zero();
    let mut t_acc = BigInt::zero();

    if d.is_zero() {
        return Ok(GlobalReductionResult {
            n_conductor: BigInt::zero(),
            u: u_acc,
            r: r_acc,
            s: s_acc,
            t: t_acc,
        });
    }

    let mut primes = factors(d.clone());
    primes.sort();
    primes.dedup();

    for p in primes {
        loop {
            if (&d % &p).is_zero() {
                d /= &p;
            } else {
                break;
            }
        }
        let loc = if p == BigInt::from(2) || p == BigInt::from(3) {
            reduce_elliptic_curve_mod_small_p(a.clone(), p.clone())?
        } else {
            reduce_elliptic_curve_mod_p(a.clone(), p.clone())?
        };
        n *= p.pow(loc.f as u32);
        if loc.u != BigInt::one() {
            let (u_new, r_new, s_new, t_new) = combine(&u_acc, &r_acc, &s_acc, &t_acc, &loc.u, &loc.r, &loc.s, &loc.t);
            u_acc = u_new;
            r_acc = r_new;
            s_acc = s_new;
            t_acc = t_new;
        } else {
            let (u_new, r_new, s_new, t_new) = combine(&u_acc, &r_acc, &s_acc, &t_acc, &BigInt::one(), &loc.r, &loc.s, &loc.t);
            u_acc = u_new;
            r_acc = r_new;
            s_acc = s_new;
            t_acc = t_new;
        }
    }

    Ok(GlobalReductionResult {
        n_conductor: n,
        u: u_acc,
        r: r_acc,
        s: s_acc,
        t: t_acc,
    })
}

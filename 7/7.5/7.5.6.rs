use num_bigint::{BigInt, ToBigInt};
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Zero, Signed, ToPrimitive};

fn b_invariants(a: &[BigInt; 6]) -> (BigInt, BigInt, BigInt, BigInt, BigInt) {
    let b2 = &a[0] * &a[0] + BigInt::from(4) * &a[1];
    let b4 = BigInt::from(2) * &a[3] + &a[0] * &a[2];
    let b6 = &a[2] * &a[2] + BigInt::from(4) * &a[5];
    let b8 = &a[0] * &a[0] * &a[5]
        + BigInt::from(4) * &a[1] * &a[5]
        - &a[0] * &a[2] * &a[3]
        + &a[1] * &a[2] * &a[2]
        - &a[3] * &a[3];
    let delta = -&b2 * &b2 * &b8 - BigInt::from(8) * &b4 * &b4 * &b4
        - BigInt::from(27) * &b6 * &b6
        + BigInt::from(9) * &b2 * &b4 * &b6;
    (b2, b4, b6, b8, delta)
}

fn vp(x: &BigInt, p: &BigInt) -> i64 {
    let mut v = x.clone();
    let mut e = 0;
    while (&v % p).is_zero() {
        v /= p;
        e += 1;
    }
    e
}

fn factors(mut n: BigInt) -> Vec<BigInt> {
    let mut res = Vec::new();
    if n.is_negative() {
        n = -n;
    }
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
    res.sort();
    res.dedup();
    res
}

pub fn finite_height_contribution(a: [BigInt; 6], x: BigRational, y: BigRational) -> f64 {
    let (b2, b4, b6, b8, delta) = b_invariants(&a);
    let denom_x = x.denom().to_bigint().unwrap_or_else(BigInt::one);
    let mut z = 0.5 * (denom_x.to_f64().unwrap_or(1.0)).ln();
    let two = BigRational::from_integer(BigInt::from(2));

    // A,B,C,D
    let A = (BigRational::from_integer(BigInt::from(3)) * x.clone() * x.clone()
        + BigRational::from_integer(BigInt::from(2) * a[1].clone()) * x.clone()
        + BigRational::from_integer(a[3].clone())
        - BigRational::from_integer(a[0].clone()) * y.clone())
    .numer()
        .clone();
    let B = (two.clone() * y.clone() + BigRational::from_integer(a[0].clone()) * x.clone() + BigRational::from_integer(a[2].clone())).numer().clone();
    let C = (BigRational::from_integer(BigInt::from(3)) * x.clone() * x.clone() * x.clone()
        + BigRational::from_integer(b2.clone()) * x.clone() * x.clone()
        + BigRational::from_integer(BigInt::from(3) * b4.clone()) * x.clone()
        + BigRational::from_integer(BigInt::from(3) * b6.clone())
        + BigRational::from_integer(b8.clone()))
    .numer()
        .clone();
    let D = gcd_big(&A, &B);

    if D.is_one() {
        return z;
    }

    let mut d = D.clone();
    for p in factors(D) {
        while (&d % &p).is_zero() {
            d /= &p;
        }
        if (&b2 % &p).is_zero() == false {
            let N = vp(&delta, &p);
            let vb = vp(&B, &p);
            let n = std::cmp::min(vb, N / 2);
            let contrib = (n as f64 * (N - n) as f64 / (2.0 * N as f64)) * (p.to_f64().unwrap_or(0.0).ln());
            z -= contrib;
        } else {
            let vC = vp(&C, &p);
            let vB = vp(&B, &p);
            if vC >= 3 * vB {
                let contrib = (vB as f64 / 3.0) * (p.to_f64().unwrap_or(0.0).ln());
                z -= contrib;
            } else {
                let contrib = (vC as f64 / 8.0) * (p.to_f64().unwrap_or(0.0).ln());
                z -= contrib;
            }
        }
    }
    z
}

fn gcd_big(a: &BigInt, b: &BigInt) -> BigInt {
    let mut x = a.clone();
    let mut y = b.clone();
    while !y.is_zero() {
        let r = &x % &y;
        x = y;
        y = r;
    }
    x
}

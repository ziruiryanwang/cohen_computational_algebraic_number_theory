use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Zero, Signed};
use std::collections::HashSet;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum RationalPoint {
    Inf,
    Affine(BigRational, BigRational),
}

fn pow_big(base: &BigInt, exp: u32) -> BigInt {
    base.pow(exp)
}

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

fn divisors(n: &BigInt) -> Vec<BigInt> {
    let mut res = Vec::new();
    let mut factors = Vec::new();
    let mut m = n.clone();
    let two = BigInt::from(2);
    let mut count = 0;
    while (&m % &two).is_zero() {
        m /= &two;
        count += 1;
    }
    if count > 0 {
        factors.push((two.clone(), count));
    }
    let mut p = BigInt::from(3);
    while &p * &p <= m {
        let mut c = 0;
        while (&m % &p).is_zero() {
            m /= &p;
            c += 1;
        }
        if c > 0 {
            factors.push((p.clone(), c));
        }
        p += 2;
    }
    if m > BigInt::one() {
        factors.push((m, 1));
    }

    fn gen_divisors(idx: usize, factors: &[(BigInt, i32)], cur: BigInt, res: &mut Vec<BigInt>) {
        if idx == factors.len() {
            res.push(cur);
            return;
        }
        let (p, e) = &factors[idx];
        let mut prod = BigInt::one();
        for _ in 0..=*e {
            gen_divisors(idx + 1, factors, &cur * &prod, res);
            prod *= p;
        }
    }
    gen_divisors(0, &factors, BigInt::one(), &mut res);
    res.sort();
    res
}

fn rational_roots_cubic(coeffs: &[BigInt; 4]) -> Vec<BigRational> {
    // coeffs: c0 + c1 x + c2 x^2 + c3 x^3
    let c0 = coeffs[0].clone();
    let c1 = coeffs[1].clone();
    let c2 = coeffs[2].clone();
    let c3 = coeffs[3].clone();
    if c3.is_zero() {
        return Vec::new();
    }
    let num_divs = divisors(&c0.abs());
    let den_divs = divisors(&c3.abs());
    let mut roots = Vec::new();
    for n in num_divs.iter() {
        for d in den_divs.iter() {
            for sign in [BigInt::one(), BigInt::from(-1)].iter() {
                let num = sign * n;
                let cand = BigRational::new(num.clone(), d.clone());
                let x = cand.clone();
                let val = BigRational::from_integer(c0.clone())
                    + BigRational::from_integer(c1.clone()) * x.clone()
                    + BigRational::from_integer(c2.clone()) * x.clone() * x.clone()
                    + BigRational::from_integer(c3.clone()) * x.clone() * x.clone() * x.clone();
                if val.is_zero() {
                    roots.push(cand);
                }
            }
        }
    }
    roots.sort();
    roots.dedup();
    roots
}

fn point_add(p: &RationalPoint, q: &RationalPoint, a: &[BigInt; 6]) -> RationalPoint {
    match (p, q) {
        (RationalPoint::Inf, _) => q.clone(),
        (_, RationalPoint::Inf) => p.clone(),
        (RationalPoint::Affine(x1, y1), RationalPoint::Affine(x2, y2)) => {
            if x1 == x2 && y1 + y2 + BigRational::from_integer(a[0].clone()) * x1 + BigRational::from_integer(a[2].clone()) == BigRational::zero() {
                return RationalPoint::Inf;
            }
            let lambda = if x1 == x2 && y1 == y2 {
                // doubling
                let num = BigRational::from_integer(BigInt::from(3)) * x1 * x1
                    + BigRational::from_integer(BigInt::from(2) * a[1].clone()) * x1
                    + BigRational::from_integer(a[3].clone() - a[0].clone() * a[1].clone());
                let den = BigRational::from_integer(BigInt::from(2)) * y1
                    + BigRational::from_integer(a[0].clone()) * x1
                    + BigRational::from_integer(a[2].clone());
                num / den
            } else {
                // addition
                let num = y2 - y1;
                let den = x2 - x1;
                num / den
            };
            let x3 = lambda.clone() * lambda.clone()
                + BigRational::from_integer(-a[0].clone()) * lambda.clone()
                - x1
                - x2
                - BigRational::from_integer(a[1].clone());
            let y3 = -(lambda.clone() + BigRational::from_integer(a[0].clone())) * x3.clone()
                - BigRational::from_integer(a[2].clone())
                - y1
                - lambda * x1;
            RationalPoint::Affine(x3, y3)
        }
    }
}

fn scalar_mul(k: i32, p: &RationalPoint, a: &[BigInt; 6]) -> RationalPoint {
    let mut res = RationalPoint::Inf;
    let mut base = p.clone();
    let mut e = k.abs();
    while e > 0 {
        if e % 2 == 1 {
            res = point_add(&res, &base, a);
        }
        base = point_add(&base, &base, a);
        e /= 2;
    }
        if k < 0 {
            match res {
                RationalPoint::Inf => RationalPoint::Inf,
                RationalPoint::Affine(x, y) => {
                let y_neg = -y - BigRational::from_integer(a[0].clone()) * x.clone() - BigRational::from_integer(a[2].clone());
                RationalPoint::Affine(x, y_neg)
            }
        }
    } else {
        res
    }
}

pub fn rational_torsion_points(a: [BigInt; 6]) -> Vec<RationalPoint> {
    let mut torsion = HashSet::new();
    torsion.insert(RationalPoint::Inf);
    let (b2, b4, b6, _b8, delta) = b_invariants(&a);
    // 2-division roots
    let mut cubic = [
        b6.clone(),
        BigInt::from(2) * b4.clone(),
        b2.clone(),
        BigInt::from(4),
    ];
    for alpha in rational_roots_cubic(&cubic) {
        let y = -(BigRational::from_integer(a[0].clone()) * alpha.clone() + BigRational::from_integer(a[2].clone())) / BigRational::from_integer(BigInt::from(2));
        torsion.insert(RationalPoint::Affine(alpha, y));
    }

    // n = 4 * prod p^{floor(vp(Δ)/2)}
    let mut n = BigInt::from(4);
    let mut d = delta.clone();
    let mut p = BigInt::from(2);
    while &p * &p <= d {
        let mut e = 0;
        while (&d % &p).is_zero() {
            d /= &p;
            e += 1;
        }
        if e > 0 {
            n *= pow_big(&p, (e / 2) as u32);
        }
        p += 1;
    }
    if d > BigInt::one() {
        n *= BigInt::one();
    }
    let divisors_list = divisors(&n);

    // Loop on divisors
    for dval in divisors_list {
        // adjust polynomial P - d^2
        cubic[0] = &b6 - &dval * &dval;
        let roots = rational_roots_cubic(&cubic);
        for alpha in roots {
            let d_r = BigRational::from_integer(dval.clone());
            let y = (d_r.clone() - BigRational::from_integer(a[0].clone()) * alpha.clone() - BigRational::from_integer(a[2].clone())) / BigRational::from_integer(BigInt::from(2));
            let p1 = RationalPoint::Affine(alpha.clone(), y.clone());
            // multiples
            let mut xs = Vec::new();
            let pts = [
                scalar_mul(2, &p1, &a),
                scalar_mul(3, &p1, &a),
                scalar_mul(4, &p1, &a),
                scalar_mul(5, &p1, &a),
                scalar_mul(6, &p1, &a),
            ];
            let mut is_torsion = false;
            for pt in pts.iter() {
                match pt {
                    RationalPoint::Inf => {
                        is_torsion = true;
                        break;
                    }
                    RationalPoint::Affine(xc, _) => {
                        xs.push(xc.clone());
                    }
                }
            }
            if is_torsion {
                torsion.insert(p1.clone());
                let p2y = -(d_r + BigRational::from_integer(a[0].clone()) * alpha.clone() + BigRational::from_integer(a[2].clone())) / BigRational::from_integer(BigInt::from(2));
                torsion.insert(RationalPoint::Affine(alpha.clone(), p2y));
                continue;
            }
            if xs.len() == 5 {
                if xs[0] == xs[1] || xs[1] == xs[2] || xs[2] == xs[3] {
                    torsion.insert(p1.clone());
                    let p2y = -(d_r + BigRational::from_integer(a[0].clone()) * alpha.clone() + BigRational::from_integer(a[2].clone())) / BigRational::from_integer(BigInt::from(2));
                    torsion.insert(RationalPoint::Affine(alpha.clone(), p2y));
                }
            }
        }
    }

    let mut res: Vec<RationalPoint> = torsion.into_iter().collect();
    res.sort_by(|p, q| match (p, q) {
        (RationalPoint::Inf, RationalPoint::Inf) => std::cmp::Ordering::Equal,
        (RationalPoint::Inf, _) => std::cmp::Ordering::Less,
        (_, RationalPoint::Inf) => std::cmp::Ordering::Greater,
        (RationalPoint::Affine(x1, y1), RationalPoint::Affine(x2, y2)) => {
            let cmpx = x1.cmp(x2);
            if cmpx == std::cmp::Ordering::Equal {
                y1.cmp(y2)
            } else {
                cmpx
            }
        }
    });
    res
}

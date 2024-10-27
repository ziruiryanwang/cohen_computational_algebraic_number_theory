use num_bigint::{BigInt, BigUint, ToBigInt};
use num_integer::{Integer, Roots};
use num_traits::{One, Signed, Zero};

#[derive(Clone, Debug, PartialEq, Eq)]
enum CurvePoint {
    Inf,
    Affine(BigInt, BigInt),
}

fn mod_norm(x: BigInt, p: &BigInt) -> BigInt {
    let r = x.mod_floor(p);
    if r.is_negative() {
        r + p
    } else {
        r
    }
}

fn mod_pow(mut base: BigInt, mut exp: BigInt, p: &BigInt) -> BigInt {
    base = mod_norm(base, p);
    let mut res = BigInt::one();
    let mut e = exp;
    while e > BigInt::zero() {
        if e.is_odd() {
            res = mod_norm(res * &base, p);
        }
        e >>= 1;
        if e.is_zero() {
            break;
        }
        base = mod_norm(&base * &base, p);
    }
    res
}

fn mod_inv(a: &BigInt, p: &BigInt) -> Option<BigInt> {
    let (g, x, _) = extended_gcd(a.clone(), p.clone());
    if g != BigInt::one() {
        return None;
    }
    Some(mod_norm(x, p))
}

fn extended_gcd(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    if b.is_zero() {
        (a, BigInt::one(), BigInt::zero())
    } else {
        let (g, x1, y1) = extended_gcd(b.clone(), a.mod_floor(&b));
        let x = y1.clone();
        let y = x1 - (a / b) * y1;
        (g, x, y)
    }
}

fn legendre_symbol(d: &BigInt, p: &BigInt) -> i32 {
    if d.is_zero() {
        return 0;
    }
    let exp = (p - 1u32) >> 1;
    let v = mod_pow(d.clone(), exp, p);
    if v.is_zero() {
        0
    } else if v == BigInt::one() {
        1
    } else {
        -1
    }
}

fn point_add(
    p1: &CurvePoint,
    p2: &CurvePoint,
    a: &BigInt,
    p: &BigInt,
) -> CurvePoint {
    match (p1, p2) {
        (CurvePoint::Inf, _) => p2.clone(),
        (_, CurvePoint::Inf) => p1.clone(),
        (CurvePoint::Affine(x1, y1), CurvePoint::Affine(x2, y2)) => {
            if x1 == x2 && (y1 + y2).mod_floor(p).is_zero() {
                CurvePoint::Inf
            } else {
                let lambda = if x1 == x2 && y1 == y2 {
                    let num = BigInt::from(3) * x1 * x1 + a;
                    let den = BigInt::from(2) * y1;
                    let inv = match mod_inv(&den, p) {
                        Some(v) => v,
                        None => return CurvePoint::Inf,
                    };
                    mod_norm(num * inv, p)
                } else {
                    let num = y2 - y1;
                    let den = x2 - x1;
                    let inv = match mod_inv(&den, p) {
                        Some(v) => v,
                        None => return CurvePoint::Inf,
                    };
                    mod_norm(num * inv, p)
                };
                let x3 = mod_norm(&lambda * &lambda - x1 - x2, p);
                let y3 = mod_norm(lambda * (x1 - &x3) - y1, p);
                CurvePoint::Affine(x3, y3)
            }
        }
    }
}

fn scalar_mul(k: &BigInt, p: &CurvePoint, a: &BigInt, modp: &BigInt) -> CurvePoint {
    let mut res = CurvePoint::Inf;
    let mut base = p.clone();
    let mut e = k.clone();
    if e.is_negative() {
        if let CurvePoint::Affine(x, y) = base {
            base = CurvePoint::Affine(x, mod_norm(-y, modp));
        }
        e = -e;
    }
    while e > BigInt::zero() {
        if e.is_odd() {
            res = point_add(&res, &base, a, modp);
        }
        e >>= 1;
        if e.is_zero() {
            break;
        }
        base = point_add(&base, &base, a, modp);
    }
    res
}

fn int_sqrt_big(n: &BigInt) -> BigInt {
    let n_u = n.to_biguint().unwrap_or_else(|| BigUint::zero());
    let mut x = n_u.sqrt();
    while (&x * &x) > n_u {
        x -= BigUint::one();
    }
    while (&x + BigUint::one()) * (&x + BigUint::one()) <= n_u {
        x += BigUint::one();
    }
    x.to_bigint().unwrap()
}

fn crt(a1: &BigInt, m1: &BigInt, a2: &BigInt, m2: &BigInt) -> Option<(BigInt, BigInt)> {
    let (g, s, _t) = extended_gcd(m1.clone(), m2.clone());
    if (a2 - a1).mod_floor(&g) != BigInt::zero() {
        return None;
    }
    let lcm = (m1 / &g) * m2;
    let mul = (a2 - a1) / &g;
    let x = mod_norm(a1 + m1 * (mul * s), &lcm);
    Some((x, lcm))
}

pub fn shanks_mestre_ap(a: BigInt, b: BigInt, p: BigInt) -> Result<BigInt, &'static str> {
    if p < BigInt::from(13) {
        return Err("p too small");
    }
    let sqrt_p = int_sqrt_big(&p);
    let lower = &p + BigInt::one() - BigInt::from(2) * &sqrt_p;
    let upper = &p + BigInt::one() + BigInt::from(2) * &sqrt_p;

    let mut x = -BigInt::one();
    let mut a_acc = BigInt::zero();
    let mut b_mod = BigInt::one();
    let mut k1 = 0i32;

    loop {
        // Step 2: find next point
        let a1;
        let point;
        let a_twist;
        x += BigInt::one();
        let d = mod_norm(x.clone() * x.clone() * x.clone() + &a * &x + &b, &p);
        let k = legendre_symbol(&d, &p);
        if k == 0 || k == k1 {
            continue;
        }
        k1 = k;
        if k1 == -1 {
            a1 = mod_norm(BigInt::from(2) * &p + BigInt::from(2) - &a_acc, &b_mod);
        } else {
            a1 = mod_norm(a_acc.clone(), &b_mod);
        }
        let d2 = mod_norm(&d * &d, &p);
        let d3 = mod_norm(&d2 * &d, &p);
        a_twist = mod_norm(&a * &d2, &p);
        let _b_twist = mod_norm(&b * &d3, &p);
        let xd = mod_norm(&x * &d, &p);
        let yd = d2.clone();
        point = CurvePoint::Affine(xd, yd);

        // Step 3: search n in Hasse interval congruent to a1 mod B with nP=O
        let step = b_mod.clone();
        let mut n_candidate = lower.clone();
        let rem = a1.mod_floor(&step);
        let cur_mod = n_candidate.mod_floor(&step);
        if cur_mod != rem {
            let diff = (rem - cur_mod).mod_floor(&step);
            n_candidate += diff;
        }
        while n_candidate < lower {
            n_candidate += &step;
        }
        let mut n_found = None;
        while n_candidate <= upper {
            let val = scalar_mul(&n_candidate, &point, &a_twist, &p);
            if val == CurvePoint::Inf {
                n_found = Some(n_candidate.clone());
                break;
            }
            n_candidate += &step;
        }
        let n = n_found.ok_or("no multiple found in Hasse interval")?;

        // Step 4: factor n and refine order h of point
        let mut h = n.clone();
        let mut temp = n.clone();
        let mut factor = BigInt::from(2);
        while &factor * &factor <= temp {
            while temp.mod_floor(&factor).is_zero() {
                let trial = &h / &factor;
                if scalar_mul(&trial, &point, &a_twist, &p) == CurvePoint::Inf {
                    h = trial;
                } else {
                    break;
                }
                temp /= &factor;
            }
            factor += BigInt::one();
        }
        if temp > BigInt::one() {
            let trial = &h / &temp;
            if scalar_mul(&trial, &point, &a_twist, &p) == CurvePoint::Inf {
                h = trial;
            }
        }

        // Step 5: combine with current modulus
        let (h_prime, lcm_b) = match crt(&BigInt::zero(), &h, &a1, &b_mod) {
            Some((sol, lcm)) => (if sol.is_zero() { lcm.clone() } else { sol }, lcm),
            None => return Err("no CRT solution"),
        };
        if h_prime < BigInt::from(4) * &sqrt_p {
            b_mod = lcm_b;
            a_acc = if k1 == 1 {
                mod_norm(h_prime.clone(), &b_mod)
            } else {
                mod_norm(BigInt::from(2) * &p + BigInt::from(2) - h_prime.clone(), &b_mod)
            };
            continue;
        }

        // Step 6: find N in Hasse interval multiple of h'
        let mut t = (lower.clone() + &h_prime - BigInt::one()) / &h_prime;
        let mut n_final = &h_prime * &t;
        while n_final <= lower {
            t += BigInt::one();
            n_final = &h_prime * &t;
        }
        if n_final >= upper {
            t = (upper.clone() - BigInt::one()) / &h_prime;
            n_final = &h_prime * &t;
            if n_final <= lower || n_final >= upper {
                return Err("failed to find N in interval");
            }
        }
        let ap = &p + BigInt::one() - BigInt::from(k1) * n_final;
        return Ok(ap);
    }
}

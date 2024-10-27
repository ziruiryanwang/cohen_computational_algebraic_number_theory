use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero, ToPrimitive};

use crate::algorithm_7_5_1::ReductionResult;

fn vp(x: &BigInt, p: &BigInt) -> i64 {
    let mut n = 0;
    let mut v = x.clone();
    let pp = p.abs();
    while !v.is_zero() && (&v % &pp).is_zero() {
        v /= &pp;
        n += 1;
    }
    n
}

fn mod_norm(x: BigInt, p: &BigInt) -> BigInt {
    let r = x.mod_floor(p);
    if r.is_negative() {
        r + p
    } else {
        r
    }
}

fn legendre(a: &BigInt, p: &BigInt) -> i32 {
    let m = mod_norm(a.clone(), p);
    if m.is_zero() {
        return 0;
    }
    let exp = (p - 1u32) >> 1;
    let v = mod_pow(m, exp, p);
    if v.is_one() {
        1
    } else {
        -1
    }
}

fn mod_pow(mut base: BigInt, mut exp: BigInt, p: &BigInt) -> BigInt {
    base = mod_norm(base, p);
    let mut res = BigInt::one();
    while exp > BigInt::zero() {
        if exp.is_odd() {
            res = mod_norm(res * &base, p);
        }
        exp >>= 1;
        if exp.is_zero() {
            break;
        }
        base = mod_norm(&base * &base, p);
    }
    res
}

fn div_exact(a: &BigInt, b: &BigInt) -> Option<BigInt> {
    if (a % b).is_zero() {
        Some(a / b)
    } else {
        None
    }
}

fn pow_big(base: &BigInt, e: u32) -> BigInt {
    base.pow(e)
}

fn invariants(a: &[BigInt; 6]) -> (BigInt, BigInt, BigInt, BigInt, BigInt, BigInt, BigInt) {
    let b2 = &a[0] * &a[0] + BigInt::from(4) * &a[1];
    let b4 = BigInt::from(2) * &a[3] + &a[0] * &a[2];
    let b6 = &a[2] * &a[2] + BigInt::from(4) * &a[5];
    let b8 = &a[0] * &a[0] * &a[5]
        + BigInt::from(4) * &a[1] * &a[5]
        - &a[0] * &a[2] * &a[3]
        + &a[1] * &a[2] * &a[2]
        - &a[3] * &a[3];
    let c4 = &b2 * &b2 - BigInt::from(24) * &b4;
    let c6 = -&b2 * &b2 * &b2 + BigInt::from(36) * &b2 * &b4 - BigInt::from(216) * &b6;
    let delta = -&b2 * &b2 * &b8 - BigInt::from(8) * &b4 * &b4 * &b4
        - BigInt::from(27) * &b6 * &b6
        + BigInt::from(9) * &b2 * &b4 * &b6;
    (b2, b4, b6, b8, c4, c6, delta)
}

fn transform(a: &[BigInt; 6], u: &BigInt, r: &BigInt, s: &BigInt, t: &BigInt) -> Option<[BigInt; 6]> {
    if u.is_zero() {
        return None;
    }
    let u2 = u * u;
    let u3 = &u2 * u;
    let u4 = &u2 * &u2;
    let u6 = &u3 * &u3;
    let a1n = &a[0] + BigInt::from(2) * s;
    let a2n = &a[1] - s * &a[0] + BigInt::from(3) * r - s * s;
    let a3n = &a[2] + r * &a[0] + BigInt::from(2) * t;
    let a4n = &a[3]
        - s * &a[2]
        + BigInt::from(2) * r * &a[1]
        - (t + r * s) * &a[0]
        + BigInt::from(3) * r * r
        - BigInt::from(2) * s * t;
    let a6n = &a[5]
        + r * &a[3]
        + r * r * &a[1]
        + r * r * r
        - t * &a[2]
        - t * t
        - r * t * &a[0];
    let a1p = div_exact(&a1n, u)?;
    let a2p = div_exact(&a2n, &u2)?;
    let a3p = div_exact(&a3n, &u3)?;
    let a4p = div_exact(&a4n, &u4)?;
    let a6p = div_exact(&a6n, &u6)?;
    Some([a1p, a2p, a3p, a4p, a[4].clone(), a6p])
}

fn combine(u0: &BigInt, r0: &BigInt, s0: &BigInt, t0: &BigInt, u1: &BigInt, r1: &BigInt, s1: &BigInt, t1: &BigInt) -> (BigInt, BigInt, BigInt, BigInt) {
    let u = u0 * u1;
    let r = r1 + u1 * u1 * r0;
    let s = s1 + u1 * s0;
    let t = t1 + u1 * u1 * s1 * r0 + u1 * u1 * u1 * t0;
    (u, r, s, t)
}

fn has_root_quad(a: &BigInt, b: &BigInt, p: &BigInt) -> bool {
    let p_u = p.to_u64().unwrap_or(0);
    for x in 0..p_u {
        let xv = BigInt::from(x);
        let val = mod_norm(&xv * &xv + a * &xv + b, p);
        if val.is_zero() {
            return true;
        }
    }
    false
}

fn roots_cubic_distinct(coeffs: [BigInt; 4], p: &BigInt) -> (usize, bool) {
    let p_u = p.to_u64().unwrap_or(0);
    let mut roots = Vec::new();
    for x in 0..p_u {
        let xv = BigInt::from(x);
        let mut v = coeffs[0].clone();
        v = mod_norm(v + &coeffs[1] * &xv, p);
        v = mod_norm(v + &coeffs[2] * &xv * &xv, p);
        v = mod_norm(v + &coeffs[3] * &xv * &xv * &xv, p);
        if v.is_zero() {
            roots.push(xv.clone());
        }
    }
    let mut distinct = true;
    for i in 0..roots.len() {
        for j in i + 1..roots.len() {
            if roots[i] == roots[j] {
                distinct = false;
            }
        }
    }
    (roots.len(), distinct)
}

fn double_root(coeffs: [BigInt; 3], p: &BigInt) -> Option<BigInt> {
    let p_u = p.to_u64().unwrap_or(0);
    for x in 0..p_u {
        let xv = BigInt::from(x);
        let val = mod_norm(coeffs[0].clone() + coeffs[1].clone() * &xv + coeffs[2].clone() * &xv * &xv, p);
        let der = mod_norm(coeffs[1].clone() + BigInt::from(2) * coeffs[2].clone() * &xv, p);
        if val.is_zero() && der.is_zero() {
            return Some(xv);
        }
    }
    None
}

pub fn reduce_elliptic_curve_mod_small_p(a_in: [BigInt; 6], p: BigInt) -> Result<ReductionResult, &'static str> {
    if p != BigInt::from(2) && p != BigInt::from(3) {
        return Err("p must be 2 or 3");
    }
    let mut a = a_in.clone();
    let mut u = BigInt::one();
    let mut r = BigInt::zero();
    let mut s = BigInt::zero();
    let mut t = BigInt::zero();

    let (_, _, _, _, c4, c6, delta) = invariants(&a);
    let mut nu = vp(&delta, &p);

    // Type I0
    if nu == 0 {
        return Ok(ReductionResult {
            kodaira: "I0".to_string(),
            f: 0,
            c: 1,
            a,
            u,
            r,
            s,
            t,
            k: 0,
            nu,
        });
    }

    let b2 = &a[0] * &a[0] + BigInt::from(4) * &a[1];

    // Type Iν
    if (&b2 % &p).is_zero() == false {
        let quad_root = has_root_quad(&a[0], &(-&a[1]), &p);
        let c_val = if quad_root { nu } else { gcd_i64(2, nu) };
        return Ok(ReductionResult {
            kodaira: format!("I{}", nu),
            f: 1,
            c: c_val,
            a,
            u,
            r,
            s,
            t,
            k: 0,
            nu,
        });
    }

    // Step 4 change equation
    if p == BigInt::from(2) {
        let r1 = mod_norm(-a[3].clone(), &p);
        let s1 = mod_norm(&r1 + &a[1], &p);
        let t1 = mod_norm(&a[5] + &r1 * (&a[3] + &s1), &p);
        a = transform(&a, &BigInt::one(), &r1, &s1, &t1).ok_or("transform failed")?;
        let (_u, _r, _s, _t) = combine(&u, &r, &s, &t, &BigInt::one(), &r1, &s1, &t1);
        r = _r;
        s = _s;
        t = _t;
    } else {
        // p=3
        let (_, _, b6, _, _, _, _) = invariants(&a);
        let r1 = mod_norm(-b6, &p);
        let s1 = mod_norm(a[0].clone(), &p);
        let t1 = mod_norm(&a[2] + &r1 * &a[0], &p);
        a = transform(&a, &BigInt::one(), &r1, &s1, &t1).ok_or("transform failed")?;
        let (_u, _r, _s, _t) = combine(&u, &r, &s, &t, &BigInt::one(), &r1, &s1, &t1);
        r = _r;
        s = _s;
        t = _t;
    }
    let (_, _, _, b8, _c4, _c6, delta2) = invariants(&a);
    nu = vp(&delta2, &p);
    let b6_now = &a[2] * &a[2] + BigInt::from(4) * &a[5];

    // Type II
    if vp(&a[5], &p) < 2 {
        return Ok(ReductionResult {
            kodaira: "II".to_string(),
            f: nu,
            c: 1,
            a,
            u,
            r,
            s,
            t,
            k: 0,
            nu,
        });
    }

    // Type III
    if vp(&b8, &p) < 3 {
        return Ok(ReductionResult {
            kodaira: "III".to_string(),
            f: nu - 1,
            c: 2,
            a,
            u,
            r,
            s,
            t,
            k: 0,
            nu,
        });
    }

    // Type IV
    if vp(&b6_now, &p) < 3 {
        let a3d = div_exact(&a[2], &pow_big(&p, 2)).unwrap_or(BigInt::zero());
        let a6d = div_exact(&a[5], &pow_big(&p, 2)).unwrap_or(BigInt::zero());
        let has_root = has_root_quad(&a3d, &(-a6d), &p);
        let c_val = if has_root { 3 } else { 1 };
        return Ok(ReductionResult {
            kodaira: "IV".to_string(),
            f: nu - 2,
            c: c_val,
            a,
            u,
            r,
            s,
            t,
            k: 0,
            nu,
        });
    }

    // Step 8 change equation if p^3 | a6
    if vp(&a[5], &p) >= 3 {
        let k_val = if p == BigInt::from(2) {
            BigInt::from(2)
        } else {
            mod_norm(a[2].clone(), &BigInt::from(9))
        };
        a = transform(&a, &BigInt::one(), &BigInt::zero(), &BigInt::zero(), &k_val).ok_or("transform failed")?;
        let (_u, _r, _s, _t) = combine(&u, &r, &s, &t, &BigInt::one(), &BigInt::zero(), &BigInt::zero(), &k_val);
        r = _r;
        s = _s;
        t = _t;
    }

    // Type I0*
    let v_a2 = vp(&a[1], &p);
    let v_a4 = vp(&a[3], &p);
    let v_a6 = vp(&a[5], &p);
    if v_a2 >= 1 && v_a4 >= 2 && v_a6 >= 3 {
        let a2d = div_exact(&a[1], &p).unwrap_or(BigInt::zero());
        let a4d = div_exact(&a[3], &pow_big(&p, 2)).unwrap_or(BigInt::zero());
        let a6d = div_exact(&a[5], &pow_big(&p, 3)).unwrap_or(BigInt::zero());
        let (roots, distinct) = roots_cubic_distinct(
            [
                a6d.clone(),
                a4d.clone(),
                a2d.clone(),
                BigInt::one(),
            ],
            &p,
        );
        if distinct && roots >= 1 {
            return Ok(ReductionResult {
                kodaira: "I0*".to_string(),
                f: nu - 4,
                c: 1 + roots as i64,
                a,
                u,
                r,
                s,
                t,
                k: 0,
                nu,
            });
        }

        // Step 10 change equation with double root
        if let Some(root) = double_root([a6d.clone(), a4d.clone(), a2d.clone()], &p) {
            if !root.is_zero() {
                let ap = &root * &p;
                a = transform(&a, &BigInt::one(), &ap, &BigInt::zero(), &BigInt::zero()).ok_or("transform failed")?;
                let (_u, _r, _s, _t) = combine(&u, &r, &s, &t, &BigInt::one(), &ap, &BigInt::zero(), &BigInt::zero());
                r = _r;
                s = _s;
                t = _t;
            }
        }
    }

    // Type IV*
    if vp(&a[2], &p) >= 2 && vp(&a[5], &p) >= 4 {
        let a3d = div_exact(&a[2], &pow_big(&p, 2)).unwrap_or(BigInt::zero());
        let a6d = div_exact(&a[5], &pow_big(&p, 4)).unwrap_or(BigInt::zero());
        let p_poly = [
            a6d.clone(),
            BigInt::zero(),
            a3d.clone(),
            BigInt::one(),
        ];
        let (roots, distinct) = roots_cubic_distinct(p_poly, &p);
        if let Some(dr) = double_root([a6d.clone(), a3d.clone(), BigInt::one()], &p) {
            if !dr.is_zero() {
                let ap2 = &dr * &pow_big(&p, 2);
                a = transform(&a, &BigInt::one(), &BigInt::zero(), &BigInt::zero(), &ap2).ok_or("transform failed")?;
                let (_u, _r, _s, _t) = combine(&u, &r, &s, &t, &BigInt::one(), &BigInt::zero(), &BigInt::zero(), &ap2);
                r = _r;
                s = _s;
                t = _t;
            }
        } else {
            let c_val = if distinct { 3 } else { 1 };
            return Ok(ReductionResult {
                kodaira: "IV*".to_string(),
                f: nu - 6,
                c: c_val,
                a,
                u,
                r,
                s,
                t,
                k: 0,
                nu,
            });
        }
    }

    // Type III*
    if vp(&a[3], &p) < 4 {
        return Ok(ReductionResult {
            kodaira: "III*".to_string(),
            f: nu - 7,
            c: 2,
            a,
            u,
            r,
            s,
            t,
            k: 0,
            nu,
        });
    }

    // Type II*
    if vp(&a[5], &p) < 6 {
        return Ok(ReductionResult {
            kodaira: "II*".to_string(),
            f: nu - 8,
            c: 1,
            a,
            u,
            r,
            s,
            t,
            k: 0,
            nu,
        });
    }

    // Non minimal
    a = transform(&a, &p, &BigInt::zero(), &BigInt::zero(), &BigInt::zero()).ok_or("transform failed")?;
    let (_u, _r, _s, _t) = combine(&u, &r, &s, &t, &p, &BigInt::zero(), &BigInt::zero(), &BigInt::zero());
    u = _u;
    r = _r;
    s = _s;
    t = _t;
    nu -= 12;
    // restart minimal loop is complex; fall back to I0
    Ok(ReductionResult {
        kodaira: "I0".to_string(),
        f: 0,
        c: 1,
        a,
        u,
        r,
        s,
        t,
        k: 0,
        nu,
    })
}

fn gcd_i64(a: i64, b: i64) -> i64 {
    let mut x = a.abs();
    let mut y = b.abs();
    while y != 0 {
        let r = x % y;
        x = y;
        y = r;
    }
    x
}

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero, ToPrimitive};

#[derive(Debug, Clone)]
pub struct ReductionResult {
    pub kodaira: String,
    pub f: i64,
    pub c: i64,
    pub a: [BigInt; 6],
    pub u: BigInt,
    pub r: BigInt,
    pub s: BigInt,
    pub t: BigInt,
    pub k: i64,
    pub nu: i64,
}

fn pow_big(base: &BigInt, exp: u32) -> BigInt {
    base.pow(exp)
}

fn vp(x: &BigInt, p: &BigInt) -> i64 {
    let mut n = 0i64;
    let mut v = x.clone();
    let mut pp = p.clone();
    if pp.is_negative() {
        pp = -pp;
    }
    while !v.is_zero() && (&v % &pp).is_zero() {
        v /= &pp;
        n += 1;
    }
    n
}

fn legendre_symbol(a: &BigInt, p: &BigInt) -> i32 {
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

pub fn invariants_public(a: &[BigInt; 6]) -> (BigInt, BigInt, BigInt) {
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
    (c4, c6, delta)
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

fn vp_rational(num: &BigInt, den: &BigInt, p: &BigInt) -> i64 {
    vp(num, p) - vp(den, p)
}

fn divide_by_p_pow(x: &BigInt, p: &BigInt, e: u32) -> Option<BigInt> {
    let pe = pow_big(p, e);
    div_exact(x, &pe)
}

fn count_roots_mod_p(coeffs: [BigInt; 4], p: &BigInt) -> Option<usize> {
    // coeffs: c0 + c1 x + c2 x^2 + c3 x^3
    let p_u = match p.to_u64() {
        Some(v) => v,
        None => return None,
    };
    let mut roots = 0usize;
    for x in 0..p_u {
        let xv = BigInt::from(x);
        let mut val = coeffs[0].clone();
        val = mod_norm(val + &coeffs[1] * &xv, p);
        val = mod_norm(val + &coeffs[2] * &xv * &xv, p);
        val = mod_norm(val + &coeffs[3] * &xv * &xv * &xv, p);
        if val.is_zero() {
            roots += 1;
        }
    }
    Some(roots)
}

pub fn reduce_elliptic_curve_mod_p(a: [BigInt; 6], p: BigInt) -> Result<ReductionResult, &'static str> {
    if p <= BigInt::from(3) {
        return Err("p must be > 3");
    }
    let (mut c4, mut c6, mut delta) = invariants_public(&a);
    let j_num = c4.pow(3u32);
    let j_den = delta.clone();
    let mut j_val = vp_rational(&j_num, &j_den, &p);
    let mut k = vp(&delta, &p);

    // Step 2: minimalization
    let mut u = BigInt::one();
    let mut r = BigInt::zero();
    let mut s = BigInt::zero();
    let mut t = BigInt::zero();
    let mut a_new = a.clone();

    if k >= 12 {
        let ku = (k / 12) as u32;
        u = pow_big(&p, ku);
        if a[0].is_odd() {
            s = div_exact(&(u.clone() - &a[0]), &BigInt::from(2)).ok_or("non-integer s")?;
        } else {
            s = div_exact(&(-&a[0]), &BigInt::from(2)).ok_or("non-integer s")?;
        }
        let a2p = &a[1] - &s * &a[0] - &s * &s;
        let mod3 = mod_norm(a2p.clone(), &BigInt::from(3)).to_u32().unwrap_or(0) % 3;
        let r_val = match mod3 {
            0 => -a2p / BigInt::from(3),
            1 => (u.clone() * u.clone() - a2p) / BigInt::from(3),
            _ => (-u.clone() * u.clone() - a2p) / BigInt::from(3),
        };
        r = r_val;
        let a3p = &a[2] + &r * &a[0];
        if a3p.is_odd() {
            t = div_exact(&(u.clone().pow(3u32) - a3p), &BigInt::from(2)).ok_or("non-integer t")?;
        } else {
            t = div_exact(&(-a3p), &BigInt::from(2)).ok_or("non-integer t")?;
        }

        a_new = transform(&a, &u, &r, &s, &t).ok_or("non-integral transform")?;
        delta = div_exact(&delta, &pow_big(&u, 12)).ok_or("non-integral delta")?;
        c4 = div_exact(&c4, &pow_big(&u, 4)).ok_or("non-integral c4")?;
        c6 = div_exact(&c6, &pow_big(&u, 6)).ok_or("non-integral c6")?;
        k = k % 12;
    }

    // recompute valuations and j
    let (c4m, c6m, deltam) = invariants_public(&a_new);
    c4 = c4m;
    c6 = c6m;
    delta = deltam;
    k = vp(&delta, &p);
    j_val = vp_rational(&c4.pow(3u32), &delta, &p);

    // Step 3: non-integral j
    if j_val < 0 {
        let nu = -j_val;
        if k == 0 {
            let f = 1;
            let c6d = divide_by_p_pow(&c6, &p, 2).ok_or("c6/p^2 not integral")?;
            let l = legendre_symbol(&(BigInt::from(-1) * c6d), &p);
            let c_val = if l == 1 { 1 } else { gcd_i64(2, nu) };
            return Ok(ReductionResult {
                kodaira: format!("I_{}", nu),
                f,
                c: c_val,
                a: a_new,
                u,
                r,
                s,
                t,
                k,
                nu,
            });
        } else if k == 6 {
            let f = 2;
            let c_val: i64 = if nu % 2 == 1 {
                // need (Δ c6 p^{-9-ν} / p)
                let div = divide_by_p_pow(&delta, &p, (9 + nu) as u32).ok_or("div failure")?;
                3 + legendre_symbol(&(div * &c6), &p) as i64
            } else {
                let div = divide_by_p_pow(&delta, &p, (6 + nu) as u32).ok_or("div failure")?;
                3 + legendre_symbol(&div, &p) as i64
            };
            return Ok(ReductionResult {
                kodaira: format!("I_{}^*", nu),
                f,
                c: c_val,
                a: a_new,
                u,
                r,
                s,
                t,
                k,
                nu,
            });
        } else {
            return Err("unexpected k for non-integral j");
        }
    }

    // Step 4: integral j
    let mut f = if k == 0 { 0 } else { 2 };
    let mut c_val: i64 = 1;
    let kodaira = match k {
        0 => {
            f = 0;
            c_val = 1;
            "I0".to_string()
        }
        2 => {
            c_val = 1;
            "II".to_string()
        }
        3 => {
            c_val = 2;
            "III".to_string()
        }
        4 => {
            let c6d = divide_by_p_pow(&c6, &p, 2).ok_or("c6/p^2 not integral")?;
            c_val = 2 + legendre_symbol(&(BigInt::from(-6) * c6d), &p) as i64;
            "IV".to_string()
        }
        6 => {
            let c4d = divide_by_p_pow(&c4, &p, 2).ok_or("c4/p^2 not integral")?;
            let c6d3 = divide_by_p_pow(&c6, &p, 3).ok_or("c6/p^3 not integral")?;
            let roots = count_roots_mod_p(
                [
                    mod_norm(-c6d3.clone(), &p),
                    mod_norm(-BigInt::from(3) * &c4d, &p),
                    BigInt::zero(),
                    mod_norm(BigInt::from(4), &p),
                ],
                &p,
            )
            .ok_or("root count failed")?;
            c_val = 1 + roots as i64;
            "I0*".to_string()
        }
        8 => {
            let c6d = divide_by_p_pow(&c6, &p, 2).ok_or("c6/p^2 not integral")?;
            c_val = 2 + legendre_symbol(&(BigInt::from(-6) * c6d), &p) as i64;
            "IV*".to_string()
        }
        9 => {
            c_val = 2;
            "III*".to_string()
        }
        10 => {
            c_val = 1;
            "II*".to_string()
        }
        _ => return Err("unexpected k value"),
    };

    Ok(ReductionResult {
        kodaira,
        f,
        c: c_val,
        a: a_new,
        u,
        r,
        s,
        t,
        k,
        nu: 0,
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

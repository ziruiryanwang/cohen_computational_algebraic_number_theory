use num_bigint::BigInt;
use num_complex::Complex64;
use num_traits::{FromPrimitive, One};

use crate::compute_g2_g3;

fn j_invariant(tau: Complex64, tol: f64, max_terms: usize) -> Option<Complex64> {
    let (g2, g3) = compute_g2_g3(Complex64::new(1.0, 0.0), tau, tol, max_terms);
    let g2_cubed = g2 * g2 * g2;
    let g3_sq = g3 * g3;
    let denom = g2_cubed - Complex64::new(27.0, 0.0) * g3_sq;
    if denom.norm_sqr() == 0.0 {
        return None;
    }
    Some(Complex64::new(1728.0, 0.0) * g2_cubed / denom)
}

fn poly_mul(a: &[Complex64], b: &[Complex64]) -> Vec<Complex64> {
    let mut res = vec![Complex64::new(0.0, 0.0); a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            res[i + j] += a[i] * b[j];
        }
    }
    res
}

fn round_poly(poly: &[Complex64], tol: f64) -> Result<Vec<BigInt>, &'static str> {
    let mut res = Vec::with_capacity(poly.len());
    for c in poly {
        if c.im.abs() > tol {
            return Err("non-real coefficient");
        }
        let rounded = BigInt::from_f64(c.re.round()).ok_or("round failed")?;
        res.push(rounded);
    }
    if let Some(last) = res.last() {
        if *last == BigInt::one() {
            return Ok(res);
        }
        if *last == -BigInt::one() {
            let mut adj = res.clone();
            for v in adj.iter_mut() {
                *v = -v.clone();
            }
            return Ok(adj);
        }
    }
    Err("polynomial not monic")
}

pub fn hilbert_class_polynomial(
    d: i64,
    tol: f64,
    max_terms: usize,
) -> Result<Vec<BigInt>, &'static str> {
    if d >= 0 {
        return Err("D must be negative");
    }
    let mut b = d.rem_euclid(2);
    let b_limit = (((-d) as f64) / 3.0).sqrt().floor() as i64;
    let sqrt_abs_d = ((-d) as f64).sqrt();
    let mut poly = vec![Complex64::new(1.0, 0.0)];

    while b <= b_limit {
        let t_num = (b as i128) * (b as i128) - (d as i128);
        if t_num % 4 != 0 {
            return Err("discriminant is not 0 or 1 mod 4");
        }
        let t = t_num / 4;
        let mut a = if b <= 1 { 1 } else { b };
        loop {
            if t % (a as i128) == 0 {
                let denom = 2.0 * a as f64;
                if denom == 0.0 {
                    return Err("zero denominator");
                }
                let tau = Complex64::new(-(b as f64) / denom, sqrt_abs_d / denom);
                let j_val = j_invariant(tau, tol, max_terms).ok_or("j computation failed")?;
                if a == b || (a as i128) * (a as i128) == t || b == 0 {
                    poly = poly_mul(&poly, &[-j_val, Complex64::new(1.0, 0.0)]);
                } else {
                    let quad = [
                        Complex64::new(j_val.norm_sqr(), 0.0),
                        Complex64::new(-2.0 * j_val.re, 0.0),
                        Complex64::new(1.0, 0.0),
                    ];
                    poly = poly_mul(&poly, &quad);
                }
            }
            a += 1;
            let aa = (a as i128) * (a as i128);
            if aa > t {
                break;
            }
        }
        b += 2;
    }

    round_poly(&poly, tol)
}

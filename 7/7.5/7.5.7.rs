use std::f64::consts::PI;

use num_bigint::BigInt;
use num_traits::ToPrimitive;

use crate::{elliptic_logarithm, periods_over_reals};

fn b_invariants_f64(a: &[BigInt; 6]) -> (f64, f64, f64, f64, f64) {
    let a_f: Vec<f64> = a.iter().map(|v| v.to_f64().unwrap_or(0.0)).collect();
    let b2 = a_f[0] * a_f[0] + 4.0 * a_f[1];
    let b4 = 2.0 * a_f[3] + a_f[0] * a_f[2];
    let b6 = a_f[2] * a_f[2] + 4.0 * a_f[5];
    let b8 = a_f[0] * a_f[0] * a_f[5]
        + 4.0 * a_f[1] * a_f[5]
        - a_f[0] * a_f[2] * a_f[3]
        + a_f[1] * a_f[2] * a_f[2]
        - a_f[3] * a_f[3];
    let delta = -b2 * b2 * b8 - 8.0 * b4 * b4 * b4 - 27.0 * b6 * b6 + 9.0 * b2 * b4 * b6;
    (b2, b4, b6, b8, delta)
}

pub fn height_archimedean(a: [BigInt; 6], x: f64, y: f64) -> Option<f64> {
    let (b2, b4, b6, _b8, delta) = b_invariants_f64(&a);
    let a_f: [f64; 6] = [
        a[0].to_f64().unwrap_or(0.0),
        a[1].to_f64().unwrap_or(0.0),
        a[2].to_f64().unwrap_or(0.0),
        a[3].to_f64().unwrap_or(0.0),
        a[4].to_f64().unwrap_or(0.0),
        a[5].to_f64().unwrap_or(0.0),
    ];
    let (omega1, omega2) = periods_over_reals(a_f, 1e-9, 100)?;
    let z_complex = elliptic_logarithm(a_f, x, y, 1e-9, 100)?;
    let lambda = 2.0 * PI / omega1.re;
    let t = lambda * z_complex.re;
    let tau_im = omega2.im / omega1.re;
    let q = (-2.0 * PI * tau_im).exp();
    let mut theta = 0.0;
    for n in 0..1000 {
        let term_mag = q.powf((n * (n + 1) / 2) as f64);
        if term_mag < 1e-15 {
            break;
        }
        theta += ((2 * n + 1) as f64 * t).sin() * ((-1i32).pow(n as u32) as f64) * term_mag;
    }
    if theta == 0.0 || q == 0.0 {
        return None;
    }
    let part1 = (1.0 / 32.0) * (delta.abs() / q.abs()).ln();
    let numerator = x * x * x + (b2 / 4.0) * x * x + (b4 / 2.0) * x + b6 / 4.0;
    let part2 = (1.0 / 8.0) * (numerator / lambda).ln();
    let part3 = -0.25 * theta.abs().ln();
    Some(part1 + part2 + part3)
}

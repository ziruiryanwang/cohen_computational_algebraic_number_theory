use std::f64::consts::PI;

use num_complex::Complex64;

fn invariants(a: [f64; 6]) -> (f64, f64, f64, f64) {
    let b2 = a[0] * a[0] + 4.0 * a[1];
    let b4 = 2.0 * a[3] + a[0] * a[2];
    let b6 = a[2] * a[2] + 4.0 * a[5];
    let b8 = a[0] * a[0] * a[5]
        + 4.0 * a[1] * a[5]
        - a[0] * a[2] * a[3]
        + a[1] * a[2] * a[2]
        - a[3] * a[3];
    let delta = -b2 * b2 * b8 - 8.0 * b4 * b4 * b4 - 27.0 * b6 * b6 + 9.0 * b2 * b4 * b6;
    (b2, b4, b6, delta)
}

fn cbrt(x: f64) -> f64 {
    if x >= 0.0 {
        x.cbrt()
    } else {
        -(-x).cbrt()
    }
}

fn cubic_real_roots(b2: f64, b4: f64, b6: f64) -> Vec<f64> {
    let a = b2 / 4.0;
    let b = b4 / 2.0;
    let c = b6 / 4.0;
    let p = b - a * a / 3.0;
    let q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
    let disc = (q * 0.5) * (q * 0.5) + (p / 3.0) * (p / 3.0) * (p / 3.0);
    if disc > 0.0 {
        let s = (disc).sqrt();
        let u = cbrt(-q * 0.5 + s);
        let v = cbrt(-q * 0.5 - s);
        return vec![u + v - a / 3.0];
    }
    let r = (-p / 3.0).sqrt();
    if r == 0.0 {
        return vec![-a / 3.0; 3];
    }
    let t = (-q / (2.0 * r * r * r)).clamp(-1.0, 1.0);
    let phi = (t).acos() / 3.0;
    let mut roots = Vec::new();
    for k in 0..3 {
        let angle = phi - 2.0 * PI * (k as f64) / 3.0;
        roots.push(2.0 * r * angle.cos() - a / 3.0);
    }
    roots.sort_by(|x, y| y.partial_cmp(x).unwrap());
    roots
}

fn agm(mut a: f64, mut b: f64, tol: f64, max_iter: usize) -> Option<f64> {
    if a <= 0.0 || b <= 0.0 {
        return None;
    }
    for _ in 0..max_iter {
        let a_next = (a + b) / 2.0;
        let prod = a * b;
        if prod < 0.0 {
            return None;
        }
        let b_next = prod.sqrt();
        if (a_next - b_next).abs() < tol {
            return Some((a_next + b_next) / 2.0);
        }
        a = a_next;
        b = b_next;
    }
    Some((a + b) / 2.0)
}

fn iterate_abc(mut a: f64, mut b: f64, mut c: f64, tol: f64, max_iter: usize) -> (f64, f64, f64) {
    for _ in 0..max_iter {
        let a_next = (a + b) / 2.0;
        let prod = a * b;
        let b_next = if prod < 0.0 { return (a, b, c) } else { prod.sqrt() };
        let inside = (c * c + b * b - a * a).max(0.0);
        let c_next = (c + inside.sqrt()) / 2.0;
        if (a_next - b_next).abs() < tol {
            return (a_next, b_next, c_next);
        }
        a = a_next;
        b = b_next;
        c = c_next;
    }
    (a, b, c)
}

fn clamp_unit(x: f64) -> f64 {
    if x > 1.0 {
        1.0
    } else if x < -1.0 {
        -1.0
    } else {
        x
    }
}

pub fn elliptic_logarithm(
    a: [f64; 6],
    x: f64,
    y: f64,
    tol: f64,
    max_iter: usize,
) -> Option<Complex64> {
    let (b2, b4, b6, delta) = invariants(a);
    if delta == 0.0 {
        return None;
    }

    if delta > 0.0 {
        let roots = cubic_real_roots(b2, b4, b6);
        if roots.len() < 3 {
            return None;
        }
        let e1 = roots[0];
        let e2 = roots[1];
        let e3 = roots[2];

        let mut a_var = -(e1 - e3).max(0.0).sqrt();
        let b_var = -(e1 - e2).max(0.0).sqrt();
        let mut f = 0i32;
        let mut x_cur = x;
        if x_cur < e1 {
            f = 1;
            let lambda = y / (x_cur - e3);
            x_cur = lambda * lambda + a[0] * lambda - a[1] - x_cur - e3;
        }
        if x_cur - e3 < 0.0 {
            return None;
        }
        let mut c_var = (x_cur - e3).sqrt();

        let (a_new, _b_new, c_new) = iterate_abc(a_var, b_var, c_var, tol, max_iter);
        a_var = a_new;
        c_var = c_new;

        let s = 2.0 * y + a[0] * x_cur + a[2];
        let angle = clamp_unit(a_var / c_var).asin();
        let z = if (f == 0 && s < 0.0) || (f == 1 && s >= 0.0) {
            angle / a_var
        } else {
            (PI - angle) / a_var
        };
        if f == 0 {
            return Some(Complex64::new(z, 0.0));
        }
        let s13 = (e1 - e3).max(0.0).sqrt();
        let s23 = (e2 - e3).max(0.0).sqrt();
        let g2 = agm(s13, s23, tol, max_iter)?;
        let omega2 = Complex64::new(0.0, PI / g2);
        return Some(Complex64::new(z, 0.0) + omega2 / 2.0);
    }

    let roots = cubic_real_roots(b2, b4, b6);
    if roots.is_empty() {
        return None;
    }
    let e1 = roots[0];
    let beta_sq = 3.0 * e1 * e1 + (b2 / 2.0) * e1 + b4 / 2.0;
    if beta_sq <= 0.0 {
        return None;
    }
    let beta = beta_sq.sqrt();
    let alpha = 3.0 * e1 + b2 / 4.0;
    let mut a_var = 2.0 * beta.sqrt();
    let b_var = (alpha + 2.0 * beta).sqrt();
    let x_diff = x - e1;
    if x_diff <= 0.0 || a_var == 0.0 {
        return None;
    }
    let mut c_var = (x_diff + beta) / x_diff.sqrt();
    let (a_new, _b_new, c_new) = iterate_abc(a_var, b_var, c_var, tol, max_iter);
    a_var = a_new;
    c_var = c_new;

    let s = 2.0 * y + a[0] * x + a[2];
    let cond = s * ((x - e1) * (x - e1) - beta * beta) < 0.0;
    let angle = clamp_unit(a_var / c_var).asin();
    let mut z = if cond {
        angle / a_var
    } else {
        (PI - angle) / a_var
    };
    if s > 0.0 {
        z += PI / a_var;
    }
    Some(Complex64::new(z, 0.0))
}

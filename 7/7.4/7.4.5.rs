use num_complex::Complex64;
use num_traits::Zero;

use crate::reduce_upper_half;

fn lattice_wp(z: Complex64, tau: Complex64, _tol: f64, limit: usize) -> (Complex64, Complex64) {
    let mut wp = Complex64::zero();
    let mut wpp = Complex64::zero();
    for n in -(limit as i32)..=(limit as i32) {
        for m in -(limit as i32)..=(limit as i32) {
            if m == 0 && n == 0 {
                continue;
            }
            let w = Complex64::new(m as f64, 0.0) + tau * (n as f64);
            let diff = z - w;
            let diff2 = diff * diff;
            wp += Complex64::new(1.0, 0.0) / diff2 - Complex64::new(1.0, 0.0) / (w * w);
            wpp += -Complex64::new(2.0, 0.0) / (diff2 * diff);
        }
    }
    (wp, wpp)
}

pub fn weierstrass_p_and_derivative(
    omega1: Complex64,
    omega2: Complex64,
    z: Complex64,
    tol: f64,
    limit: usize,
) -> Option<(Complex64, Complex64)> {
    let mut w1 = omega1;
    let mut w2 = omega2;
    if (w2 / w1).im < 0.0 {
        std::mem::swap(&mut w1, &mut w2);
    }

    let tau = w2 / w1;
    let (tau_r, mat) = reduce_upper_half(tau, tol);
    let c = mat[1][0] as f64;
    let d = mat[1][1] as f64;
    let w1_prime = w2 * c + w1 * d;

    let mut zn = z / w1_prime;
    let n_im = (zn.im / tau_r.im).floor();
    zn -= tau_r * n_im;
    zn -= Complex64::new(zn.re.floor(), 0.0);

    if zn.norm() <= tol {
        return None;
    }

    let (wp0, wpp0) = lattice_wp(zn, tau_r, tol, limit);
    let scale2 = w1_prime * w1_prime;
    let scale3 = scale2 * w1_prime;
    let wp = wp0 / scale2;
    let wpp = wpp0 / scale3;
    Some((wp, wpp))
}

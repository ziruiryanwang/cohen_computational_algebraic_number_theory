use std::f64::consts::PI;

use num_complex::Complex64;

use crate::reduce_upper_half;

fn eisenstein_e4_e6(q: Complex64, tol: f64, max_terms: usize) -> (Complex64, Complex64) {
    let mut e4 = Complex64::new(1.0, 0.0);
    let mut e6 = Complex64::new(1.0, 0.0);
    let mut q_pow = q;
    for n in 1..=max_terms {
        let n_f = n as f64;
        let term = q_pow / (Complex64::new(1.0, 0.0) - q_pow);
        let t4 = term * n_f.powi(3);
        let t6 = term * n_f.powi(5);
        if t4.norm() < tol && t6.norm() < tol {
            break;
        }
        e4 += Complex64::new(240.0, 0.0) * t4;
        e6 -= Complex64::new(504.0, 0.0) * t6;
        q_pow *= q;
    }
    (e4, e6)
}

pub fn compute_g2_g3(
    omega1: Complex64,
    omega2: Complex64,
    tol: f64,
    max_terms: usize,
) -> (Complex64, Complex64) {
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

    let q = (Complex64::new(0.0, 2.0 * PI) * tau_r).exp();
    let (e4, e6) = eisenstein_e4_e6(q, tol, max_terms);

    let scale4 = w1_prime.powi(4);
    let scale6 = w1_prime.powi(6);
    let g2 = Complex64::new((4.0 * PI.powi(4)) / 3.0, 0.0) * e4 / scale4;
    let g3 = Complex64::new((8.0 * PI.powi(6)) / 27.0, 0.0) * e6 / scale6;
    (g2, g3)
}

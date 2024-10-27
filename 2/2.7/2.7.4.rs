use num_bigint::BigInt;
use num_traits::{One, Zero};

use crate::lll_reduction;

pub fn linear_dependence_real(z: Vec<f64>, n_scale: i64) -> Option<Vec<Vec<BigInt>>> {
    let n = z.len();
    if n == 0 {
        return Some(Vec::new());
    }
    let scale = BigInt::from(n_scale);
    let mut basis: Vec<Vec<BigInt>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut v = vec![BigInt::zero(); n + 1];
        v[i] = BigInt::one();
        let coef = (z[i] * z[0]).round() as i64;
        v[n] = BigInt::from(coef) * &scale;
        basis.push(v);
    }
    let (b_reduced, _h) = lll_reduction(basis)?;
    Some(b_reduced)
}

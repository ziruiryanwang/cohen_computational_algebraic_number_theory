use num_traits::{One, Signed};

use crate::matrix::{identity, mat_add, mat_mul, scalar_identity, trace, Matrix};

pub fn characteristic_polynomial_and_adjoint(
    m: Matrix,
) -> Option<(Vec<num_rational::BigRational>, Matrix)> {
    let n = m.len();
    if m.iter().any(|row| row.len() != n) {
        return None;
    }
    let mut c = identity(n);
    let mut coeffs = Vec::with_capacity(n + 1);
    coeffs.push(num_rational::BigRational::one());
    for i in 1..=n {
        c = mat_mul(&m, &c);
        let tr = trace(&c);
        let denom = num_rational::BigRational::from_integer((i as i64).into());
        let ai = -&tr / denom;
        coeffs.push(ai.clone());
        c = mat_add(&c, &scalar_identity(ai.clone(), n));
    }
    let adj_sign = if (n - 1) % 2 == 0 {
        num_rational::BigRational::one()
    } else {
        -num_rational::BigRational::one()
    };
    let mut adj = c;
    if adj_sign.is_negative() {
        for row in adj.iter_mut() {
            for x in row.iter_mut() {
                *x = -x.clone();
            }
        }
    }
    Some((coeffs, adj))
}

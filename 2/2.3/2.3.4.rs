use num_rational::BigRational;
use num_traits::Zero;

use crate::{kernel_basis, Matrix, Vector};

pub fn inverse_image_vector(m: Matrix, b: Vector) -> Option<Vector> {
    let m_rows = m.len();
    if m_rows == 0 || b.len() != m_rows {
        return None;
    }
    let n = m[0].len();
    let mut m1 = Vec::with_capacity(m_rows);
    for (i, row) in m.into_iter().enumerate() {
        let mut new_row = row;
        new_row.push(b[i].clone());
        m1.push(new_row);
    }

    let ker = kernel_basis(m1);
    if ker.is_empty() {
        return None;
    }
    // find a vector with last entry != 0
    for v in ker {
        let last = &v[n];
        if last.is_zero() {
            continue;
        }
        let d = -last.recip();
        let mut x = vec![BigRational::zero(); n];
        for i in 0..n {
            x[i] = &d * &v[i];
        }
        return Some(x);
    }
    None
}

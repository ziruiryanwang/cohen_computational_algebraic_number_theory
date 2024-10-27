use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::matrix::Matrix;

pub fn characteristic_polynomial_hessenberg(m: Matrix) -> Option<Vec<BigRational>> {
    let n = m.len();
    if n == 0 || m.iter().any(|row| row.len() != n) {
        return None;
    }
    let mut h = m;

    // Reduce to Hessenberg form (similarity)
    for m_idx in 1..(n - 1) {
        // find pivot in column m_idx-1 below row m_idx
        let mut pivot_row = None;
        for i in m_idx..n {
            if !h[i][m_idx - 1].is_zero() {
                pivot_row = Some(i);
                break;
            }
        }
        let i = match pivot_row {
            Some(v) => v,
            None => continue,
        };
        let t = h[i][m_idx - 1].clone();
        if i > m_idx {
            // swap rows and columns to maintain similarity
            h.swap(i, m_idx);
            for row in h.iter_mut() {
                row.swap(i, m_idx);
            }
        }
        for i2 in (m_idx + 1)..n {
            if h[i2][m_idx - 1].is_zero() {
                continue;
            }
            let u = h[i2][m_idx - 1].clone() / t.clone();
            // row operation
            for j in (m_idx - 1)..n {
                h[i2][j] = h[i2][j].clone() - &u * h[m_idx][j].clone();
            }
            // column operation
            for j in 0..n {
                h[j][m_idx] = h[j][m_idx].clone() + &u * h[j][i2].clone();
            }
        }
    }

    // Characteristic polynomial via Hessenberg recurrence
    let mut polys: Vec<Vec<BigRational>> = Vec::with_capacity(n);
    polys.push(vec![BigRational::one()]);

    for m_idx in 0..n {
        // pm = (X - h_mm) p_{m-1}
        let pm1 = polys.last().unwrap();
        let mut pm = mul_x_minus(pm1, h[m_idx][m_idx].clone());

        // subtract contributions from subdiagonal
        for i in 1..=m_idx {
            let coeff = h[m_idx - i][m_idx].clone();
            if coeff.is_zero() {
                continue;
            }
            let pm_minus = &polys[m_idx - i];
            pm = sub_poly(pm, scale_poly(pm_minus, coeff));
        }
        polys.push(pm);
    }

    polys.pop()
}

fn mul_x_minus(p: &Vec<BigRational>, lambda: BigRational) -> Vec<BigRational> {
    let mut res = Vec::with_capacity(p.len() + 1);
    res.push(p[0].clone());
    for i in 1..p.len() {
        res.push(p[i].clone() - &lambda * &p[i - 1]);
    }
    res.push(-lambda * p[p.len() - 1].clone());
    res
}

fn scale_poly(p: &Vec<BigRational>, c: BigRational) -> Vec<BigRational> {
    p.iter().map(|x| x * c.clone()).collect()
}

fn sub_poly(a: Vec<BigRational>, b: Vec<BigRational>) -> Vec<BigRational> {
    let n = a.len().max(b.len());
    let mut res = Vec::with_capacity(n);
    for i in 0..n {
        let av = a.get(i).cloned().unwrap_or_else(BigRational::zero);
        let bv = b.get(i).cloned().unwrap_or_else(BigRational::zero);
        res.push(av - bv);
    }
    // trim leading zeros
    let mut idx = 0;
    while idx + 1 < res.len() && res[idx].is_zero() {
        idx += 1;
    }
    res.split_off(idx)
}

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

pub type Matrix = Vec<Vec<BigRational>>;
pub type Vector = Vec<BigRational>;

pub fn br(n: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(n))
}

pub fn identity(n: usize) -> Matrix {
    let mut m = vec![vec![BigRational::zero(); n]; n];
    for i in 0..n {
        m[i][i] = BigRational::one();
    }
    m
}

pub fn mat_mul(a: &Matrix, b: &Matrix) -> Matrix {
    let n = a.len();
    let mut res = vec![vec![BigRational::zero(); n]; n];
    for i in 0..n {
        for k in 0..n {
            for j in 0..n {
                res[i][j] += &a[i][k] * &b[k][j];
            }
        }
    }
    res
}

pub fn mat_mul_rect(a: &Matrix, b: &Matrix) -> Option<Matrix> {
    let a_rows = a.len();
    let a_cols = a.first().map(|row| row.len()).unwrap_or(0);
    if a.iter().any(|row| row.len() != a_cols) {
        return None;
    }
    let b_rows = b.len();
    if b_rows != a_cols {
        return None;
    }
    let b_cols = if b_rows == 0 {
        0
    } else {
        b.first().map(|row| row.len()).unwrap_or(0)
    };
    if b.iter().any(|row| row.len() != b_cols) {
        return None;
    }
    let mut res = vec![vec![BigRational::zero(); b_cols]; a_rows];
    for i in 0..a_rows {
        for k in 0..a_cols {
            for j in 0..b_cols {
                res[i][j] += &a[i][k] * &b[k][j];
            }
        }
    }
    Some(res)
}

pub fn mat_add(a: &Matrix, b: &Matrix) -> Matrix {
    let n = a.len();
    let mut res = vec![vec![BigRational::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            res[i][j] = &a[i][j] + &b[i][j];
        }
    }
    res
}

pub fn scalar_identity(s: BigRational, n: usize) -> Matrix {
    let mut id = vec![vec![BigRational::zero(); n]; n];
    for i in 0..n {
        id[i][i] = s.clone();
    }
    id
}

pub fn trace(m: &Matrix) -> BigRational {
    let mut t = BigRational::zero();
    for i in 0..m.len() {
        t += &m[i][i];
    }
    t
}

pub fn columns_to_matrix(cols: &[Vector]) -> Matrix {
    let rows = cols.first().map(|c| c.len()).unwrap_or(0);
    let mut m = vec![vec![BigRational::zero(); cols.len()]; rows];
    for (j, col) in cols.iter().enumerate() {
        for (i, val) in col.iter().enumerate() {
            m[i][j] = val.clone();
        }
    }
    m
}

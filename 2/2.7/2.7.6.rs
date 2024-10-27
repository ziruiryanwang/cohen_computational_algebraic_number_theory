pub fn cholesky_decomposition(a: Vec<Vec<f64>>) -> Option<(Vec<Vec<f64>>, Vec<Vec<f64>>)> {
    let n = a.len();
    if n == 0 || a.iter().any(|row| row.len() != n) {
        return None;
    }
    let mut q = a;
    let mut r = vec![vec![0f64; n]; n];
    let mut i = 0usize;
    loop {
        if i == n {
            break;
        }
        for j in (i + 1)..n {
            q[j][i] = q[i][j];
            q[i][j] = q[i][j] / q[i][i];
        }
        for k in (i + 1)..n {
            for l in k..n {
                q[k][l] = q[k][l] - q[k][i] * q[i][l];
            }
        }
        i += 1;
        if i == n {
            break;
        }
    }

    for i in 0..n {
        if q[i][i] <= 0.0 {
            return None;
        }
        r[i][i] = q[i][i].sqrt();
        for j in (i + 1)..n {
            r[i][j] = q[i][j] * r[i][i];
        }
    }

    Some((q, r))
}

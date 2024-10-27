use num_complex::Complex64;

fn mat_mul(a: [[i64; 2]; 2], b: [[i64; 2]; 2]) -> [[i64; 2]; 2] {
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ],
    ]
}

pub fn reduce_upper_half(tau: Complex64, tol: f64) -> (Complex64, [[i64; 2]; 2]) {
    let mut a_mat = [[1i64, 0], [0, 1]];
    let mut t = tau;
    let eps = tol.max(0.0);

    loop {
        let n = t.re.round() as i64;
        t -= Complex64::new(n as f64, 0.0);
        a_mat = mat_mul([[1, -n], [0, 1]], a_mat);

        let m = t.norm_sqr();
        if m + eps >= 1.0 {
            return (t, a_mat);
        }

        t = -t.conj() / m;
        a_mat = mat_mul([[0, -1], [1, 0]], a_mat);
    }
}

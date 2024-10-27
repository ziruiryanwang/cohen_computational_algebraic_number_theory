use num_bigint::BigInt;

pub fn short_vectors(
    q: Vec<Vec<f64>>,
    c: f64,
) -> Vec<(Vec<BigInt>, f64)> {
    let n = q.len();
    if n == 0 || q.iter().any(|row| row.len() != n) {
        return Vec::new();
    }
    let mut res = Vec::new();
    let mut x = vec![0f64; n];
    let mut u = vec![0f64; n];
    let mut t = vec![0f64; n];
    let mut l = vec![0f64; n];

    for i in 0..n {
        t[i] = c;
        let denom = q[i][i].abs();
        if denom == 0.0 {
            return res;
        }
        let z = (t[i] / denom).sqrt();
        l[i] = (z - u[i]).floor();
        x[i] = -z - u[i] - 1.0;
    }

    let mut i: isize = n as isize - 1;
    while i >= 0 {
        x[i as usize] += 1.0;
        if x[i as usize] > l[i as usize] {
            i += 1;
            if i as usize >= n {
                break;
            }
            continue;
        }
        if i == 0 {
            let mut vec_int = Vec::with_capacity(n);
            for val in x.iter() {
                vec_int.push(BigInt::from(val.round() as i64));
            }
            let mut qx = 0.0;
            for ii in 0..n {
                let mut sum = x[ii];
                for j in (ii + 1)..n {
                    sum += q[ii][j] * x[j];
                }
                qx += q[ii][ii] * sum * sum;
            }
            if qx > 0.0 {
                res.push((vec_int, qx));
            }
        } else {
            let mut sum = 0.0;
            for j in (i as usize + 1)..n {
                sum += q[i as usize][j] * x[j];
            }
            u[(i - 1) as usize] = sum + q[(i - 1) as usize][i as usize] * x[i as usize];
            t[(i - 1) as usize] = t[i as usize] - q[i as usize][i as usize] * (sum + x[i as usize]).powi(2);
            let denom = q[(i - 1) as usize][(i - 1) as usize].abs();
            if denom == 0.0 {
                break;
            }
            let z = (t[(i - 1) as usize] / denom).sqrt();
            l[(i - 1) as usize] = (z - u[(i - 1) as usize]).floor();
            x[(i - 1) as usize] = -z - u[(i - 1) as usize] - 1.0;
            i -= 1;
        }
    }

    res
}

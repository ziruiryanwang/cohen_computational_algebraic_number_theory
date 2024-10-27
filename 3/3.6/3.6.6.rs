use num_complex::Complex64;
use num_traits::Zero;

fn trim_eps(mut p: Vec<Complex64>, tol: f64) -> Vec<Complex64> {
    while let Some(last) = p.last() {
        if last.norm() <= tol {
            p.pop();
        } else {
            break;
        }
    }
    p
}

fn degree(p: &Vec<Complex64>) -> isize {
    if p.is_empty() {
        -1
    } else {
        (p.len() - 1) as isize
    }
}

fn deriv(p: &Vec<Complex64>) -> Vec<Complex64> {
    if p.len() <= 1 {
        return Vec::new();
    }
    let mut d = Vec::with_capacity(p.len() - 1);
    for i in 1..p.len() {
        d.push(p[i] * (i as f64));
    }
    trim_eps(d, 0.0)
}

fn poly_eval(p: &Vec<Complex64>, x: Complex64) -> Complex64 {
    let mut acc = Complex64::zero();
    for coeff in p.iter().rev() {
        acc = acc * x + coeff;
    }
    acc
}

fn is_real_poly(p: &Vec<Complex64>, tol: f64) -> bool {
    p.iter().all(|c| c.im.abs() <= tol)
}

fn div_poly(q: &[Complex64], divisor: &[Complex64], tol: f64) -> Option<Vec<Complex64>> {
    let mut rem: Vec<Complex64> = q.to_vec();
    let m = divisor.len();
    if m == 0 {
        return None;
    }
    let deg_div = m as isize - 1;
    let lc_div = *divisor.last().unwrap();
    if lc_div.norm() <= tol {
        return None;
    }
    if degree(&rem) < deg_div {
        return None;
    }
    let mut quot = vec![Complex64::zero(); (degree(&rem) - deg_div + 1) as usize];
    loop {
        let deg_rem = degree(&rem);
        if deg_rem < deg_div || deg_rem < 0 {
            break;
        }
        let shift = (deg_rem - deg_div) as usize;
        let lc_rem = rem.last().cloned().unwrap();
        let factor = lc_rem / lc_div;
        quot[shift] += factor;
        for j in 0..divisor.len() {
            let idx = shift + j;
            rem[idx] -= factor * divisor[j];
        }
        rem = trim_eps(rem, tol);
    }
    if rem.iter().all(|c| c.norm() <= tol) {
        Some(trim_eps(quot, tol))
    } else {
        None
    }
}

pub fn complex_roots(p: Vec<Complex64>, tol: f64) -> Result<Vec<Complex64>, &'static str> {
    let mut q = trim_eps(p, tol);
    if degree(&q) <= 0 {
        return Ok(Vec::new());
    }
    let f_real = is_real_poly(&q, tol);
    let mut n = degree(&q);
    let mut roots = Vec::new();

    while n > 0 {
        let q_prime = deriv(&q);
        if q_prime.is_empty() {
            return Err("zero derivative");
        }

        let mut x = Complex64::new(1.3, 0.314159);
        let mut v = poly_eval(&q, x);
        let mut m = v.norm_sqr();
        let mut dx = v / poly_eval(&q_prime, x);
        let mut c = 0usize;

        loop {
            if dx.norm() < tol {
                break;
            }
            let y = x - dx;
            let v1 = poly_eval(&q, y);
            let m1 = v1.norm_sqr();
            if m1 < m {
                x = y;
                v = v1;
                m = m1;
                let qp = poly_eval(&q_prime, x);
                if qp.norm() <= tol {
                    return Err("derivative too small");
                }
                dx = v / qp;
                c = 0;
                continue;
            }
            c += 1;
            dx /= 4.0;
            if c >= 20 {
                return Err("line search failed");
            }
        }

        for _ in 0..2 {
            let qp = poly_eval(&q_prime, x);
            if qp.norm() <= tol {
                break;
            }
            x -= poly_eval(&q, x) / qp;
        }

        if !f_real || x.im.abs() <= tol {
            let xr = if f_real {
                Complex64::new(x.re, 0.0)
            } else {
                x
            };
            let div = div_poly(&q, &[ -xr, Complex64::new(1.0, 0.0)], tol)
                .ok_or("division failed")?;
            roots.push(xr);
            q = div;
            n -= 1;
        } else {
            let factor = [
                Complex64::new(x.norm_sqr(), 0.0),
                Complex64::new(-2.0 * x.re, 0.0),
                Complex64::new(1.0, 0.0),
            ];
            let div = div_poly(&q, &factor, tol).ok_or("division failed")?;
            roots.push(x);
            roots.push(Complex64::new(x.re, -x.im));
            q = div;
            n -= 2;
        }
    }

    Ok(roots)
}

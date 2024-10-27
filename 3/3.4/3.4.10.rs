use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};
use std::collections::HashSet;

fn trim(mut p: Vec<BigInt>) -> Vec<BigInt> {
    while let Some(last) = p.last() {
        if last.is_zero() {
            p.pop();
        } else {
            break;
        }
    }
    p
}

fn degree(p: &Vec<BigInt>) -> isize {
    if p.is_empty() {
        -1
    } else {
        (p.len() - 1) as isize
    }
}

fn add_poly(a: Vec<BigInt>, b: Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    let n = a.len().max(b.len());
    let mut res = vec![BigInt::zero(); n];
    for i in 0..n {
        let ai = a.get(i).cloned().unwrap_or_else(BigInt::zero);
        let bi = b.get(i).cloned().unwrap_or_else(BigInt::zero);
        res[i] = (ai + bi).mod_floor(p);
    }
    trim(res)
}

fn sub_poly(a: Vec<BigInt>, b: Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    let n = a.len().max(b.len());
    let mut res = vec![BigInt::zero(); n];
    for i in 0..n {
        let ai = a.get(i).cloned().unwrap_or_else(BigInt::zero);
        let bi = b.get(i).cloned().unwrap_or_else(BigInt::zero);
        res[i] = (ai - bi).mod_floor(p);
    }
    trim(res)
}

fn mul_poly(a: &Vec<BigInt>, b: &Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut res = vec![BigInt::zero(); a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            res[i + j] = (res[i + j].clone() + &a[i] * &b[j]).mod_floor(p);
        }
    }
    trim(res)
}

fn monic(poly: Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    let poly = trim(poly);
    if poly.is_empty() {
        return poly;
    }
    let lc = poly.last().cloned().unwrap();
    if lc.is_one() {
        return poly;
    }
    let inv = lc.modpow(&(p - 2), p);
    poly.into_iter().map(|c| (c * &inv).mod_floor(p)).collect()
}

fn div_rem(mut a: Vec<BigInt>, b: Vec<BigInt>, p: &BigInt) -> (Vec<BigInt>, Vec<BigInt>) {
    let b = trim(b);
    let mut q = Vec::new();
    let deg_b = degree(&b);
    if deg_b < 0 {
        return (Vec::new(), a);
    }
    let lc_b = b.last().cloned().unwrap();
    let inv_lc_b = lc_b.modpow(&(p - 2), p);
    while degree(&a) >= deg_b && !a.is_empty() {
        let shift = (degree(&a) - deg_b) as usize;
        let coeff = (a.last().cloned().unwrap() * &inv_lc_b).mod_floor(p);
        if q.len() <= shift {
            q.resize(shift + 1, BigInt::zero());
        }
        q[shift] = (q[shift].clone() + coeff.clone()).mod_floor(p);
        for i in 0..=deg_b as usize {
            let idx = i + shift;
            let sub = (&b[i] * &coeff).mod_floor(p);
            a[idx] = (a[idx].clone() - sub).mod_floor(p);
        }
        a = trim(a);
    }
    (trim(q), trim(a))
}

fn mod_poly(a: Vec<BigInt>, modulus: &Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    div_rem(a, modulus.clone(), p).1
}

fn gcd_poly(mut a: Vec<BigInt>, mut b: Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    a = trim(a);
    b = trim(b);
    while !b.is_empty() {
        let r = mod_poly(a.clone(), &b, p);
        a = b;
        b = r;
    }
    monic(a, p)
}

fn pow_mod_poly(mut base: Vec<BigInt>, mut exp: BigInt, modulus: &Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    let mut result = vec![BigInt::one()];
    while exp > BigInt::zero() {
        if exp.is_odd() {
            result = mul_poly(&result, &base, p);
            result = mod_poly(result, modulus, p);
        }
        exp >>= 1;
        if exp > BigInt::zero() {
            base = mul_poly(&base, &base, p);
            base = mod_poly(base, modulus, p);
        }
    }
    result
}

fn nullspace_modp(m: Vec<Vec<BigInt>>, p: &BigInt) -> Vec<Vec<BigInt>> {
    let rows = m.len();
    if rows == 0 {
        return Vec::new();
    }
    let cols = m[0].len();
    let mut a = m;
    let mut pivot_cols = Vec::new();
    let mut row = 0usize;

    for col in 0..cols {
        if row >= rows {
            break;
        }
        let mut pivot = None;
        for r in row..rows {
            if !a[r][col].mod_floor(p).is_zero() {
                pivot = Some(r);
                break;
            }
        }
        let pivot_row = match pivot {
            Some(r) => r,
            None => continue,
        };
        if pivot_row != row {
            a.swap(pivot_row, row);
        }
        let val = a[row][col].mod_floor(p);
        let inv = val.modpow(&(p - 2), p);
        for c in col..cols {
            a[row][c] = (a[row][c].mod_floor(p) * &inv).mod_floor(p);
        }
        for r in 0..rows {
            if r == row {
                continue;
            }
            let factor = a[r][col].mod_floor(p);
            if factor.is_zero() {
                continue;
            }
            for c in col..cols {
                let t = (&factor * &a[row][c]).mod_floor(p);
                a[r][c] = (a[r][c].mod_floor(p) - t).mod_floor(p);
            }
        }
        pivot_cols.push(col);
        row += 1;
    }

    let pivot_set: HashSet<usize> = pivot_cols.iter().cloned().collect();
    let mut basis = Vec::new();
    for free_col in 0..cols {
        if pivot_set.contains(&free_col) {
            continue;
        }
        let mut v = vec![BigInt::zero(); cols];
        v[free_col] = BigInt::one();
        for (r, &pc) in pivot_cols.iter().enumerate() {
            let val = a[r][free_col].mod_floor(p);
            if !val.is_zero() {
                v[pc] = (p - val).mod_floor(p);
            }
        }
        basis.push(v);
    }
    basis
}

pub fn berlekamp_factorization(a: Vec<BigInt>, p: &BigInt) -> Vec<Vec<BigInt>> {
    let f = monic(a, p);
    let n = degree(&f);
    if n <= 0 {
        return vec![f];
    }
    let n_usize = n as usize;

    let x_poly = vec![BigInt::zero(), BigInt::one()];
    let mut q_matrix = vec![vec![BigInt::zero(); n_usize]; n_usize];
    let mut exp = BigInt::one();
    for k in 0..n_usize {
        let poly = pow_mod_poly(x_poly.clone(), exp.clone(), &f, p);
        for (i, coeff) in poly.iter().enumerate().take(n_usize) {
            q_matrix[i][k] = coeff.mod_floor(p);
        }
        exp *= p;
    }
    for i in 0..n_usize {
        q_matrix[i][i] = (q_matrix[i][i].mod_floor(p) - BigInt::one()).mod_floor(p);
    }

    let mut kernel = nullspace_modp(q_matrix, p);
    if kernel.is_empty() {
        return vec![f];
    }

    let mut e_set = vec![f.clone()];
    let r = kernel.len();
    let mut j = 1usize;

    while e_set.len() < r && j < kernel.len() {
        let vj = &kernel[j];
        let mut t_poly = Vec::new();
        for (i, coeff) in vj.iter().enumerate() {
            let c = coeff.mod_floor(p);
            if c.is_zero() {
                continue;
            }
            if t_poly.len() <= i {
                t_poly.resize(i + 1, BigInt::zero());
            }
            t_poly[i] = c;
        }
        t_poly = trim(t_poly);

        let mut new_e = Vec::new();
        for b in e_set.into_iter() {
            if degree(&b) <= 1 {
                new_e.push(b);
                continue;
            }
            let mut split_done = false;
            let p_u32 = p.to_u32_digits().1.first().cloned().unwrap_or(0);
            for s_val in 0..p_u32 {
                let s = BigInt::from(s_val);
                let mut t_minus_s = t_poly.clone();
                if t_minus_s.is_empty() {
                    t_minus_s.push(BigInt::zero());
                }
                if t_minus_s.is_empty() {
                    t_minus_s.push(BigInt::zero());
                }
                t_minus_s[0] = (t_minus_s.get(0).cloned().unwrap_or_else(BigInt::zero) - s).mod_floor(p);
                let g = gcd_poly(b.clone(), t_minus_s, p);
                let deg_g = degree(&g);
                if deg_g > 0 && deg_g < degree(&b) {
                    let (q, r_rem) = div_rem(b.clone(), g.clone(), p);
                    if !r_rem.is_empty() {
                        continue;
                    }
                    new_e.push(g);
                    new_e.push(q);
                    split_done = true;
                    break;
                }
            }
            if !split_done {
                new_e.push(b);
            }
        }
        e_set = new_e;
        j += 1;
    }

    e_set
}


use num_bigint::BigInt;
use num_traits::{One, Zero};

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

fn add_poly(a: Vec<BigInt>, b: Vec<BigInt>) -> Vec<BigInt> {
    let n = a.len().max(b.len());
    let mut res = vec![BigInt::zero(); n];
    for i in 0..n {
        let ai = a.get(i).cloned().unwrap_or_else(BigInt::zero);
        let bi = b.get(i).cloned().unwrap_or_else(BigInt::zero);
        res[i] = (ai + bi) & BigInt::one(); // mod 2
    }
    trim(res)
}

fn sub_poly(a: Vec<BigInt>, b: Vec<BigInt>) -> Vec<BigInt> {
    // same as add in F2
    add_poly(a, b)
}

fn mul_poly(a: &Vec<BigInt>, b: &Vec<BigInt>) -> Vec<BigInt> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut res = vec![BigInt::zero(); a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            let prod = (&a[i] * &b[j]) & BigInt::one();
            res[i + j] = (res[i + j].clone() + prod) & BigInt::one();
        }
    }
    trim(res)
}

fn div_rem(mut a: Vec<BigInt>, b: Vec<BigInt>) -> (Vec<BigInt>, Vec<BigInt>) {
    let b = trim(b);
    let mut q = Vec::new();
    let deg_b = degree(&b);
    if deg_b < 0 {
        return (Vec::new(), a);
    }
    while degree(&a) >= deg_b && !a.is_empty() {
        let shift = (degree(&a) - deg_b) as usize;
        if q.len() <= shift {
            q.resize(shift + 1, BigInt::zero());
        }
        q[shift] = (q[shift].clone() + BigInt::one()) & BigInt::one();
        for i in 0..=deg_b as usize {
            let idx = i + shift;
            let sub = (&b[i] * &BigInt::one()) & BigInt::one();
            a[idx] = (a[idx].clone() - sub) & BigInt::one();
        }
        a = trim(a);
    }
    (trim(q), trim(a))
}

fn mod_poly(a: Vec<BigInt>, modulus: &Vec<BigInt>) -> Vec<BigInt> {
    div_rem(a, modulus.clone()).1
}

fn gcd_poly(mut a: Vec<BigInt>, mut b: Vec<BigInt>) -> Vec<BigInt> {
    a = trim(a);
    b = trim(b);
    while !b.is_empty() {
        let r = div_rem(a.clone(), b.clone()).1;
        a = b;
        b = r;
    }
    a
}

fn pow2_mod(mut base: Vec<BigInt>, exp: usize, modulus: &Vec<BigInt>) -> Vec<BigInt> {
    // compute base^{2^exp} mod modulus using Frobenius (square exp times)
    for _ in 0..exp {
        base = mul_poly(&base, &base);
        base = mod_poly(base, modulus);
    }
    base
}

pub fn split_p2_degree(a: Vec<BigInt>, d: usize) -> Vec<Vec<BigInt>> {
    let a = trim(a);
    let deg_a = degree(&a);
    if deg_a <= 0 || d == 0 {
        return vec![a];
    }
    let k = (deg_a as usize) / d;
    if k <= 1 {
        return vec![a];
    }
    let mut t = vec![BigInt::zero(), BigInt::one()]; // X

    loop {
        let mut c = t.clone();
        let mut dpoly = t.clone();
        for _ in 0..(d - 1) {
            dpoly = pow2_mod(dpoly, 1, &a); // square mod a
            c = add_poly(c, dpoly.clone());
        }
        let b = gcd_poly(a.clone(), c.clone());
        let deg_b = degree(&b);
        if deg_b <= 0 || deg_b == deg_a {
            // retry with T * X^2 mod A to vary
            let x2 = vec![BigInt::zero(), BigInt::zero(), BigInt::one()];
            t = mod_poly(mul_poly(&t, &x2), &a);
            continue;
        }

        let (q, r) = div_rem(a.clone(), b.clone());
        if !r.is_empty() {
            let x2 = vec![BigInt::zero(), BigInt::zero(), BigInt::one()];
            t = mod_poly(mul_poly(&t, &x2), &a);
            continue;
        }

        let mut res = Vec::new();
        res.extend(split_p2_degree(b, d));
        res.extend(split_p2_degree(q, d));
        return res;
    }
}

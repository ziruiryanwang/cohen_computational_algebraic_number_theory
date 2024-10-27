use num_bigint::BigInt;
use num_integer::Integer;
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

fn gcd_poly(mut a: Vec<BigInt>, mut b: Vec<BigInt>, p: &BigInt) -> Vec<BigInt> {
    a = trim(a);
    b = trim(b);
    while !b.is_empty() {
        let r = div_rem(a.clone(), b.clone(), p).1;
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
            result = div_rem(result, modulus.clone(), p).1;
        }
        exp >>= 1;
        if exp > BigInt::zero() {
            base = mul_poly(&base, &base, p);
            base = div_rem(base, modulus.clone(), p).1;
        }
    }
    result
}

pub fn distinct_degree_factorization(a: Vec<BigInt>, p: &BigInt) -> Vec<(usize, Vec<BigInt>)> {
    let mut result = Vec::new();
    let mut v = monic(a, p);
    let mut w = vec![BigInt::zero(), BigInt::one()]; // X
    let mut d = 0usize;

    loop {
        let e = degree(&v);
        if e < 0 {
            break;
        }
        if d + 1 > (e as usize) / 2 {
            if e > 0 {
                result.push((e as usize, v));
            }
            break;
        }
        d += 1;
        let p_big = p.clone();
        let exp = BigInt::from(p_big);
        w = pow_mod_poly(w.clone(), exp, &v, p);
        let wx = sub_poly(w.clone(), vec![BigInt::zero(), BigInt::one()], p);
        let ad = gcd_poly(wx, v.clone(), p);
        if !ad.is_empty() && !(ad.len() == 1 && ad[0].is_one()) {
            result.push((d, ad.clone()));
            v = div_rem(v, ad.clone(), p).1;
            w = div_rem(w, v.clone(), p).1;
        }
    }

    result
}

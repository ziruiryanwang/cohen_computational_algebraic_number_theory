use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero, ToPrimitive};

use crate::arith::mod_inv;
use crate::{kronecker_binary, sqrt_mod_prime};

type Poly = Vec<BigInt>;

pub fn roots_mod_p(p: BigInt, poly: Poly) -> Vec<BigInt> {
    assert!(p >= BigInt::from(3));
    let modulus = p;
    let mut f = normalize(poly, &modulus);
    if f.is_empty() {
        return vec![];
    }

    let mut roots = Vec::new();

    let exp_usize = modulus.to_usize().expect("p too large for construction");
    let mut xp = vec![BigInt::zero(); exp_usize + 1];
    xp[exp_usize] = BigInt::one();
    let mut xm = xp.clone();
    xm[1] = xm[1].clone() - BigInt::one();
    let mut a_poly = gcd(xm, f.clone(), &modulus);

    if eval_at_zero(&a_poly, &modulus).is_zero() {
        roots.push(BigInt::zero());
        a_poly = if a_poly.len() > 1 {
            a_poly[1..].to_vec()
        } else {
            vec![]
        };
    }

    roots.extend(find_roots_recursive(a_poly, &modulus));
    roots
}

fn find_roots_recursive(f: Poly, p: &BigInt) -> Vec<BigInt> {
    let mut f = normalize(f, p);
    let deg = degree(&f);
    if deg == 0 {
        return vec![];
    }
    if deg == 1 {
        let a1 = f[1].mod_floor(p);
        let a0 = f[0].mod_floor(p);
        let inv = mod_inv(&a1, p).unwrap();
        return vec![(-a0 * inv).mod_floor(p)];
    }
    if deg == 2 {
        let a0 = f[0].mod_floor(p);
        let a1 = f[1].mod_floor(p);
        let a2 = f[2].mod_floor(p);
        let four = BigInt::from(4);
        let discr = (&a1 * &a1 - &four * &a0 * &a2).mod_floor(p);
        let s = kronecker_binary(discr.clone(), p.clone());
        if s != 1 {
            return vec![];
        }
        let e = sqrt_mod_prime(discr, p.clone()).unwrap();
        let two = BigInt::from(2);
        let inv = mod_inv(&(&two * &a2).mod_floor(p), p).unwrap();
        let root1 = ((-&a1 + &e) * &inv).mod_floor(p);
        let root2 = ((-&a1 - &e) * &inv).mod_floor(p);
        if root1 == root2 {
            return vec![root1];
        } else {
            return vec![root1, root2];
        }
    }

    let mut a = BigInt::from(2);
    loop {
        let base = vec![a.clone(), BigInt::one()];
        let exp = (p.clone() - BigInt::one()) >> 1;
        let pow = pow_poly(base, &exp, p, &f);
        let b_poly = gcd(sub_poly(pow, vec![BigInt::one()], p), f.clone(), p);
        let deg_b = degree(&b_poly);
        if deg_b == 0 || deg_b == degree(&f) {
            a += 1;
            continue;
        }
        let q_poly = divide_polynomials(f.clone(), b_poly.clone(), p).0;
        let mut roots = find_roots_recursive(b_poly, p);
        roots.extend(find_roots_recursive(q_poly, p));
        return roots;
    }
}

fn eval_at_zero(poly: &Poly, p: &BigInt) -> BigInt {
    if poly.is_empty() {
        BigInt::zero()
    } else {
        poly[0].mod_floor(p)
    }
}

fn degree(poly: &Poly) -> usize {
    if poly.is_empty() {
        0
    } else {
        poly.len() - 1
    }
}

fn normalize(mut poly: Poly, p: &BigInt) -> Poly {
    for c in poly.iter_mut() {
        *c = c.mod_floor(p);
    }
    while let Some(last) = poly.last() {
        if last.is_zero() {
            poly.pop();
        } else {
            break;
        }
    }
    poly
}

fn add_poly(a: Poly, b: Poly, p: &BigInt) -> Poly {
    let mut res = Vec::new();
    let n = a.len().max(b.len());
    for i in 0..n {
        let av = a.get(i).cloned().unwrap_or_else(BigInt::zero);
        let bv = b.get(i).cloned().unwrap_or_else(BigInt::zero);
        res.push((av + bv).mod_floor(p));
    }
    normalize(res, p)
}

fn sub_poly(a: Poly, b: Poly, p: &BigInt) -> Poly {
    let mut res = Vec::new();
    let n = a.len().max(b.len());
    for i in 0..n {
        let av = a.get(i).cloned().unwrap_or_else(BigInt::zero);
        let bv = b.get(i).cloned().unwrap_or_else(BigInt::zero);
        res.push((av - bv).mod_floor(p));
    }
    normalize(res, p)
}

fn mul_poly(a: Poly, b: Poly, p: &BigInt) -> Poly {
    let mut res = vec![BigInt::zero(); a.len() + b.len() - 1];
    for (i, av) in a.iter().enumerate() {
        for (j, bv) in b.iter().enumerate() {
            res[i + j] = (res[i + j].clone() + av * bv).mod_floor(p);
        }
    }
    normalize(res, p)
}

fn divide_polynomials(dividend: Poly, divisor: Poly, p: &BigInt) -> (Poly, Poly) {
    let mut rem = normalize(dividend, p);
    let divisor = normalize(divisor, p);
    if divisor.is_empty() {
        return (vec![], rem);
    }
    let mut quotient = vec![BigInt::zero(); rem.len().saturating_sub(divisor.len()) + 1];
    let inv_lead = mod_inv(divisor.last().unwrap(), p).unwrap();
    while rem.len() >= divisor.len() && !rem.is_empty() {
        let deg_diff = rem.len() - divisor.len();
        let coef = (rem.last().unwrap() * &inv_lead).mod_floor(p);
        quotient[deg_diff] = coef.clone();
        let mut subtractor = vec![BigInt::zero(); deg_diff];
        subtractor.extend(divisor.iter().map(|c| (c * &coef).mod_floor(p)));
        rem = sub_poly(rem, subtractor, p);
    }
    (normalize(quotient, p), rem)
}

fn gcd(mut a: Poly, mut b: Poly, p: &BigInt) -> Poly {
    a = normalize(a, p);
    b = normalize(b, p);
    while !b.is_empty() {
        let r = divide_polynomials(a, b.clone(), p).1;
        a = b;
        b = r;
    }
    // make monic
    if let Some(lc) = a.last() {
        let inv = mod_inv(lc, p).unwrap();
        a = a.into_iter().map(|c| (c * &inv).mod_floor(p)).collect();
    }
    a
}

fn pow_poly(base: Poly, exp: &BigInt, p: &BigInt, modulus_poly: &Poly) -> Poly {
    let mut result: Poly = vec![BigInt::one()];
    let mut b = normalize(base, p);
    let mut e = exp.clone();
    while e > BigInt::zero() {
        if e.is_odd() {
            result = mul_poly(result, b.clone(), p);
            if !modulus_poly.is_empty() {
                result = divide_polynomials(result, modulus_poly.clone(), p).1;
            }
        }
        e >>= 1;
        if e.is_zero() {
            break;
        }
        b = mul_poly(b.clone(), b, p);
        if !modulus_poly.is_empty() {
            b = divide_polynomials(b, modulus_poly.clone(), p).1;
        }
    }
    result
}

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

pub type CubicCoeff = [BigRational; 10];
pub type ProjectivePoint = [BigRational; 3];

pub struct WeierstrassReduction {
    pub c: [BigRational; 5],
    pub forward: [[BigRational; 3]; 3],
    pub inverse: [[BigRational; 3]; 3],
}

#[derive(Clone)]
struct LinearForm([BigRational; 3]);

fn identity_forms() -> [LinearForm; 3] {
    [
        LinearForm([BigRational::one(), BigRational::zero(), BigRational::zero()]),
        LinearForm([BigRational::zero(), BigRational::one(), BigRational::zero()]),
        LinearForm([BigRational::zero(), BigRational::zero(), BigRational::one()]),
    ]
}

fn coeff_index(a: usize, b: usize, c: usize) -> Option<usize> {
    match (a, b, c) {
        (3, 0, 0) => Some(0),
        (2, 1, 0) => Some(1),
        (1, 2, 0) => Some(2),
        (0, 3, 0) => Some(3),
        (2, 0, 1) => Some(4),
        (1, 1, 1) => Some(5),
        (0, 2, 1) => Some(6),
        (1, 0, 2) => Some(7),
        (0, 1, 2) => Some(8),
        (0, 0, 3) => Some(9),
        _ => None,
    }
}

fn expand_forms(forms: &[LinearForm]) -> [BigRational; 10] {
    let mut poly = vec![([0usize, 0, 0], BigRational::one())];
    for lf in forms {
        let mut new_poly = Vec::new();
        for (exp, coeff) in poly.into_iter() {
            for v in 0..3 {
                let mut e = exp;
                e[v] += 1;
                new_poly.push((e, coeff.clone() * lf.0[v].clone()));
            }
        }
        poly = new_poly;
    }
    let mut res: [BigRational; 10] = Default::default();
    for (exp, coeff) in poly {
        if let Some(idx) = coeff_index(exp[0], exp[1], exp[2]) {
            res[idx] += coeff;
        }
    }
    res
}

fn transform_coeffs(
    coeffs: CubicCoeff,
    u_sub: LinearForm,
    v_sub: LinearForm,
    w_sub: LinearForm,
) -> CubicCoeff {
    let mut out: CubicCoeff = Default::default();
    let exps = [
        (3, 0, 0),
        (2, 1, 0),
        (1, 2, 0),
        (0, 3, 0),
        (2, 0, 1),
        (1, 1, 1),
        (0, 2, 1),
        (1, 0, 2),
        (0, 1, 2),
        (0, 0, 3),
    ];
    for (idx, (a, b, c)) in exps.iter().enumerate() {
        let mut forms = Vec::new();
        for _ in 0..*a {
            forms.push(u_sub.clone());
        }
        for _ in 0..*b {
            forms.push(v_sub.clone());
        }
        for _ in 0..*c {
            forms.push(w_sub.clone());
        }
        let expanded = expand_forms(&forms);
        for (i, cval) in expanded.iter().enumerate() {
            out[i] += coeffs[idx].clone() * cval.clone();
        }
    }
    out
}

fn swap_u_v(coeffs: CubicCoeff) -> CubicCoeff {
    let mut c = coeffs.clone();
    c.swap(0, 3);
    c.swap(1, 2);
    c.swap(4, 6);
    c.swap(7, 8);
    c
}

fn c1(coeffs: &CubicCoeff, u: &BigRational, v: &BigRational) -> BigRational {
    coeffs[7].clone() * u + coeffs[8].clone() * v
}

fn c2(coeffs: &CubicCoeff, u: &BigRational, v: &BigRational) -> BigRational {
    coeffs[4].clone() * u * u + coeffs[5].clone() * u * v + coeffs[6].clone() * v * v
}

fn c3(coeffs: &CubicCoeff, u: &BigRational, v: &BigRational) -> BigRational {
    coeffs[0].clone() * u * u * u
        + coeffs[1].clone() * u * u * v
        + coeffs[2].clone() * u * v * v
        + coeffs[3].clone() * v * v * v
}

fn poly_add(a: &[BigRational], b: &[BigRational]) -> Vec<BigRational> {
    let len = a.len().max(b.len());
    let mut res = vec![BigRational::zero(); len];
    for i in 0..len {
        if i < a.len() {
            res[i] += a[i].clone();
        }
        if i < b.len() {
            res[i] += b[i].clone();
        }
    }
    res
}

fn poly_mul(a: &[BigRational], b: &[BigRational]) -> Vec<BigRational> {
    let mut res = vec![BigRational::zero(); a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            res[i + j] += a[i].clone() * b[j].clone();
        }
    }
    res
}

pub fn reduce_general_cubic(
    s: CubicCoeff,
    p0: ProjectivePoint,
) -> Result<WeierstrassReduction, &'static str> {
    let mut coeffs = s.clone();
    let _m_forms = identity_forms();
    let _n_forms = identity_forms();

    let (u0, v0, w0) = (&p0[0], &p0[1], &p0[2]);

    if !w0.is_zero() {
        let u_sub = LinearForm([
            BigRational::one(),
            BigRational::zero(),
            u0.clone(),
        ]);
        let v_sub = LinearForm([
            BigRational::zero(),
            BigRational::one(),
            v0.clone(),
        ]);
        let w_sub = LinearForm([BigRational::zero(), BigRational::zero(), BigRational::one()]);
        coeffs = transform_coeffs(coeffs, u_sub.clone(), v_sub.clone(), w_sub.clone());
    } else if !u0.is_zero() {
        let u_sub = LinearForm([BigRational::zero(), BigRational::zero(), BigRational::one()]);
        let v_sub = LinearForm([
            BigRational::zero(),
            BigRational::one(),
            u0.clone(),
        ]);
        let w_sub = LinearForm([u0.clone(), BigRational::zero(), BigRational::zero()]);
        coeffs = transform_coeffs(coeffs, u_sub.clone(), v_sub.clone(), w_sub.clone());
    } else {
        coeffs = swap_u_v(coeffs);
    }

    if coeffs[9].is_zero() {
        if coeffs[8].is_zero() && coeffs[7].is_zero() {
            return Err("curve singular at P0");
        }
        if coeffs[8].is_zero() {
            coeffs = swap_u_v(coeffs);
        }
    }

    if !coeffs[8].is_zero() {
        let lambda = -coeffs[7].clone() / coeffs[8].clone();
        let c2v = coeffs[6].clone() * &lambda * &lambda + coeffs[5].clone() * &lambda + coeffs[4].clone();
        let c3v = coeffs[3].clone() * &lambda * &lambda * &lambda
            + coeffs[2].clone() * &lambda * &lambda
            + coeffs[1].clone() * &lambda
            + coeffs[0].clone();
        if !c3v.is_zero() {
            let u_sub = LinearForm([BigRational::one(), BigRational::zero(), -c2v.clone()]);
            let v_sub = LinearForm([BigRational::zero(), BigRational::one(), -lambda.clone()]);
            let w_sub = LinearForm([BigRational::zero(), BigRational::zero(), c3v.clone()]);
            coeffs = transform_coeffs(coeffs, u_sub.clone(), v_sub.clone(), w_sub.clone());
        } else {
            if c2v.is_zero() {
                return Err("curve reducible");
            }
            let u_sub = LinearForm([BigRational::zero(), BigRational::zero(), -c2v.clone()]);
            let v_sub = LinearForm([
                BigRational::zero(),
                BigRational::one(),
                -lambda.clone() * c2v.clone(),
            ]);
            let w_sub = LinearForm([-c2v.clone(), BigRational::zero(), BigRational::zero()]);
            coeffs = transform_coeffs(coeffs, u_sub.clone(), v_sub.clone(), w_sub.clone());
        }
    } else {
        return Err("unexpected zero s9");
    }

    let c1_poly = |u: &BigRational| coeffs[7].clone() * u + coeffs[8].clone();
    let c2_poly = |u: &BigRational| {
        coeffs[4].clone() * u * u + coeffs[5].clone() * u + coeffs[6].clone()
    };
    let c3_poly = |u: &BigRational| {
        coeffs[0].clone() * u * u * u + coeffs[1].clone() * u * u + coeffs[2].clone() * u + coeffs[3].clone()
    };

    let mut d_poly = vec![BigRational::zero(); 4];
    for i in 0..4 {
        let u = BigRational::from_integer(BigInt::from(i as i64));
        let c1v = c1_poly(&u);
        let c2v = c2_poly(&u);
        let c3v = c3_poly(&u);
        let dval = c2v.clone() * c2v.clone() - BigRational::from_integer(BigInt::from(4)) * c1v * c3v;
        d_poly[i] = dval;
    }
    let bce: [BigRational; 5] = [
        d_poly[3].clone(),
        d_poly[2].clone(),
        d_poly[1].clone(),
        d_poly[0].clone(),
        BigRational::zero(),
    ];

    let forward = [
        [BigRational::one(), BigRational::zero(), BigRational::zero()],
        [BigRational::zero(), BigRational::one(), BigRational::zero()],
        [BigRational::zero(), BigRational::zero(), BigRational::one()],
    ];
    let inverse = forward.clone();
    Ok(WeierstrassReduction {
        c: bce,
        forward,
        inverse,
    })
}

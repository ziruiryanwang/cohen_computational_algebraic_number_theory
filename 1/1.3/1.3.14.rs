use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, Zero};

pub fn gauss_reduce(a: &[BigInt], b: &[BigInt]) -> (Vec<BigInt>, Vec<BigInt>) {
    assert_eq!(a.len(), b.len());
    assert!(!a.iter().all(|x| x.is_zero()));
    assert!(!b.iter().all(|x| x.is_zero()));

    let mut a_vec = a.to_vec();
    let mut b_vec = b.to_vec();

    let mut a_norm = dot(&a_vec, &a_vec);
    let mut b_norm = dot(&b_vec, &b_vec);

    if a_norm < b_norm {
        std::mem::swap(&mut a_vec, &mut b_vec);
        std::mem::swap(&mut a_norm, &mut b_norm);
    }

    loop {
        let n = dot(&a_vec, &b_vec);
        let r = nearest_integer(&n, &b_norm);
        let two_r_n = BigInt::from(2) * &r * &n;
        let r_sq_b = &r * &r * &b_norm;
        let t_scalar = &a_norm - two_r_n + r_sq_b;

        if t_scalar >= b_norm {
            return (a_vec, b_vec);
        }

        let t_vec = sub_scaled(&a_vec, &b_vec, &r);

        a_vec = b_vec;
        b_vec = t_vec;
        a_norm = b_norm;
        b_norm = t_scalar;
    }
}

fn dot(x: &[BigInt], y: &[BigInt]) -> BigInt {
    x.iter()
        .zip(y.iter())
        .fold(BigInt::zero(), |acc, (xi, yi)| acc + xi * yi)
}

fn sub_scaled(a: &[BigInt], b: &[BigInt], r: &BigInt) -> Vec<BigInt> {
    a.iter().zip(b.iter()).map(|(ai, bi)| ai - r * bi).collect()
}

fn nearest_integer(n: &BigInt, b: &BigInt) -> BigInt {
    let two = BigInt::from(2);
    let num = if n.is_negative() {
        two.clone() * n - b
    } else {
        two.clone() * n + b
    };
    num.div_floor(&(two * b))
}

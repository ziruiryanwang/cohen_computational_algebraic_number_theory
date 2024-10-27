use num_integer::Integer;
use num_traits::Signed;

use crate::group::GroupElement;

pub fn left_right_binary_bits_power<G, I>(g: G, n: I) -> G
where
    G: GroupElement,
    I: Integer + Signed + Clone,
{
    if n.is_zero() {
        return G::identity();
    }

    let n_abs = n.abs();
    let z = if n.is_negative() { g.inverse() } else { g };
    let mut y = z.clone();
    let mut f = msb_index(&n_abs);

    loop {
        if f == 0 {
            return y;
        }
        f -= 1;
        y = y.mul(&y);
        if bit_at(&n_abs, f) {
            y = y.mul(&z);
        }
    }
}

fn msb_index<I>(n: &I) -> u64
where
    I: Integer + Signed + Clone,
{
    let two = num_traits::one::<I>() + num_traits::one::<I>();
    let mut temp = n.clone();
    let mut e = 0u64;
    while temp >= two {
        temp = temp.div_floor(&two);
        e += 1;
    }
    e
}

fn bit_at<I>(n: &I, f: u64) -> bool
where
    I: Integer + Signed + Clone,
{
    let two = num_traits::one::<I>() + num_traits::one::<I>();
    let mut q = n.clone();
    for _ in 0..f {
        q = q.div_floor(&two);
    }
    q.is_odd()
}

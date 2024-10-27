use num_integer::Integer;
use num_traits::Signed;

use crate::group::GroupElement;

pub fn left_right_binary_power<G, I>(g: G, n: I) -> G
where
    G: GroupElement,
    I: Integer + Signed + Clone,
{
    if n.is_zero() {
        return G::identity();
    }

    let mut n_abs = n.abs();
    let z = if n.is_negative() { g.inverse() } else { g };
    let mut y = z.clone();
    let one: I = num_traits::one();
    let two: I = one.clone() + one.clone();

    let mut e = one.clone();
    loop {
        let next = e.clone() + e.clone();
        if next <= n_abs {
            e = next;
        } else {
            break;
        }
    }

    n_abs = n_abs - e.clone();

    loop {
        if e == one {
            return y;
        }
        e = e.div_floor(&two);
        y = y.mul(&y);
        if n_abs >= e {
            n_abs = n_abs - e.clone();
            y = y.mul(&z);
        }
    }
}

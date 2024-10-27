use num_integer::Integer;
use num_traits::Signed;

use crate::group::GroupElement;

pub fn right_left_binary_power<G, I>(g: G, n: I) -> G
where
    G: GroupElement,
    I: Integer + Signed + Clone,
{
    if n.is_zero() {
        return G::identity();
    }

    let mut y = G::identity();
    let mut n_abs = n.abs();
    let mut z = if n.is_negative() { g.inverse() } else { g };
    let two: I = num_traits::one::<I>() + num_traits::one::<I>();

    loop {
        if n_abs.is_odd() {
            y = z.mul(&y);
        }
        n_abs = n_abs.div_floor(&two);
        if n_abs.is_zero() {
            return y;
        }
        z = z.mul(&z);
    }
}

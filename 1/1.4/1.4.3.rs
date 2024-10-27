use num_bigint::BigInt;

use crate::{right_left_binary_power, GroupElement};

pub fn order_of_element<G>(g: G, h: &BigInt, factors: &[(BigInt, u32)]) -> BigInt
where
    G: GroupElement,
{
    let mut e = h.clone();
    for (p, v) in factors {
        let p_power: BigInt = p.pow(*v);
        e /= &p_power;
        let mut g1 = right_left_binary_power(g.clone(), e.clone());
        while g1 != G::identity() {
            g1 = right_left_binary_power(g1, p.clone());
            e *= p;
        }
    }
    e
}

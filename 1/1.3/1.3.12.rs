use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::One;

use crate::extended_euclid;

pub fn chinese_remainder_inductive(mi: &[BigInt], xi: &[BigInt]) -> Option<BigInt> {
    if mi.len() != xi.len() || mi.is_empty() {
        return None;
    }

    let mut m = mi[0].clone();
    let mut x = xi[0].mod_floor(&m);

    for i in 1..mi.len() {
        let m_i = mi[i].clone();
        let (u, v, d) = extended_euclid(m.clone(), m_i.clone());
        if !d.is_one() {
            return None;
        }
        let term = &u * m.clone() * xi[i].mod_floor(&m_i) + v * m_i.clone() * x;
        m *= &m_i;
        x = term.mod_floor(&m);
    }

    Some(x.mod_floor(&m))
}

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::One;

use crate::extended_euclid;

pub fn chinese_remainder_pairwise(mi: &[BigInt], xi: &[BigInt]) -> Option<BigInt> {
    if mi.len() != xi.len() || mi.is_empty() {
        return None;
    }

    let mut c: Vec<BigInt> = Vec::with_capacity(mi.len());
    let mut p = mi[0].clone();
    c.push(BigInt::one());

    for j in 1..mi.len() {
        let (u, _v, d) = extended_euclid(p.clone(), mi[j].clone());
        if !d.is_one() {
            return None;
        }
        c.push(u);
        p *= &mi[j];
    }

    let mut y = xi[0].mod_floor(&mi[0]);
    let mut m_prod = mi[0].clone();

    for j in 1..mi.len() {
        let m_j = &mi[j];
        let c_j = &c[j];
        let delta = (xi[j].mod_floor(m_j) - y.mod_floor(m_j)) * c_j;
        let yj = delta.mod_floor(m_j);
        y += &m_prod * yj;
        m_prod *= m_j;
    }

    Some(y.mod_floor(&m_prod))
}

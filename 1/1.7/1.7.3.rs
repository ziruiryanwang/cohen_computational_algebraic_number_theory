use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, ToPrimitive, Zero};

use crate::arith::int_sqrt;

fn residue_table(modulus: u32) -> Vec<bool> {
    let m = modulus as usize;
    let mut table = vec![false; m];
    for x in 0..modulus {
        let r = ((x as u64 * x as u64) % modulus as u64) as usize;
        table[r] = true;
    }
    table
}

pub fn square_test(n: BigInt) -> Option<BigInt> {
    if n.is_zero() {
        return Some(BigInt::zero());
    }
    if n.is_negative() {
        return None;
    }

    let q64 = residue_table(64);
    let t = n.mod_floor(&BigInt::from(64u32)).to_usize()?;
    if !q64[t] {
        return None;
    }
    let r = n.mod_floor(&BigInt::from(45045u32)).to_usize()?;

    let q63 = residue_table(63);
    if !q63[r % 63] {
        return None;
    }
    let q65 = residue_table(65);
    if !q65[r % 65] {
        return None;
    }
    let q11 = residue_table(11);
    if !q11[r % 11] {
        return None;
    }

    let q = int_sqrt(&n);
    if &q * &q == n {
        Some(q)
    } else {
        None
    }
}

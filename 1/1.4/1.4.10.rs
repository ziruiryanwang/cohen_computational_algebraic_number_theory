use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

pub fn kronecker(a: BigInt, b: BigInt) -> i32 {
    let mut a = a;
    let mut b = b;
    let mut k: i32 = 1;

    if b.is_zero() {
        return if a.abs() == BigInt::one() { 1 } else { 0 };
    }

    if a.is_zero() {
        return if b.abs() == BigInt::one() { 1 } else { 0 };
    }

    loop {
        // remove powers of 2 from b
        let mut v = 0usize;
        while b.is_even() {
            v += 1;
            b >>= 1;
        }
        if v % 2 == 1 {
            let a_mod8 = (a.mod_floor(&BigInt::from(8))).to_i32().unwrap();
            if a_mod8 == 3 || a_mod8 == 5 {
                k = -k;
            }
        }
        if b.is_negative() {
            b = -b;
            if a.is_negative() {
                k = -k;
            }
        }

        if b == BigInt::one() {
            return k;
        }
        if a.is_zero() {
            return 0;
        }

        // remove powers of 2 from a
        let mut v = 0usize;
        while a.is_even() {
            v += 1;
            a >>= 1;
        }
        if v % 2 == 1 {
            let b_mod8 = (b.mod_floor(&BigInt::from(8))).to_i32().unwrap();
            if b_mod8 == 3 || b_mod8 == 5 {
                k = -k;
            }
        }

        // reciprocity
        if (a.mod_floor(&BigInt::from(4))).to_i32().unwrap() == 3
            && (b.mod_floor(&BigInt::from(4))).to_i32().unwrap() == 3
        {
            k = -k;
        }

        let r = a.abs();
        let new_a = b.mod_floor(&r);
        b = r;
        a = new_a;
    }
}

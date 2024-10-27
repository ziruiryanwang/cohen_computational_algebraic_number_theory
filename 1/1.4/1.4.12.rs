use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

pub fn kronecker_binary(a: BigInt, b: BigInt) -> i32 {
    let mut a = a;
    let mut b = b;

    if b.is_zero() {
        return if a.abs() == BigInt::one() { 1 } else { 0 };
    }

    if a.is_even() && b.is_even() {
        return 0;
    }

    let mut k: i32 = 1;

    loop {
        // Step 2: remove factors of 2 from b
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

        // Step 3: reduce size once (b is now odd and positive)
        a = a.mod_floor(&b);

        // Step 4: finished?
        if a.is_zero() {
            return if b == BigInt::one() { k } else { 0 };
        }

        // Step 5: remove powers of 2 from a
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

        // Step 6: subtract and apply reciprocity (a and b odd)
        let r = &b - &a;
        if r.is_positive() {
            let a_mod4 = (a.mod_floor(&BigInt::from(4))).to_i32().unwrap();
            let b_mod4 = (b.mod_floor(&BigInt::from(4))).to_i32().unwrap();
            if a_mod4 == 3 && b_mod4 == 3 {
                k = -k;
            }
            b = a;
            a = r;
        } else {
            a = -r;
        }
    }
}

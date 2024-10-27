use num_bigint::BigInt;

use crate::arith::int_sqrt;

pub fn integer_sqrt(n: BigInt) -> BigInt {
    int_sqrt(&n)
}

use num_integer::Integer;
use num_traits::Signed;

pub fn euclid_gcd<I>(a: I, b: I) -> I
where
    I: Integer + Signed + Clone,
{
    let mut a = a.abs();
    let mut b = b.abs();

    while !b.is_zero() {
        let r = a.mod_floor(&b);
        a = b;
        b = r;
    }

    a
}

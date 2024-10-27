use num_integer::Integer;
use num_traits::{Signed, ToPrimitive};

use crate::group::GroupElement;

pub fn left_right_base_2k_power<G, I>(g: G, n: I, k: usize) -> G
where
    G: GroupElement,
    I: Integer + Signed + Clone + ToPrimitive,
{
    assert!(k > 0);

    if n.is_zero() {
        return G::identity();
    }

    let n_abs = n.abs();
    let z = if n.is_negative() { g.inverse() } else { g };
    let mut base = num_traits::one::<I>();
    for _ in 0..k {
        base = base.clone() + base.clone();
    }

    let e = msd_index_base(&n_abs, &base);
    let mut f = e;
    let table = precompute_odd_powers(&z, k);
    let mut y = G::identity();

    loop {
        let alpha = digit(&n_abs, &base, f);
        if alpha == 0 {
            for _ in 0..k {
                y = y.mul(&y);
            }
        } else {
            let mut alpha_val = alpha;
            let mut t = 0usize;
            while alpha_val % 2 == 0 {
                alpha_val /= 2;
                t += 1;
            }
            let b = alpha_val;
            if f != e {
                for _ in 0..(k - t) {
                    y = y.mul(&y);
                }
                y = y.mul(table[b].as_ref().unwrap());
            } else {
                y = table[b].as_ref().unwrap().clone();
            }
            for _ in 0..t {
                y = y.mul(&y);
            }
        }

        if f == 0 {
            return y;
        }
        f -= 1;
    }
}

fn msd_index_base<I>(n: &I, base: &I) -> u64
where
    I: Integer + Signed + Clone,
{
    let mut temp = n.clone();
    let mut e = 0u64;
    while temp >= base.clone() {
        temp = temp.div_floor(base);
        e += 1;
    }
    e
}

fn digit<I>(n: &I, base: &I, f: u64) -> usize
where
    I: Integer + Signed + Clone + ToPrimitive,
{
    let mut q = n.clone();
    for _ in 0..f {
        q = q.div_floor(base);
    }
    let rem = q.mod_floor(base);
    rem.to_usize().unwrap()
}

fn precompute_odd_powers<G>(z: &G, k: usize) -> Vec<Option<G>>
where
    G: GroupElement,
{
    let size = 1usize << k;
    let mut table: Vec<Option<G>> = vec![None; size];
    let z2 = z.mul(z);
    let mut current = z.clone();

    for b in (1..size).step_by(2) {
        if b == 1 {
            table[b] = Some(current.clone());
        } else {
            current = current.mul(&z2);
            table[b] = Some(current.clone());
        }
    }

    table
}

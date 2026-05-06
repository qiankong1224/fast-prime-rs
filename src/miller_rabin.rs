use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, Zero, ToPrimitive};

pub fn miller_rabin(n: &BigUint, base: &BigUint) -> bool {
    if n <= &BigUint::one() { return false; }
    if n.to_u32() == Some(2) { return true; }
    if n.is_even() { return false; }

    let n_minus_1 = n - 1u32;
    // trailing_zeros 返回 Option<u64>
    let s = n_minus_1.trailing_zeros().unwrap_or(0);
    let d = &n_minus_1 >> s;

    let mut b = base % n;
    if b.is_zero() { return true; }

    check_mr(n, &n_minus_1, s, &d, &b)
}

pub fn miller_rabin_multi_round(n: &BigUint, rounds: usize) -> bool {
    if n <= &BigUint::one() { return false; }
    let n_u32 = n.to_u32();
    if n_u32 == Some(2) || n_u32 == Some(3) { return true; }
    if n.is_even() { return false; }

    let n_minus_1 = n - 1u32;
    let s = n_minus_1.trailing_zeros().unwrap_or(0);
    let d = &n_minus_1 >> s;

    const BASES: [u32; 24] = [
        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
        43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97
    ];

    for i in 0..rounds {
        let base_biguint = if i < BASES.len() {
            let b = BASES[i];
            if n_u32 == Some(b) { return true; }
            BigUint::from(b)
        } else {
            BigUint::from(97 + (i - BASES.len() + 1) as u32 * 2)
        };

        let b = &base_biguint % n;
        if b.is_zero() { continue; }

        if !check_mr(n, &n_minus_1, s, &d, &b) {
            return false;
        }
    }
    true
}


#[inline]
fn check_mr(n: &BigUint, n_minus_1: &BigUint, s: u64, d: &BigUint, base: &BigUint) -> bool {
    let mut x = base.modpow(d, n);

    if x.is_one() || x == *n_minus_1 {
        return true;
    }


    for _ in 0..s.saturating_sub(1) {
        x = (&x * &x) % n;

        if x == *n_minus_1 {
            return true;
        }
        if x.is_one() {
            return false;
        }
    }
    false
}
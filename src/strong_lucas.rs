use num_bigint::{BigUint};
use num_integer::Integer;
use num_traits::{One, Zero, ToPrimitive};

pub fn extra_strong_lucas(n: &BigUint) -> bool {
    // 参数选择：寻找第一个 Jacobi(D, n) = -1
    let mut p = 3u32;
    let d = loop {
        let d_val = p * p - 4;
        let j = jacobi_u(d_val, n);
        if j == -1 { break p; }
        if j == 0 {
            let big_p = BigUint::from(p);
            return n == &(&big_p + 2u32);
        }
        p += 1;
        if p == 43 {
            let s = n.sqrt();
            if &s * &s == *n { return false; }
        }
    };

    let n_plus_1 = n + 1u32;
    let r = n_plus_1.trailing_zeros().unwrap_or(0);
    let s_val = &n_plus_1 >> r;

    let big_p = BigUint::from(d);
    let mut vk = BigUint::from(2u32);
    let mut vk1 = big_p.clone();
    let two = BigUint::from(2u32);

    // 核心梯子循环：全部使用 BigUint 减法优化
    for i in (0..s_val.bits()).rev() {
        if s_val.bit(i) {
            // vk = (vk * vk1 - big_p) % n
            vk = sub_mod(&(&vk * &vk1), &big_p, n);
            // vk1 = (vk1 * vk1 - 2) % n
            vk1 = sub_mod(&(&vk1 * &vk1), &two, n);
        } else {
            // vk1 = (vk * vk1 - big_p) % n
            vk1 = sub_mod(&(&vk * &vk1), &big_p, n);
            // vk = (vk * vk - 2) % n
            vk = sub_mod(&(&vk * &vk), &two, n);
        }
    }

    if vk == two || vk == (n - 2u32) {
        // P*vk - 2*vk1 ≡ 0 (mod n)
        let left = (&big_p * &vk) % n;
        let right = (&vk1 * 2u32) % n;
        if left == right { return true; }
    }

    for _ in 0..r - 1 {
        if vk.is_zero() { return true; }
        if vk == two { return false; }
        vk = sub_mod(&(&vk * &vk), &two, n);
    }
    false
}

// 专门处理 (a - b) % n，避开 BigInt 转换
#[inline]
fn sub_mod(a: &BigUint, b: &BigUint, n: &BigUint) -> BigUint {
    let a_m = a % n;
    let b_m = b % n;
    if a_m >= b_m {
        a_m - b_m
    } else {
        a_m + n - b_m
    }
}

fn jacobi_u(a: u32, n: &BigUint) -> i32 {
    let mut a = BigUint::from(a) % n;
    let mut n = n.clone();
    let mut t = 1;
    while !a.is_zero() {
        while a.is_even() {
            a >>= 1;
            let n_mod_8 = (n.clone() % 8u32).to_u32().unwrap();
            if n_mod_8 == 3 || n_mod_8 == 5 { t = -t; }
        }
        std::mem::swap(&mut a, &mut n);
        if (a.clone() % 4u32).to_u32().unwrap() == 3 && (n.clone() % 4u32).to_u32().unwrap() == 3 {
            t = -t;
        }
        a %= &n;
    }
    if n.is_one() { t } else { 0 }
}
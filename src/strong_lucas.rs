use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, Zero, ToPrimitive};

pub fn extra_strong_lucas(n: &BigUint) -> bool {
    // 1. 参数选择优化
    let mut p = 3u32;
    let d_val = loop {
        let d = p * p - 4;
        let j = jacobi_fast(d, n);
        if j == -1 { break p; }
        if j == 0 {
            // 如果 j == 0 且 n > d+2，则 n 必为合数.n 是 d 的因子
            // 注意n == p+2 的情况处理
            let big_p_plus_2 = BigUint::from(p + 2);
            return n == &big_p_plus_2;
        }
        p += 1;
        // 安全限制
        if p == 1000 {
            let s = n.sqrt();
            if &s * &s == *n { return false; }
        }
    };

    let n_plus_1 = n + 1u32;
    let r = n_plus_1.trailing_zeros().unwrap_or(0);
    let s_val = &n_plus_1 >> r;

    let big_p = BigUint::from(d_val);
    let mut vk = BigUint::from(2u32);
    let mut vk1 = big_p.clone();
    let two = BigUint::from(2u32);


    // 提前获取位长度，避免在循环中重复计算
    let bits = s_val.bits();
    for i in (0..bits).rev() {
        if s_val.bit(i) {
            // vk = (vk * vk1 - big_p) % n
            vk = mul_sub_mod(&vk, &vk1, &big_p, n);
            // vk1 = (vk1 * vk1 - 2) % n
            vk1 = mul_sub_mod(&vk1, &vk1, &two, n);
        } else {
            // vk1 = (vk * vk1 - big_p) % n
            vk1 = mul_sub_mod(&vk, &vk1, &big_p, n);
            // vk = (vk * vk - 2) % n
            vk = mul_sub_mod(&vk, &vk, &two, n);
        }
    }


    if vk == two || vk == (n - 2u32) {
        let left = (&big_p * &vk) % n;
        let right = (&vk1 * 2u32) % n;
        if left == right { return true; }
    }

    for _ in 0..r - 1 {
        if vk.is_zero() { return true; }
        if vk == two { return false; }
        vk = mul_sub_mod(&vk, &vk, &two, n);
    }
    false
}


#[inline]
fn mul_sub_mod(a: &BigUint, b: &BigUint, c: &BigUint, n: &BigUint) -> BigUint {
    let prod = (a * b) % n;
    if prod >= *c {
        prod - c
    } else {
        prod + n - c
    }
}


fn jacobi_fast(a_u32: u32, n: &BigUint) -> i32 {
    let mut a = (a_u32 % n).to_u32().unwrap_or_else(|| {
        // 如果 n 小于 u32，直接取模转换
        (n % a_u32).to_u32().unwrap()
    });

    // 获取 n % 8 用于计算性质
    let mut n_low = n.to_u32_digits().first().cloned().unwrap_or(0);
    let mut n_big = n.clone();
    let mut t = 1;

    // 处理 BigUint
    if a == 0 { return 0; }

    while a != 0 {
        // 消除 a 中的 2
        let tz = a.trailing_zeros();
        a >>= tz;
        if tz % 2 == 1 {
            let mod8 = n_low & 7;
            if mod8 == 3 || mod8 == 5 {
                t = -t;
            }
        }

        // 互反律Jacobi(a, n)

        let old_a = a;

        a = (n_big % old_a).to_u32().unwrap();

        // 更新 n_low 为 old_a
        if (old_a & 3) == 3 && (n_low & 3) == 3 {
            t = -t;
        }
        n_low = old_a;
        // 之后所有的 n_big 运算都在 u32 范围内
        n_big = BigUint::from(old_a);
    }

    if n_low == 1 { t } else { 0 }
}


fn _jacobi_u32(mut a: u32, mut n: u32) -> i32 {
    let mut t = 1;
    a %= n;
    while a != 0 {
        while a % 2 == 0 {
            a /= 2;
            if n % 8 == 3 || n % 8 == 5 { t = -t; }
        }
        std::mem::swap(&mut a, &mut n);
        if a % 4 == 3 && n % 4 == 3 { t = -t; }
        a %= n;
    }
    if n == 1 { t } else { 0 }
}
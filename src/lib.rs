pub mod miller_rabin;
pub mod strong_lucas;

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{ToPrimitive, Zero};
use crate::miller_rabin::{miller_rabin, miller_rabin_multi_round};
use crate::strong_lucas::extra_strong_lucas;

// 预计算的素数乘积，用于快速预筛选
const P_3_53: u64 = 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47 * 53;
const P_59_97: u64 = 59 * 61 * 67 * 71 * 73 * 79 * 83 * 89 * 97;

/// 内部核心逻辑：判断是否通过 Baillie-PSW 测试
fn bpsw_internal(n: &BigUint) -> bool {
    // 1. MR Base 2 (强伪素数测试)
    // 这里的 miller_rabin 内部应尽量使用静态 BigUint::from(2)
    if !miller_rabin(n, &BigUint::from(2u32)) {
        return false;
    }

    // 2. Extra Strong Lucas 测试
    // 这是 Baillie-PSW 的第二部分，能排除掉绝大多数 MR Base 2 的伪素数
    if !extra_strong_lucas(n) {
        return false;
    }

    true
}

pub fn is_prime(n: &BigUint) -> bool {
    // 小规模数字处理
    let bits = n.bits();
    if bits <= 64 {
        let w = n.to_u64().unwrap_or(0);
        return is_prime_u64(w);
    }


    if n.is_even() { return false; }


    let r1 = (n % P_3_53).to_u64().unwrap();
    if r1 % 3 == 0 || r1 % 5 == 0 || r1 % 7 == 0 || r1 % 11 == 0 || r1 % 13 == 0 ||
        r1 % 17 == 0 || r1 % 19 == 0 || r1 % 23 == 0 || r1 % 29 == 0 || r1 % 31 == 0 ||
        r1 % 37 == 0 || r1 % 41 == 0 || r1 % 43 == 0 || r1 % 47 == 0 || r1 % 53 == 0 {
        return false;
    }

    let r2 = (n % P_59_97).to_u64().unwrap();
    if r2 % 59 == 0 || r2 % 61 == 0 || r2 % 67 == 0 || r2 % 71 == 0 || r2 % 73 == 0 ||
        r2 % 79 == 0 || r2 % 83 == 0 || r2 % 89 == 0 || r2 % 97 == 0 {
        return false;
    }


    bpsw_internal(n)
}

pub fn is_prime_mr25(n: &BigUint) -> bool {
    // 先走标准的 is_prime (BPSW)
    if !is_prime(n) {
        return false;
    }

    // 对于极高安全性需求，追加 24 轮随机底数的 MR 测试
    if !miller_rabin_multi_round(n, 24) {
        return false;
    }

    true
}

#[inline]
fn is_prime_u64(n: u64) -> bool {
    if n < 64 {
        const PRIME_BIT_MASK: u64 = 0x28208a20a08a28ac;
        return (PRIME_BIT_MASK & (1 << n)) != 0;
    }
    if n % 2 == 0 || n % 3 == 0 { return false; }

    // 使用蒙哥马利加速的 Miller-Rabin
    let mont = Montgomery::new(n);
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;

    // u64 确定性底数
    const BASES: [u64; 7] = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];
    for &base in &BASES {
        if n <= base { break; }
        if !miller_rabin_montgomery(n, base, s, d, &mont) {
            return false;
        }
    }
    true
}

struct Montgomery {
    n: u64,
    n_inv: u64,
    r2: u64,
}

impl Montgomery {
    #[inline]
    fn new(n: u64) -> Self {
        let mut inv = n;
        for _ in 0..5 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
        }
        let n_inv = inv.wrapping_neg();
        let r = (1u128.wrapping_shl(64) % n as u128) as u64;
        let r2 = ((r as u128 * r as u128) % n as u128) as u64;
        Montgomery { n, n_inv, r2 }
    }

    #[inline]
    fn redc(&self, t: u128) -> u64 {
        let m = (t as u64).wrapping_mul(self.n_inv);
        let mn = m as u128 * self.n as u128;

        // 优化高位进位处理
        let (res_low, carry) = (t as u64).overflowing_add(mn as u64);
        let mut res_hi = (t >> 64) + (mn >> 64) + (if carry { 1 } else { 0 });

        if res_hi >= self.n as u128 {
            res_hi -= self.n as u128;
        }
        res_hi as u64
    }

    #[inline]
    fn transform(&self, x: u64) -> u64 {
        self.redc(x as u128 * self.r2 as u128)
    }

    #[inline]
    fn mul(&self, a_prime: u64, b_prime: u64) -> u64 {
        self.redc(a_prime as u128 * b_prime as u128)
    }

    #[inline]
    fn pow(&self, mut a_prime: u64, mut exp: u64) -> u64 {
        let mut res_prime = self.transform(1);
        while exp > 0 {
            if exp & 1 == 1 { res_prime = self.mul(res_prime, a_prime); }
            a_prime = self.mul(a_prime, a_prime);
            exp >>= 1;
        }
        res_prime
    }
}

#[inline]
fn miller_rabin_montgomery(n: u64, a: u64, s: u32, d: u64, mont: &Montgomery) -> bool {
    let a_prime = mont.transform(a);
    let mut x_prime = mont.pow(a_prime, d);
    let one_prime = mont.transform(1);
    let n_minus_1_prime = mont.transform(n - 1);

    if x_prime == one_prime || x_prime == n_minus_1_prime {
        return true;
    }
    for _ in 0..s - 1 {
        x_prime = mont.mul(x_prime, x_prime);
        if x_prime == n_minus_1_prime {
            return true;
        }
    }
    false
}
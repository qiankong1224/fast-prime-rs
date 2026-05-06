pub mod miller_rabin;
pub mod strong_lucas;

use num_bigint::{BigUint};
use num_integer::Integer;
use num_traits::{One, ToPrimitive};
use crate::miller_rabin::miller_rabin;
use crate::strong_lucas::extra_strong_lucas;
use crate::miller_rabin::miller_rabin_multi_round;

#[inline]
fn is_prime_u64(n: u64) -> bool {
    // 处理 n < 64 的情况
    if n < 64 {
        const PRIME_BIT_MASK: u64 = 0x28208a20a08a28ac;
        return (PRIME_BIT_MASK & (1 << n)) != 0;
    }

    //排除偶数和 3 的倍数
    if n % 2 == 0 || n % 3 == 0 { return false; }

    //原生积取模筛法 (已经通过 2, 3，现在筛 5-97)
    const P_A: u64 = 5 * 7 * 11 * 13 * 17 * 19 * 23 * 37;
    const P_B: u64 = 29 * 31 * 41 * 43 * 47 * 53;
    const P_C: u64 = 59 * 61 * 67 * 71 * 73 * 79 * 83 * 89 * 97;

    if n % P_A == 0 && n <= 37 { return true; } // 排除本身
    if n % P_A % 5 == 0 || n % P_A % 7 == 0 || n % P_A % 11 == 0 || n % P_A % 13 == 0 ||
        n % P_A % 17 == 0 || n % P_A % 19 == 0 || n % P_A % 23 == 0 || n % P_A % 37 == 0 { return false; }

    if n % P_B == 0 && n <= 53 { return true; }
    if n % P_B % 29 == 0 || n % P_B % 31 == 0 || n % P_B % 41 == 0 || n % P_B % 43 == 0 ||
        n % P_B % 47 == 0 || n % P_B % 53 == 0 { return false; }

    if n % P_C == 0 && n <= 97 { return true; }
    if n % P_C % 59 == 0 || n % P_C % 61 == 0 || n % P_C % 67 == 0 || n % P_C % 71 == 0 ||
        n % P_C % 73 == 0 || n % P_C % 79 == 0 || n % P_C % 83 == 0 || n % P_C % 89 == 0 || n % P_C % 97 == 0 { return false; }

    //蒙哥马利 Miller-Rabin
    let mont = Montgomery::new(n);
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;

    // 对于 u64，这 7 个底数是确定的
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
    n_inv: u64, // n * n_inv ≡ -1 (mod 2^64)
    r2: u64,    // R^2 mod n
}

impl Montgomery {
    #[inline]
    fn new(n: u64) -> Self {
        // 牛顿迭代法计算模逆 n_inv
        let mut inv = n;
        for _ in 0..5 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
        }
        let n_inv = inv.wrapping_neg();

        // 计算 R^2 mod n，其中 R = 2^64
        let r = (1u128.wrapping_shl(64) % n as u128) as u64;
        let r2 = ((r as u128 * r as u128) % n as u128) as u64;

        Montgomery { n, n_inv, r2 }
    }

    #[inline]
    fn redc(&self, t: u128) -> u64 {
        let m = (t as u64).wrapping_mul(self.n_inv);
        let mn = m as u128 * self.n as u128;


        let (_lo, carry) = (t as u64).overflowing_add(mn as u64);
        let hi = (t >> 64) + (mn >> 64) + (if carry { 1 } else { 0 });


        let res = if hi >= self.n as u128 {
            (hi - self.n as u128) as u64
        } else {
            hi as u64
        };
        res
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
    // 整个过程都在蒙哥马利空间完成，避开除法指令
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
pub fn is_prime(n: &BigUint) -> bool {
    let bits = n.bits();

    if bits <= 64 {
        let w = n.to_u64().unwrap_or(0);
        return is_prime_u64(w);
    }


    if n.is_even() { return false; }

    // 增强预筛选筛掉 100 以内所有素数的倍数
    const P3_53: u64 = 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47 * 53;
    if n.gcd(&BigUint::from(P3_53)) > BigUint::one() {
        return false;
    }

    // 继续筛 59-100 之间的素数
    const P59_97: u64 = 59 * 61 * 67 * 71 * 73 * 79 * 83 * 89 * 97;
    if n.gcd(&BigUint::from(P59_97)) > BigUint::one() {
        return false;
    }

    // Miller-Rabin
    if !miller_rabin(n, &BigUint::from(2u32)) {
        return false;
    }

    // Extra Strong Lucas
    if !extra_strong_lucas(n) {
        return false;
    }

    true
}



pub fn is_prime_mr25(n: &BigUint) -> bool {
    let bits = n.bits();

    if bits <= 64 {
        let w = n.to_u64().unwrap_or(0);
        return is_prime_u64(w);
    }


    if n.is_even() { return false; }

    const P3_53: u64 = 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47 * 53;
    if n.gcd(&BigUint::from(P3_53)) > BigUint::one() {
        return false;
    }

    const P59_97: u64 = 59 * 61 * 67 * 71 * 73 * 79 * 83 * 89 * 97;
    if n.gcd(&BigUint::from(P59_97)) > BigUint::one() {
        return false;
    }


    if !miller_rabin(n, &BigUint::from(2u32)) {
        return false;
    }


    if !extra_strong_lucas(n) {
        return false;
    }

    //上面一样，结束后再进行24次mille-rabin
    if !miller_rabin_multi_round(n, 24) {
        return false;
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;

    // 清洗字符串逻辑
    fn parse_big(s: &str) -> BigUint {
        let cleaned = s.replace(|c: char| c.is_whitespace(), "");
        BigUint::parse_bytes(cleaned.as_bytes(), 10).expect("Invalid test data")
    }

    const PRIMES: &[&str] = &[
        "2", "3", "5", "7", "11", "13", "17", "19", "23",
        "13756265695458089029",
        "13496181268022124907",
        // ECC: Curve25519 (2^255 - 19)
        "57896044618658097711785492504343953926634992332820282019728792003956564819949",
        // E-521 (2^521 - 1)
        "6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151",
    ];

    const COMPOSITES: &[&str] = &[
        "0", "1", "4", "6", "8", "9", "10", "15", "21", "25",
        "561",   // Carmichael number
        "2047",  // MR-Base2 pseudoprime: 23 * 89
        "3277",  // MR-Base2 pseudoprime: 29 * 113
        "8321",  // MR-Base2 pseudoprime: 53 * 157
        "1194649", // Strong Lucas pseudoprime
        // Arnault's tricky composite (passes MR 2, 3, 5, 7... up to 29)
        "1195068768795265792518361315725116351898245581",
        // OEIS A217719: Extra-strong Lucas pseudoprimes
        "989", "3239", "5777", "10877", "27971", "30739",
    ];

    #[test]
    fn test_is_prime_on_primes() {
        for s in PRIMES {
            let n = parse_big(s);
            assert!(is_prime(&n), "Should be prime: {}", s);
        }
    }

    #[test]
    fn test_is_prime_on_composites() {
        for s in COMPOSITES {
            let n = parse_big(s);
            assert!(!is_prime(&n), "Should be composite: {}", s);
        }
    }

    /// 特殊测试：验证 BPSW 的必要性
    /// 验证那些能骗过 Miller-Rabin Base 2 但骗不过 Lucas 的数字
    #[test]
    fn test_miller_rabin_pseudoprimes() {
        // 这些数字在 Base 2 下会返回 true (错误)，但 BPSW 整体应该返回 false (正确)
        let mr_pseudoprimes = ["2047", "3277", "4033", "4681", "8321"];
        for s in &mr_pseudoprimes {
            let n = parse_big(s);
            // 单独测 MR 应该失败 (返回 true 表示它以为是素数)
            assert!(miller_rabin(&n, &BigUint::from(2u32)), "MR base 2 should be fooled by {}", s);
            // 但 BPSW 整体应该识破它
            assert!(!is_prime(&n), "BPSW should have caught MR pseudoprime {}", s);
        }
    }

    /// 特殊测试：验证 BPSW 的必要性
    /// 验证那些能骗过 Lucas 测试但骗不过 Miller-Rabin Base 2 的数字
    #[test]
    fn test_lucas_pseudoprimes() {
        let lucas_pseudoprimes = ["989", "3239", "5777", "10877"];
        for s in &lucas_pseudoprimes {
            let n = parse_big(s);
            // 单独测 Lucas 应该失败 (返回 true)
            assert!(extra_strong_lucas(&n), "Lucas should be fooled by {}", s);
            // 但 BPSW 整体应该识破它
            assert!(!is_prime(&n), "BPSW should have caught Lucas pseudoprime {}", s);
        }
    }
}
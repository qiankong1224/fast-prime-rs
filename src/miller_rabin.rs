use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, Zero};

pub fn miller_rabin(n: &BigUint, base: &BigUint) -> bool {
    if n <= &BigUint::one() { return false; }
    if n == &BigUint::from(2u32) { return true; }
    if n.is_even() { return false; }
    if base % n == BigUint::zero() { return true; }

    let n_minus_1 = n - 1u32;

    // 快速分解 n-1 = 2^s * d
    // trailing_zeros() 在 BigUint 上非常快，它直接数二进制末尾的 0
    let s = n_minus_1.trailing_zeros().unwrap_or(0);
    let d = &n_minus_1 >> s;

    // x = base^d mod n
    let mut x = base.modpow(&d, n);

    // 如果 base^d ≡ 1 (mod n) 或 base^d ≡ n-1 (mod n)，则可能是素数
    if x.is_one() || x == n_minus_1 {
        return true;
    }

    // 进行 s-1 次平方
    for _ in 0..s - 1 {
        // x = x^2 mod n
        // 这里的 &x * &x 是关键，避免了 x 的所有权移动
        x = (&x * &x) % n;

        if x == n_minus_1 {
            return true;
        }
        if x.is_one() {
            return false;
        }
    }
    false
}

pub fn miller_rabin_multi_round(n: &BigUint, rounds: usize) -> bool {
    // 1基础边界检查
    if n <= &BigUint::one() { return false; }
    if n == &BigUint::from(2u32) || n == &BigUint::from(3u32) { return true; }
    if n.is_even() { return false; }

    // 分解 n-1 = 2^s * d
    let n_minus_1 = n - 1u32;
    let s = n_minus_1.trailing_zeros().unwrap_or(0);
    let d = &n_minus_1 >> s;


    const BASES: [u32; 24] = [
        3, 5, 7, 11, 13, 17, 19, 23, 29,
        31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
        73, 79, 83, 89, 97
    ];

    // 开始多轮迭代
    for i in 0..rounds {
        let base_val = if i < 25 {
            BigUint::from(BASES[i])
        } else {

            BigUint::from(BASES[23] + (i - 24) as u32)
        };

        // 如果 n 本身就是底数中的一个，返回 true
        if n == &base_val { return true; }
        // 如果底数大于 n，要做取模
        let base = &base_val % n;
        if base.is_zero() { continue; }

        // 计算 x = base^d mod n
        let mut x = base.modpow(&d, n);

        // 如果 x == 1 或 x == n-1，通过这一轮
        if x.is_one() || x == n_minus_1 {
            continue;
        }

        // 进行 s-1 次平方
        let mut composite = true;
        for _ in 0..s - 1 {
            // x = x^2 mod n
            x = (&x * &x) % n;

            if x == n_minus_1 {
                composite = false;
                break;
            }
            if x.is_one() {
                return false; // 绝对是合数
            }
        }

        if composite {
            return false;
        }
    }

    true
}
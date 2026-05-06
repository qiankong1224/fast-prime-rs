use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use num_bigint::BigUint;
use rug::Integer;
use num_prime::PrimalityTestConfig;

fn bench_showdown(c: &mut Criterion) {
    let test_cases = vec![
        ("127-bit", "170141183460469231731687303715884105727"),
        ("521-bit", "6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151"),
        ("1024-bit", "32317006071311007300714876688669951960444102669715484032130345427524655138867890893197201411522913463688717960921898019494119559150490921095088152386448283120630877367300996091750197750389652106796057638384067568276792218642619756161838094338476170470581645852036305042887575891541065808607552399123930385521914333389668342420684974786564569494856176035326322058077805659331026192708460314150258592864177116725943603718461857357598351152301645904403697613233287231227125684710820209725157101726931323469678542580656697935045997268352998638215525166389647960126939249806625440700685819469589938384356951833568218188663"),
    ];

    let mut group = c.benchmark_group("BPSW");

    for (label, n_str) in test_cases {
        let n_num = BigUint::parse_bytes(n_str.as_bytes(), 10).unwrap();
        let n_rug = Integer::from_str_radix(n_str, 10).unwrap();

        // 你的实现
        group.bench_with_input(BenchmarkId::new("My-BPSW", label), &n_num, |b, n| {
            b.iter(|| bpsw::is_prime(black_box(n)))
        });

        // num-prime 实现 25 轮 MR + Lucas
        let mut np_config = PrimalityTestConfig::default();

        // 设定为 25 轮确定性底数 MR 测试 (从底数 2 开始)
        np_config.sprp_trials = 1;
        // 随机底数设为 0，确保测试是确定性的
        np_config.sprp_random_trials = 0;
        // 关闭普通强 Lucas 测试
        np_config.slprp_test = false;
        // 开启“额外强 Lucas 测试”，以匹配你的算法 (Extra Strong Lucas)
        np_config.eslprp_test = true;

        group.bench_with_input(BenchmarkId::new("num-prime", label), &n_num, |b, n| {
            b.iter(|| num_prime::nt_funcs::is_prime(black_box(n), Some(np_config)))
        });

        //rug-GMP 实现 (等效于 BPSW + 追加轮次)
        group.bench_with_input(BenchmarkId::new("rug-GMP", label), &n_rug, |b, n| {
            b.iter(|| n.is_probably_prime(25))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_showdown);
criterion_main!(benches);


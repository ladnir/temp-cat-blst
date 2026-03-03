#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#if defined(BLST_AVX512)

#include <array>
#include <cstddef>
#include <random>
#include <vector>

#include "blst/fields/Fp2.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/Fp2_avx.h"
#include "blst/fields/Fp6_avx.h"
#include "blst/fields/Fp12_avx.h"

using namespace blst;

namespace
{
void random_redc_fp(Fp& v, std::mt19937_64& gen)
{
    for (u64 i = 0; i < v.num_limb(); ++i) {
        v[i] = gen();
    }
    while (Fp::is_redc(v) == 0) {
        Fp::rshift(v, v, 1);
    }
}

void random_redc_fp2(Fp2& v, std::mt19937_64& gen)
{
    random_redc_fp(v[0], gen);
    random_redc_fp(v[1], gen);
}

void random_redc_fp6(Fp6& v, std::mt19937_64& gen)
{
    random_redc_fp2(v[0], gen);
    random_redc_fp2(v[1], gen);
    random_redc_fp2(v[2], gen);
}

void random_redc_fp12(Fp12& v, std::mt19937_64& gen)
{
    for (std::size_t i = 0; i < 2; ++i) {
        random_redc_fp6(v[i], gen);
    }
}

void mul_fp6_scalar(Fp6& r, const Fp6& a, const Fp6& b)
{
    Fp2 t0{}, t1{}, t2{}, t3{}, t4{}, t5{};

    Fp2::mul(t0, a[0], b[0]);
    Fp2::mul(t1, a[1], b[1]);
    Fp2::mul(t2, a[2], b[2]);

    Fp2::add(t3, a[1], a[2]);
    Fp2::add(t4, b[1], b[2]);
    Fp2::mul(t5, t3, t4);
    Fp2::sub(t5, t5, t1);
    Fp2::sub(t5, t5, t2);
    Fp2::mul_by_u_plus_1(t5, t5);
    Fp2::add(r[0], t5, t0);

    Fp2::add(t3, a[0], a[1]);
    Fp2::add(t4, b[0], b[1]);
    Fp2::mul(t5, t3, t4);
    Fp2::sub(t5, t5, t0);
    Fp2::sub(t5, t5, t1);
    Fp2::mul_by_u_plus_1(t3, t2);
    Fp2::add(r[1], t5, t3);

    Fp2::add(t3, a[0], a[2]);
    Fp2::add(t4, b[0], b[2]);
    Fp2::mul(t5, t3, t4);
    Fp2::sub(t5, t5, t0);
    Fp2::sub(t5, t5, t2);
    Fp2::add(r[2], t5, t1);
}

void sqr_fp6_scalar(Fp6& r, const Fp6& a)
{
    mul_fp6_scalar(r, a, a);
}

void sqr_fp4_pack_bench(
    Fp2Pack8x1x1& r0,
    Fp2Pack8x1x1& r1,
    const Fp2Pack8x1x1& a0,
    const Fp2Pack8x1x1& a1)
{
    Fp2Pack8x1x1 t0{};
    Fp2Pack8x1x1 t1{};
    Fp2Pack8x1x1 t2{};

    Fp2Pack8x1x1::sqr(t0, a0);
    Fp2Pack8x1x1::sqr(t1, a1);

    Fp2Pack8x1x1::add(r1, a0, a1);
    Fp2Pack8x1x1::mul_by_u_plus_1(r0, t1);
    Fp2Pack8x1x1::add(r0, r0, t0);

    Fp2Pack8x1x1::sqr(t2, r1);
    Fp2Pack8x1x1::sub(t2, t2, t0);
    Fp2Pack8x1x1::sub(r1, t2, t1);
}

void ufp2_sqr_v2(UFp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
{
    UFp2Pack8x1x1::sqr(r, a);
}

void sqr_fp4_pack_bench_v2(
    Fp2Pack8x1x1& r0,
    Fp2Pack8x1x1& r1,
    const Fp2Pack8x1x1& a0,
    const Fp2Pack8x1x1& a1)
{
    UFp2Pack8x1x1 t0{};
    UFp2Pack8x1x1 t1{};
    UFp2Pack8x1x1 m01{};
    UFp2Pack8x1x1 t3{};

    ufp2_sqr_v2(t0, a0);         // a0^2
    ufp2_sqr_v2(t1, a1);         // a1^2
    UFp2Pack8x1x1::mul(m01, a0, a1);

    // r0 = a0^2 + (u+1)*a1^2
    UFp2Pack8x1x1::mul_by_u_plus_1(t3, t1);
    UFp2Pack8x1x1::add(t3, t3, t0);
    UFp2Pack8x1x1::redc(r0, t3);

    // r1 = 2*a0*a1
    UFp2Pack8x1x1::add(m01, m01, m01);
    UFp2Pack8x1x1::redc(r1, m01);
}

void sqr_fp4_pack_bench_v3(
    Fp2Pack8x1x1& r0,
    Fp2Pack8x1x1& r1,
    const Fp2Pack8x1x1& a0,
    const Fp2Pack8x1x1& a1)
{
    Fp2Pack8x1x1 s01{};
    UFp2Pack8x1x1 t0{};
    UFp2Pack8x1x1 t1{};
    UFp2Pack8x1x1 t2{};
    UFp2Pack8x1x1 t3{};

    Fp2Pack8x1x1::add(s01, a0, a1);
    UFp2Pack8x1x1::sqr(t0, a0);
    UFp2Pack8x1x1::sqr(t1, a1);
    UFp2Pack8x1x1::sqr(t2, s01);

    // r0 = a0^2 + (u+1)*a1^2
    UFp2Pack8x1x1::mul_by_u_plus_1(t3, t1);
    UFp2Pack8x1x1::add(t3, t3, t0);
    UFp2Pack8x1x1::redc(r0, t3);

    // r1 = (a0+a1)^2 - a0^2 - a1^2 = 2*a0*a1
    UFp2Pack8x1x1::sub(t2, t2, t0);
    UFp2Pack8x1x1::sub(t2, t2, t1);
    UFp2Pack8x1x1::redc(r1, t2);
}

void ufp2_mul_schoolbook(UFp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
{
    UFpPack8x1 a0b0{};
    UFpPack8x1 a1b1{};
    UFpPack8x1 a0b1{};
    UFpPack8x1 a1b0{};

    FpPack8x1::mul(a0b0, a[0], b[0]);
    FpPack8x1::mul(a1b1, a[1], b[1]);
    FpPack8x1::mul(a0b1, a[0], b[1]);
    FpPack8x1::mul(a1b0, a[1], b[0]);

    UFpPack8x1::sub(r[0], a0b0, a1b1);
    UFpPack8x1::add(r[1], a0b1, a1b0);
}

void sqr_fp4_pack_bench_v4(
    Fp2Pack8x1x1& r0,
    Fp2Pack8x1x1& r1,
    const Fp2Pack8x1x1& a0,
    const Fp2Pack8x1x1& a1)
{
    UFp2Pack8x1x1 t0{};
    UFp2Pack8x1x1 t1{};
    UFp2Pack8x1x1 m01{};
    UFp2Pack8x1x1 t3{};

    UFp2Pack8x1x1::sqr(t0, a0); // a0^2
    UFp2Pack8x1x1::sqr(t1, a1); // a1^2
    ufp2_mul_schoolbook(m01, a0, a1);

    // r0 = a0^2 + (u+1)*a1^2
    UFp2Pack8x1x1::mul_by_u_plus_1(t3, t1);
    UFp2Pack8x1x1::add(t3, t3, t0);
    UFp2Pack8x1x1::redc(r0, t3);

    // r1 = 2*a0*a1
    UFp2Pack8x1x1::add(m01, m01, m01);
    UFp2Pack8x1x1::redc(r1, m01);
}

struct CyclotomicSqrTerms8x1
{
    Fp2Pack8x1x1 t0_0{};
    Fp2Pack8x1x1 t0_1{};
    Fp2Pack8x1x1 t1_0{};
    Fp2Pack8x1x1 t1_1{};
    Fp2Pack8x1x1 t2_0{};
    Fp2Pack8x1x1 t2_1{};
};

void cyclotomic_sqr_terms_only(CyclotomicSqrTerms8x1& t, const Fp12Pack8x1& a)
{
    sqr_fp4_pack_bench(t.t0_0, t.t0_1, a.c[0].c[0], a.c[1].c[1]);
    sqr_fp4_pack_bench(t.t2_0, t.t2_1, a.c[0].c[1], a.c[1].c[2]);
    sqr_fp4_pack_bench(t.t1_0, t.t1_1, a.c[1].c[0], a.c[0].c[2]);
}

void cyclotomic_sqr_terms_only_v2(CyclotomicSqrTerms8x1& t, const Fp12Pack8x1& a)
{
    sqr_fp4_pack_bench_v2(t.t0_0, t.t0_1, a.c[0].c[0], a.c[1].c[1]);
    sqr_fp4_pack_bench_v2(t.t2_0, t.t2_1, a.c[0].c[1], a.c[1].c[2]);
    sqr_fp4_pack_bench_v2(t.t1_0, t.t1_1, a.c[1].c[0], a.c[0].c[2]);
}

void cyclotomic_sqr_terms_only_v3(CyclotomicSqrTerms8x1& t, const Fp12Pack8x1& a)
{
    sqr_fp4_pack_bench_v3(t.t0_0, t.t0_1, a.c[0].c[0], a.c[1].c[1]);
    sqr_fp4_pack_bench_v3(t.t2_0, t.t2_1, a.c[0].c[1], a.c[1].c[2]);
    sqr_fp4_pack_bench_v3(t.t1_0, t.t1_1, a.c[1].c[0], a.c[0].c[2]);
}

void cyclotomic_sqr_terms_only_v4(CyclotomicSqrTerms8x1& t, const Fp12Pack8x1& a)
{
    sqr_fp4_pack_bench_v4(t.t0_0, t.t0_1, a.c[0].c[0], a.c[1].c[1]);
    sqr_fp4_pack_bench_v4(t.t2_0, t.t2_1, a.c[0].c[1], a.c[1].c[2]);
    sqr_fp4_pack_bench_v4(t.t1_0, t.t1_1, a.c[1].c[0], a.c[0].c[2]);
}

void cyclotomic_sqr_combine_only(Fp12Pack8x1& r, const Fp12Pack8x1& a, const CyclotomicSqrTerms8x1& t)
{
    Fp2Pack8x1x1 tmp{};

    Fp2Pack8x1x1::sub(r.c[0].c[0], t.t0_0, a.c[0].c[0]);
    Fp2Pack8x1x1::add(r.c[0].c[0], r.c[0].c[0], r.c[0].c[0]);
    Fp2Pack8x1x1::add(r.c[0].c[0], r.c[0].c[0], t.t0_0);

    Fp2Pack8x1x1::add(r.c[1].c[1], t.t0_1, a.c[1].c[1]);
    Fp2Pack8x1x1::add(r.c[1].c[1], r.c[1].c[1], r.c[1].c[1]);
    Fp2Pack8x1x1::add(r.c[1].c[1], r.c[1].c[1], t.t0_1);

    Fp2Pack8x1x1::sub(r.c[0].c[1], t.t1_0, a.c[0].c[1]);
    Fp2Pack8x1x1::add(r.c[0].c[1], r.c[0].c[1], r.c[0].c[1]);
    Fp2Pack8x1x1::add(r.c[0].c[1], r.c[0].c[1], t.t1_0);

    Fp2Pack8x1x1::add(r.c[1].c[2], t.t1_1, a.c[1].c[2]);
    Fp2Pack8x1x1::add(r.c[1].c[2], r.c[1].c[2], r.c[1].c[2]);
    Fp2Pack8x1x1::add(r.c[1].c[2], r.c[1].c[2], t.t1_1);

    Fp2Pack8x1x1::sub(r.c[0].c[2], t.t2_0, a.c[0].c[2]);
    Fp2Pack8x1x1::add(r.c[0].c[2], r.c[0].c[2], r.c[0].c[2]);
    Fp2Pack8x1x1::add(r.c[0].c[2], r.c[0].c[2], t.t2_0);

    Fp2Pack8x1x1::mul_by_u_plus_1(tmp, t.t2_1);
    Fp2Pack8x1x1::add(r.c[1].c[0], tmp, a.c[1].c[0]);
    Fp2Pack8x1x1::add(r.c[1].c[0], r.c[1].c[0], r.c[1].c[0]);
    Fp2Pack8x1x1::add(r.c[1].c[0], r.c[1].c[0], tmp);
}

void cyclotomic_sqr_pack_legacy(Fp12Pack8x1& r, const Fp12Pack8x1& a)
{
    CyclotomicSqrTerms8x1 t{};
    cyclotomic_sqr_terms_only(t, a);
    cyclotomic_sqr_combine_only(r, a, t);
}

void cyclotomic_sqr_pack_v4(Fp12Pack8x1& r, const Fp12Pack8x1& a)
{
    CyclotomicSqrTerms8x1 t{};
    cyclotomic_sqr_terms_only_v4(t, a);
    cyclotomic_sqr_combine_only(r, a, t);
}

void mul_fp12_pack_legacy(Fp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b)
{
    Fp6Pack8x1 t0{};
    Fp6Pack8x1 t1{};
    Fp6Pack8x1 t2{};
    Fp6Pack8x1 t3{};

    Fp6Pack8x1::mul(t0, a.c[0], b.c[0]);
    Fp6Pack8x1::mul(t1, a.c[1], b.c[1]);

    Fp6Pack8x1::add(t2, a.c[0], a.c[1]);
    Fp6Pack8x1::add(t3, b.c[0], b.c[1]);
    Fp6Pack8x1::mul(r.c[1], t2, t3);
    Fp6Pack8x1::sub(r.c[1], r.c[1], t0);
    Fp6Pack8x1::sub(r.c[1], r.c[1], t1);

    Fp6Pack8x1::mul_by_nonresidue(t2, t1);
    Fp6Pack8x1::add(r.c[0], t0, t2);
}

void sqr_fp12_pack_legacy(Fp12Pack8x1& r, const Fp12Pack8x1& a)
{
    Fp6Pack8x1 s0{};
    Fp6Pack8x1 s1{};
    Fp6Pack8x1 m01{};
    Fp6Pack8x1 t0{};

    Fp6Pack8x1::sqr(s0, a.c[0]);
    Fp6Pack8x1::sqr(s1, a.c[1]);
    Fp6Pack8x1::mul(m01, a.c[0], a.c[1]);
    Fp6Pack8x1::add(r.c[1], m01, m01);
    Fp6Pack8x1::mul_by_nonresidue(t0, s1);
    Fp6Pack8x1::add(r.c[0], s0, t0);
}

struct Fp2BenchData
{
    explicit Fp2BenchData(const std::size_t packs)
        : n_packs(packs)
        , n_tuples(packs * 8)
        , a(n_tuples)
        , b(n_tuples)
        , out_scalar(n_tuples)
        , a_pack(n_packs)
        , b_pack(n_packs)
        , out_pack(n_packs)
    {
        std::mt19937_64 gen(0xF2001234ULL + static_cast<u64>(packs));
        for (std::size_t i = 0; i < n_tuples; ++i) {
            random_redc_fp2(a[i], gen);
            random_redc_fp2(b[i], gen);
        }
        for (std::size_t p = 0; p < n_packs; ++p) {
            Fp2Pack8x1x1::from_aos(a_pack[p], &a[p * 8]);
            Fp2Pack8x1x1::from_aos(b_pack[p], &b[p * 8]);
        }
    }

    const std::size_t n_packs;
    const std::size_t n_tuples;
    std::vector<Fp2> a;
    std::vector<Fp2> b;
    std::vector<Fp2> out_scalar;
    std::vector<Fp2Pack8x1x1> a_pack;
    std::vector<Fp2Pack8x1x1> b_pack;
    std::vector<Fp2Pack8x1x1> out_pack;
};

struct Fp6BenchData
{
    explicit Fp6BenchData(const std::size_t packs)
        : n_packs(packs)
        , n_tuples(packs * 8)
        , a(n_tuples)
        , b(n_tuples)
        , out_scalar(n_tuples)
        , a_pack(n_packs)
        , b_pack(n_packs)
        , out_pack(n_packs)
    {
        std::mt19937_64 gen(0xF6001234ULL + static_cast<u64>(packs));
        for (std::size_t i = 0; i < n_tuples; ++i) {
            random_redc_fp6(a[i], gen);
            random_redc_fp6(b[i], gen);
        }
        for (std::size_t p = 0; p < n_packs; ++p) {
            Fp6Pack8x1::from_aos(a_pack[p], &a[p * 8]);
            Fp6Pack8x1::from_aos(b_pack[p], &b[p * 8]);
        }
    }

    const std::size_t n_packs;
    const std::size_t n_tuples;
    std::vector<Fp6> a;
    std::vector<Fp6> b;
    std::vector<Fp6> out_scalar;
    std::vector<Fp6Pack8x1> a_pack;
    std::vector<Fp6Pack8x1> b_pack;
    std::vector<Fp6Pack8x1> out_pack;
};

struct Fp12BenchData
{
    explicit Fp12BenchData(const std::size_t packs)
        : n_packs(packs)
        , n_tuples(packs * 8)
        , a(n_tuples)
        , b(n_tuples)
        , cyclo(n_tuples)
        , out_scalar(n_tuples)
        , a_pack(n_packs)
        , b_pack(n_packs)
        , cyclo_pack(n_packs)
        , wide_mul_pack(n_packs)
        , wide_sqr_pack(n_packs)
        , cyclo_terms(n_packs)
        , out_pack(n_packs)
    {
        std::mt19937_64 gen(0xF1201234ULL + static_cast<u64>(packs));
        for (std::size_t i = 0; i < n_tuples; ++i) {
            random_redc_fp12(a[i], gen);
            random_redc_fp12(b[i], gen);
            Fp12 seed{};
            random_redc_fp12(seed, gen);
            Fp12::final_exp(cyclo[i], seed);
        }
        for (std::size_t p = 0; p < n_packs; ++p) {
            Fp12Pack8x1::from_aos(a_pack[p], &a[p * 8]);
            Fp12Pack8x1::from_aos(b_pack[p], &b[p * 8]);
            Fp12Pack8x1::from_aos(cyclo_pack[p], &cyclo[p * 8]);
            Fp12Pack8x1::mul(wide_mul_pack[p], a_pack[p], b_pack[p]);
            Fp12Pack8x1::sqr(wide_sqr_pack[p], a_pack[p]);
            cyclotomic_sqr_terms_only(cyclo_terms[p], cyclo_pack[p]);
        }
    }

    const std::size_t n_packs;
    const std::size_t n_tuples;
    std::vector<Fp12> a;
    std::vector<Fp12> b;
    std::vector<Fp12> cyclo;
    std::vector<Fp12> out_scalar;
    std::vector<Fp12Pack8x1> a_pack;
    std::vector<Fp12Pack8x1> b_pack;
    std::vector<Fp12Pack8x1> cyclo_pack;
    std::vector<UFp12Pack8x1> wide_mul_pack;
    std::vector<UFp12Pack8x1> wide_sqr_pack;
    std::vector<CyclotomicSqrTerms8x1> cyclo_terms;
    std::vector<Fp12Pack8x1> out_pack;
};

void BM_Fp2Mul_Scalar8(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp2BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t i = 0; i < data.n_tuples; ++i) {
            Fp2::mul(data.out_scalar[i], data.a[i], data.b[i]);
        }
        benchmark::DoNotOptimize(data.out_scalar.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp2Mul_Scalar8)->Arg(256);

void BM_Fp2Mul_Pack8x1(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp2BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            Fp2Pack8x1x1::mul(data.out_pack[p], data.a_pack[p], data.b_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp2Mul_Pack8x1)->Arg(256);

void BM_Fp2Sqr_Scalar8(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp2BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t i = 0; i < data.n_tuples; ++i) {
            Fp2::sqr(data.out_scalar[i], data.a[i]);
        }
        benchmark::DoNotOptimize(data.out_scalar.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp2Sqr_Scalar8)->Arg(256);

void BM_Fp2Sqr_Pack8x1(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp2BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            Fp2Pack8x1x1::sqr(data.out_pack[p], data.a_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp2Sqr_Pack8x1)->Arg(256);

void BM_Fp6Mul_Scalar8(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp6BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t i = 0; i < data.n_tuples; ++i) {
            mul_fp6_scalar(data.out_scalar[i], data.a[i], data.b[i]);
        }
        benchmark::DoNotOptimize(data.out_scalar.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp6Mul_Scalar8)->Arg(256);

void BM_Fp6Mul_Pack8x1(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp6BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            Fp6Pack8x1::mul(data.out_pack[p], data.a_pack[p], data.b_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp6Mul_Pack8x1)->Arg(256);

void BM_Fp6Sqr_Scalar8(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp6BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t i = 0; i < data.n_tuples; ++i) {
            sqr_fp6_scalar(data.out_scalar[i], data.a[i]);
        }
        benchmark::DoNotOptimize(data.out_scalar.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp6Sqr_Scalar8)->Arg(256);

void BM_Fp6Sqr_Pack8x1(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp6BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            Fp6Pack8x1::sqr(data.out_pack[p], data.a_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp6Sqr_Pack8x1)->Arg(256);

void BM_Fp12Mul_Scalar8(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t i = 0; i < data.n_tuples; ++i) {
            Fp12::mul_scalar(data.out_scalar[i], data.a[i], data.b[i]);
        }
        benchmark::DoNotOptimize(data.out_scalar.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Mul_Scalar8)->Arg(128);

void BM_Fp12Mul_Pack8x1(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            Fp12Pack8x1::mul(data.out_pack[p], data.a_pack[p], data.b_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Mul_Pack8x1)->Arg(128);

void BM_Fp12Mul_Pack8x1_Legacy(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            mul_fp12_pack_legacy(data.out_pack[p], data.a_pack[p], data.b_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Mul_Pack8x1_Legacy)->Arg(128);

void BM_Fp12Mul_Pack8x1_WideOnly(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            Fp12Pack8x1::mul(data.wide_mul_pack[p], data.a_pack[p], data.b_pack[p]);
        }
        benchmark::DoNotOptimize(data.wide_mul_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Mul_Pack8x1_WideOnly)->Arg(128);

void BM_Fp12Mul_Pack8x1_RedcOnly(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            UFp12Pack8x1::redc(data.out_pack[p], data.wide_mul_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Mul_Pack8x1_RedcOnly)->Arg(128);

void BM_Fp12Sqr_Scalar8(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t i = 0; i < data.n_tuples; ++i) {
            Fp12::sqr_scalar(data.out_scalar[i], data.a[i]);
        }
        benchmark::DoNotOptimize(data.out_scalar.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Sqr_Scalar8)->Arg(128);

void BM_Fp12Sqr_Pack8x1(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            Fp12Pack8x1::sqr(data.out_pack[p], data.a_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Sqr_Pack8x1)->Arg(128);

void BM_Fp12Sqr_Pack8x1_Legacy(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            sqr_fp12_pack_legacy(data.out_pack[p], data.a_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Sqr_Pack8x1_Legacy)->Arg(128);

void BM_Fp12Sqr_Pack8x1_WideOnly(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            Fp12Pack8x1::sqr(data.wide_sqr_pack[p], data.a_pack[p]);
        }
        benchmark::DoNotOptimize(data.wide_sqr_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Sqr_Pack8x1_WideOnly)->Arg(128);

void BM_Fp12Sqr_Pack8x1_RedcOnly(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            UFp12Pack8x1::redc(data.out_pack[p], data.wide_sqr_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Sqr_Pack8x1_RedcOnly)->Arg(128);

void BM_Fp12CyclotomicSqr_Scalar8(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t i = 0; i < data.n_tuples; ++i) {
            Fp12::cyclotomic_sqr_scalar(data.out_scalar[i], data.cyclo[i]);
        }
        benchmark::DoNotOptimize(data.out_scalar.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12CyclotomicSqr_Scalar8)->Arg(128);

void BM_Fp12CyclotomicSqr_Pack8x1(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            Fp12Pack8x1::cyclotomic_sqr(data.out_pack[p], data.cyclo_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12CyclotomicSqr_Pack8x1)->Arg(128);

void BM_Fp12CyclotomicSqr_Pack8x1_V4(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            cyclotomic_sqr_pack_v4(data.out_pack[p], data.cyclo_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12CyclotomicSqr_Pack8x1_V4)->Arg(128);

void BM_Fp12CyclotomicSqr_Pack8x1_Legacy(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            cyclotomic_sqr_pack_legacy(data.out_pack[p], data.cyclo_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12CyclotomicSqr_Pack8x1_Legacy)->Arg(128);

void BM_Fp12CyclotomicSqr_Pack8x1_SqrFp4Only(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);
    std::vector<CyclotomicSqrTerms8x1> work_terms(data.n_packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            cyclotomic_sqr_terms_only(work_terms[p], data.cyclo_pack[p]);
        }
        benchmark::DoNotOptimize(work_terms.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12CyclotomicSqr_Pack8x1_SqrFp4Only)->Arg(128);

void BM_Fp12CyclotomicSqr_Pack8x1_SqrFp4Only_V2(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);
    std::vector<CyclotomicSqrTerms8x1> work_terms(data.n_packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            cyclotomic_sqr_terms_only_v2(work_terms[p], data.cyclo_pack[p]);
        }
        benchmark::DoNotOptimize(work_terms.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12CyclotomicSqr_Pack8x1_SqrFp4Only_V2)->Arg(128);

void BM_Fp12CyclotomicSqr_Pack8x1_SqrFp4Only_V3(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);
    std::vector<CyclotomicSqrTerms8x1> work_terms(data.n_packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            cyclotomic_sqr_terms_only_v3(work_terms[p], data.cyclo_pack[p]);
        }
        benchmark::DoNotOptimize(work_terms.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12CyclotomicSqr_Pack8x1_SqrFp4Only_V3)->Arg(128);

void BM_Fp12CyclotomicSqr_Pack8x1_SqrFp4Only_V4(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);
    std::vector<CyclotomicSqrTerms8x1> work_terms(data.n_packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            cyclotomic_sqr_terms_only_v4(work_terms[p], data.cyclo_pack[p]);
        }
        benchmark::DoNotOptimize(work_terms.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12CyclotomicSqr_Pack8x1_SqrFp4Only_V4)->Arg(128);

void BM_Fp12CyclotomicSqr_Pack8x1_CombineOnly(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            cyclotomic_sqr_combine_only(data.out_pack[p], data.cyclo_pack[p], data.cyclo_terms[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12CyclotomicSqr_Pack8x1_CombineOnly)->Arg(128);

} // namespace

#endif // BLST_AVX512

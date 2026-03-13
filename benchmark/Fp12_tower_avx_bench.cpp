#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#if defined(BLST_AVX512)

#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <random>
#include <vector>

#include "blst/fields/Fp2.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/Fp2_avx.h"
#include "blst/fields/Fp6_avx.h"
#include "blst/fields/Fp12_avx.h"
#include "benchmark/internal/Fp_avx_lufp_v2_impl.h"

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

bool fp2_equal(const Fp2& a, const Fp2& b)
{
    return Fp::equal(a[0], b[0]) && Fp::equal(a[1], b[1]);
}

void random_redc_fp6(Fp6& v, std::mt19937_64& gen)
{
    random_redc_fp2(v[0], gen);
    random_redc_fp2(v[1], gen);
    random_redc_fp2(v[2], gen);
}

bool fp6_equal(const Fp6& a, const Fp6& b)
{
    return fp2_equal(a[0], b[0]) && fp2_equal(a[1], b[1]) && fp2_equal(a[2], b[2]);
}

bool fp12_equal(const Fp12& a, const Fp12& b)
{
    return fp6_equal(a[0], b[0]) && fp6_equal(a[1], b[1]);
}

void dump_fp(const char* label, const Fp& v)
{
    std::fprintf(
        stderr,
        "%s = [%016llx %016llx %016llx %016llx %016llx %016llx]\n",
        label,
        static_cast<unsigned long long>(v[5]),
        static_cast<unsigned long long>(v[4]),
        static_cast<unsigned long long>(v[3]),
        static_cast<unsigned long long>(v[2]),
        static_cast<unsigned long long>(v[1]),
        static_cast<unsigned long long>(v[0]));
}

void dump_fp2(const char* label, const Fp2& v)
{
    std::fprintf(stderr, "%s.re\n", label);
    dump_fp("  ", v[0]);
    std::fprintf(stderr, "%s.im\n", label);
    dump_fp("  ", v[1]);
}

void dump_fp6(const char* label, const Fp6& v)
{
    std::fprintf(stderr, "%s.c0\n", label);
    dump_fp2("  ", v[0]);
    std::fprintf(stderr, "%s.c1\n", label);
    dump_fp2("  ", v[1]);
    std::fprintf(stderr, "%s.c2\n", label);
    dump_fp2("  ", v[2]);
}

void dump_fp12(const char* label, const Fp12& v)
{
    std::fprintf(stderr, "%s.c0\n", label);
    dump_fp6("  ", v[0]);
    std::fprintf(stderr, "%s.c1\n", label);
    dump_fp6("  ", v[1]);
}

void random_redc_fp12(Fp12& v, std::mt19937_64& gen)
{
    for (std::size_t i = 0; i < 2; ++i) {
        random_redc_fp6(v[i], gen);
    }
}

struct BenchLUFp2Pack8x1x1
{
    std::array<LUFpPack8x1, 2> c;

    BLST_FORCEINLINE LUFpPack8x1& operator[](std::size_t i)
    {
        return c[i];
    }

    BLST_FORCEINLINE const LUFpPack8x1& operator[](std::size_t i) const
    {
        return c[i];
    }
};

void mul_fp6_pack_baseline_no_repair_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b);
void sqr_fp6_pack_baseline_no_repair_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a);

BLST_FORCEINLINE void bench_lufp8x1_bias_p(LUFpPack8x1& r, const LUFpPack8x1& a)
{
    for (std::size_t i = 0; i < NWORDS; ++i) {
        r[i] = a[i];
    }
    r[8] = VADD(a[8], VSET1(static_cast<i64>(P48[0])));
    r[9] = VADD(a[9], VSET1(static_cast<i64>(P48[1])));
    r[10] = VADD(a[10], VSET1(static_cast<i64>(P48[2])));
    r[11] = VADD(a[11], VSET1(static_cast<i64>(P48[3])));
    r[12] = VADD(a[12], VSET1(static_cast<i64>(P48[4])));
    r[13] = VADD(a[13], VSET1(static_cast<i64>(P48[5])));
    r[14] = VADD(a[14], VSET1(static_cast<i64>(P48[6])));
    r[15] = VADD(a[15], VSET1(static_cast<i64>(P48[7])));
}

BLST_FORCEINLINE void bench_lufp2_add(
    BenchLUFp2Pack8x1x1& r,
    const BenchLUFp2Pack8x1x1& a,
    const BenchLUFp2Pack8x1x1& b)
{
    LUFpPack8x1::add(r[0], a[0], b[0]);
    LUFpPack8x1::add(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_lufp2_add(
    BenchLUFp2Pack8x1x1& r,
    const BenchLUFp2Pack8x1x1& a,
    const UFp2Pack8x1x1& b)
{
    LUFpPack8x1::add(r[0], a[0], b[0]);
    LUFpPack8x1::add(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_lufp2_add(
    BenchLUFp2Pack8x1x1& r,
    const UFp2Pack8x1x1& a,
    const BenchLUFp2Pack8x1x1& b)
{
    LUFpPack8x1::add(r[0], a[0], b[0]);
    LUFpPack8x1::add(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_lufp2_add(
    BenchLUFp2Pack8x1x1& r,
    const UFp2Pack8x1x1& a,
    const UFp2Pack8x1x1& b)
{
    LUFpPack8x1::add(r[0], a[0], b[0]);
    LUFpPack8x1::add(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_lufp2_sub(
    BenchLUFp2Pack8x1x1& r,
    const BenchLUFp2Pack8x1x1& a,
    const BenchLUFp2Pack8x1x1& b)
{
    LUFpPack8x1::sub(r[0], a[0], b[0]);
    LUFpPack8x1::sub(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_lufp2_sub(
    BenchLUFp2Pack8x1x1& r,
    const BenchLUFp2Pack8x1x1& a,
    const UFp2Pack8x1x1& b)
{
    LUFpPack8x1::sub(r[0], a[0], b[0]);
    LUFpPack8x1::sub(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_lufp2_sub(
    BenchLUFp2Pack8x1x1& r,
    const UFp2Pack8x1x1& a,
    const BenchLUFp2Pack8x1x1& b)
{
    LUFpPack8x1::sub(r[0], a[0], b[0]);
    LUFpPack8x1::sub(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_lufp2_sub(
    BenchLUFp2Pack8x1x1& r,
    const UFp2Pack8x1x1& a,
    const UFp2Pack8x1x1& b)
{
    LUFpPack8x1::sub(r[0], a[0], b[0]);
    LUFpPack8x1::sub(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_lufp2_bias_p(BenchLUFp2Pack8x1x1& r, const BenchLUFp2Pack8x1x1& a)
{
    bench_lufp8x1_bias_p(r[0], a[0]);
    bench_lufp8x1_bias_p(r[1], a[1]);
}

BLST_FORCEINLINE void bench_lufp2_normalize(UFp2Pack8x1x1& r, const BenchLUFp2Pack8x1x1& a)
{
    LUFpPack8x1::normalize<LUFpNormalizeBound::HiNeg1Pos2>(r[0], a[0]);
    LUFpPack8x1::normalize<LUFpNormalizeBound::HiNeg1Pos2>(r[1], a[1]);
}

BLST_FORCEINLINE void bench_lufp2_copy_as_ufp(UFp2Pack8x1x1& r, const BenchLUFp2Pack8x1x1& a)
{
    for (std::size_t c = 0; c < 2; ++c) {
        for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
            r[c][i] = a[c][i];
        }
    }
}

BLST_FORCEINLINE void bench_lufp8x1_copy_as_ufp(UFpPack8x1& r, const LUFpPack8x1& a)
{
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
        r[i] = a[i];
    }
}

BLST_FORCEINLINE void bench_ufp8x1_add_no_repair(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
        r[i] = VADD(a[i], b[i]);
    }
}

BLST_FORCEINLINE void bench_ufp8x1_sub_no_repair(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
        r[i] = VSUB(a[i], b[i]);
    }
}

BLST_FORCEINLINE void bench_ufp2_add_no_repair(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    bench_ufp8x1_add_no_repair(r[0], a[0], b[0]);
    bench_ufp8x1_add_no_repair(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_ufp2_sub_no_repair(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    bench_ufp8x1_sub_no_repair(r[0], a[0], b[0]);
    bench_ufp8x1_sub_no_repair(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_ufp2_mul_by_u_plus_1_no_repair(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
{
    bench_ufp8x1_sub_no_repair(r[0], a[0], a[1]);
    bench_ufp8x1_add_no_repair(r[1], a[0], a[1]);
}

BLST_FORCEINLINE void bench_lufp8x1_add_no_normalize(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
    LUFpPack8x1 lazy{};
    LUFpPack8x1::add(lazy, a, b);
    bench_lufp8x1_copy_as_ufp(r, lazy);
}

BLST_FORCEINLINE void bench_lufp8x1_sub_no_normalize(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
    LUFpPack8x1 lazy{};
    LUFpPack8x1::sub(lazy, a, b);
    bench_lufp8x1_copy_as_ufp(r, lazy);
}

BLST_FORCEINLINE void bench_lufp8x1_add_renorm_no_normalize(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
    LUFpPack8x1 lazy{};
    LUFpPack8x1 tmp{};
    LUFpPack8x1::add(lazy, a, b);
    LUFpPack8x1::renorm<LUFpRenormClass::Loose>(tmp, lazy);
    bench_lufp8x1_copy_as_ufp(r, tmp);
}

BLST_FORCEINLINE void bench_lufp8x1_sub_renorm_no_normalize(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
    LUFpPack8x1 lazy{};
    LUFpPack8x1 tmp{};
    LUFpPack8x1::sub(lazy, a, b);
    LUFpPack8x1::renorm<LUFpRenormClass::Loose>(tmp, lazy);
    bench_lufp8x1_copy_as_ufp(r, tmp);
}

BLST_FORCEINLINE void bench_lufp2_add_no_normalize(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    bench_lufp8x1_add_no_normalize(r[0], a[0], b[0]);
    bench_lufp8x1_add_no_normalize(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_lufp2_sub_no_normalize(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    bench_lufp8x1_sub_no_normalize(r[0], a[0], b[0]);
    bench_lufp8x1_sub_no_normalize(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_lufp2_mul_by_u_plus_1_no_normalize(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
{
    bench_lufp8x1_sub_no_normalize(r[0], a[0], a[1]);
    bench_lufp8x1_add_no_normalize(r[1], a[0], a[1]);
}

BLST_FORCEINLINE void bench_lufp2_add_renorm_no_normalize(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    bench_lufp8x1_add_renorm_no_normalize(r[0], a[0], b[0]);
    bench_lufp8x1_add_renorm_no_normalize(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_lufp2_sub_renorm_no_normalize(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    bench_lufp8x1_sub_renorm_no_normalize(r[0], a[0], b[0]);
    bench_lufp8x1_sub_renorm_no_normalize(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void bench_lufp2_mul_by_u_plus_1_renorm_no_normalize(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
{
    bench_lufp8x1_sub_renorm_no_normalize(r[0], a[0], a[1]);
    bench_lufp8x1_add_renorm_no_normalize(r[1], a[0], a[1]);
}

BLST_FORCEINLINE void bench_ufp2_mul_no_repair(UFp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
{
    FpPack8x1 t0, t1;
    UFpPack8x1 tt0, tt1, tt2;

    FpPack8x1::mul(tt0, a[0], b[0]);
    FpPack8x1::mul(tt1, a[1], b[1]);
    FpPack8x1::add(t0, a[0], a[1]);
    FpPack8x1::add(t1, b[0], b[1]);
    FpPack8x1::mul(tt2, t0, t1);
    bench_ufp8x1_sub_no_repair(r[0], tt0, tt1);
    bench_ufp8x1_sub_no_repair(tt2, tt2, tt0);
    bench_ufp8x1_sub_no_repair(r[1], tt2, tt1);
}

BLST_FORCEINLINE void bench_ufp2_sqr_no_repair(UFp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
{
    FpPack8x1 sum{};
    FpPack8x1 diff{};
    UFpPack8x1 t1{};

    FpPack8x1::add(sum, a[0], a[1]);
    FpPack8x1::sub(diff, a[0], a[1]);
    FpPack8x1::mul(r[0], sum, diff);

    FpPack8x1::mul(t1, a[0], a[1]);
    bench_ufp8x1_add_no_repair(r[1], t1, t1);
}

BLST_FORCEINLINE void bench_lufp2_mul_no_normalize(UFp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
{
    FpPack8x1 t0, t1;
    UFpPack8x1 tt0, tt1, tt2;

    FpPack8x1::mul(tt0, a[0], b[0]);
    FpPack8x1::mul(tt1, a[1], b[1]);
    FpPack8x1::add(t0, a[0], a[1]);
    FpPack8x1::add(t1, b[0], b[1]);
    FpPack8x1::mul(tt2, t0, t1);

    bench_lufp8x1_sub_no_normalize(r[0], tt0, tt1);
    bench_lufp8x1_sub_no_normalize(tt2, tt2, tt0);
    bench_lufp8x1_sub_no_normalize(r[1], tt2, tt1);
}

BLST_FORCEINLINE void bench_lufp2_sqr_no_normalize(UFp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
{
    FpPack8x1 sum{};
    FpPack8x1 diff{};
    UFpPack8x1 t1{};

    FpPack8x1::add(sum, a[0], a[1]);
    FpPack8x1::sub(diff, a[0], a[1]);
    FpPack8x1::mul(r[0], sum, diff);

    FpPack8x1::mul(t1, a[0], a[1]);
    bench_lufp8x1_add_no_normalize(r[1], t1, t1);
}

BLST_FORCEINLINE void bench_lufp2_mul_renorm_no_normalize(UFp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
{
    FpPack8x1 t0, t1;
    UFpPack8x1 tt0, tt1, tt2;

    FpPack8x1::mul(tt0, a[0], b[0]);
    FpPack8x1::mul(tt1, a[1], b[1]);
    FpPack8x1::add(t0, a[0], a[1]);
    FpPack8x1::add(t1, b[0], b[1]);
    FpPack8x1::mul(tt2, t0, t1);

    bench_lufp8x1_sub_renorm_no_normalize(r[0], tt0, tt1);
    bench_lufp8x1_sub_renorm_no_normalize(tt2, tt2, tt0);
    bench_lufp8x1_sub_renorm_no_normalize(r[1], tt2, tt1);
}

BLST_FORCEINLINE void bench_lufp2_sqr_renorm_no_normalize(UFp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
{
    FpPack8x1 sum{};
    FpPack8x1 diff{};
    UFpPack8x1 t1{};

    FpPack8x1::add(sum, a[0], a[1]);
    FpPack8x1::sub(diff, a[0], a[1]);
    FpPack8x1::mul(r[0], sum, diff);

    FpPack8x1::mul(t1, a[0], a[1]);
    bench_lufp8x1_add_renorm_no_normalize(r[1], t1, t1);
}

void mul_fp6_pack_lufp_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    Fp2Pack8x1x1 aa{}, bb{};
    UFp2Pack8x1x1 t0{}, t1{}, t2{}, wide{}, strict0{}, strict1{};
    BenchLUFp2Pack8x1x1 lazy0{}, lazy1{}, lazy2{};

    UFp2Pack8x1x1::mul(t0, a.c[0], b.c[0]);
    UFp2Pack8x1x1::mul(t1, a.c[1], b.c[1]);
    UFp2Pack8x1x1::mul(t2, a.c[2], b.c[2]);

    // c0 = ((a1+a2)*(b1+b2) - a1*b1 - a2*b2) * (u+1) + a0*b0
    Fp2Pack8x1x1::add(aa, a.c[1], a.c[2]);
    Fp2Pack8x1x1::add(bb, b.c[1], b.c[2]);
    UFp2Pack8x1x1::mul(wide, aa, bb);
    bench_lufp2_sub(lazy0, wide, t1);
    bench_lufp2_sub(lazy1, lazy0, t2);
    bench_lufp2_bias_p(lazy2, lazy1);
    bench_lufp2_normalize(strict0, lazy2);
    UFp2Pack8x1x1::mul_by_u_plus_1(strict1, strict0);
    UFp2Pack8x1x1::add(r.c[0], strict1, t0);

    // c1 = (a0+a1)*(b0+b1) - a0*b0 - a1*b1 + (u+1)*a2*b2
    Fp2Pack8x1x1::add(aa, a.c[0], a.c[1]);
    Fp2Pack8x1x1::add(bb, b.c[0], b.c[1]);
    UFp2Pack8x1x1::mul(wide, aa, bb);
    bench_lufp2_sub(lazy0, wide, t0);
    bench_lufp2_sub(lazy1, lazy0, t1);
    bench_lufp2_bias_p(lazy2, lazy1);
    bench_lufp2_normalize(strict0, lazy2);
    UFp2Pack8x1x1::mul_by_u_plus_1(strict1, t2);
    UFp2Pack8x1x1::add(r.c[1], strict0, strict1);

    // c2 = (a0+a2)*(b0+b2) - a0*b0 - a2*b2 + a1*b1
    Fp2Pack8x1x1::add(aa, a.c[0], a.c[2]);
    Fp2Pack8x1x1::add(bb, b.c[0], b.c[2]);
    UFp2Pack8x1x1::mul(wide, aa, bb);
    bench_lufp2_sub(lazy0, wide, t0);
    bench_lufp2_sub(lazy1, lazy0, t2);
    bench_lufp2_bias_p(lazy2, lazy1);
    bench_lufp2_normalize(strict0, lazy2);
    UFp2Pack8x1x1::add(r.c[2], strict0, t1);
}

void mul_fp6_pack_lufp(Fp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    UFp6Pack8x1 wide{};
    mul_fp6_pack_lufp_wide(wide, a, b);
    UFp6Pack8x1::redc(r, wide);
}

void sqr_fp6_pack_lufp_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a)
{
    Fp2Pack8x1x1 sum{};
    UFp2Pack8x1x1 s0{}, p01{}, p12{}, s2{}, wide{}, m01{}, m12{}, strict0{}, strict1{};
    BenchLUFp2Pack8x1x1 lazy0{}, lazy1{}, lazy2{};

    UFp2Pack8x1x1::mul(s0, a.c[0], a.c[0]);
    UFp2Pack8x1x1::mul(p01, a.c[0], a.c[1]);
    UFp2Pack8x1x1::mul(p12, a.c[1], a.c[2]);
    UFp2Pack8x1x1::mul(s2, a.c[2], a.c[2]);

    // Double the cross terms without re-entering the strict widened contract.
    bench_lufp2_add(lazy0, p01, p01);
    bench_lufp2_normalize(m01, lazy0);
    bench_lufp2_add(lazy0, p12, p12);
    bench_lufp2_normalize(m12, lazy0);

    // c2 = (a0+a1+a2)^2 - a0^2 - a2^2 - 2*a0*a1 - 2*a1*a2
    Fp2Pack8x1x1::add(sum, a.c[2], a.c[1]);
    Fp2Pack8x1x1::add(sum, sum, a.c[0]);
    UFp2Pack8x1x1::mul(wide, sum, sum);
    bench_lufp2_sub(lazy0, wide, s0);
    bench_lufp2_sub(lazy1, lazy0, s2);
    bench_lufp2_bias_p(lazy2, lazy1);
    bench_lufp2_normalize(strict0, lazy2);
    bench_lufp2_sub(lazy0, strict0, m01);
    bench_lufp2_sub(lazy1, lazy0, m12);
    bench_lufp2_bias_p(lazy2, lazy1);
    bench_lufp2_normalize(r.c[2], lazy2);

    // c0 = a0^2 + 2*a1*a2*(u+1)
    UFp2Pack8x1x1::mul_by_u_plus_1(strict1, m12);
    UFp2Pack8x1x1::add(r.c[0], strict1, s0);

    // c1 = a2^2*(u+1) + 2*a0*a1
    UFp2Pack8x1x1::mul_by_u_plus_1(strict1, s2);
    UFp2Pack8x1x1::add(r.c[1], strict1, m01);
}

void sqr_fp6_pack_lufp(Fp6Pack8x1& r, const Fp6Pack8x1& a)
{
    UFp6Pack8x1 wide{};
    sqr_fp6_pack_lufp_wide(wide, a);
    UFp6Pack8x1::redc(r, wide);
}

void mul_fp6_pack_lufp_no_normalize_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    // Parity experiment: with renorm/normalize disabled, the lazy path should
    // collapse to the same no-repair widened schedule as the strict baseline.
    mul_fp6_pack_baseline_no_repair_wide(r, a, b);
}

void mul_fp6_pack_lufp_no_normalize(Fp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    UFp6Pack8x1 wide{};
    mul_fp6_pack_lufp_no_normalize_wide(wide, a, b);
    UFp6Pack8x1::redc(r, wide);
}

void sqr_fp6_pack_lufp_no_normalize_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a)
{
    // Parity experiment: with renorm/normalize disabled, the lazy path should
    // collapse to the same no-repair widened schedule as the strict baseline.
    sqr_fp6_pack_baseline_no_repair_wide(r, a);
}

void sqr_fp6_pack_lufp_no_normalize(Fp6Pack8x1& r, const Fp6Pack8x1& a)
{
    UFp6Pack8x1 wide{};
    sqr_fp6_pack_lufp_no_normalize_wide(wide, a);
    UFp6Pack8x1::redc(r, wide);
}

void mul_fp6_pack_lufp_renorm_no_normalize_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    Fp2Pack8x1x1 aa{}, bb{};
    UFp2Pack8x1x1 t0{}, t1{}, t2{}, t3{};

    bench_lufp2_mul_renorm_no_normalize(t0, a.c[0], b.c[0]);
    bench_lufp2_mul_renorm_no_normalize(t1, a.c[1], b.c[1]);
    bench_lufp2_mul_renorm_no_normalize(t2, a.c[2], b.c[2]);

    Fp2Pack8x1x1::add(aa, a.c[1], a.c[2]);
    Fp2Pack8x1x1::add(bb, b.c[1], b.c[2]);
    bench_lufp2_mul_renorm_no_normalize(r.c[0], aa, bb);
    bench_lufp2_sub_renorm_no_normalize(r.c[0], r.c[0], t1);
    bench_lufp2_sub_renorm_no_normalize(r.c[0], r.c[0], t2);
    bench_lufp2_mul_by_u_plus_1_renorm_no_normalize(t3, r.c[0]);
    bench_lufp2_add_renorm_no_normalize(r.c[0], t3, t0);

    Fp2Pack8x1x1::add(aa, a.c[0], a.c[1]);
    Fp2Pack8x1x1::add(bb, b.c[0], b.c[1]);
    bench_lufp2_mul_renorm_no_normalize(r.c[1], aa, bb);
    bench_lufp2_sub_renorm_no_normalize(r.c[1], r.c[1], t0);
    bench_lufp2_sub_renorm_no_normalize(r.c[1], r.c[1], t1);
    bench_lufp2_mul_by_u_plus_1_renorm_no_normalize(t3, t2);
    bench_lufp2_add_renorm_no_normalize(r.c[1], r.c[1], t3);

    Fp2Pack8x1x1::add(aa, a.c[0], a.c[2]);
    Fp2Pack8x1x1::add(bb, b.c[0], b.c[2]);
    bench_lufp2_mul_renorm_no_normalize(r.c[2], aa, bb);
    bench_lufp2_sub_renorm_no_normalize(r.c[2], r.c[2], t0);
    bench_lufp2_sub_renorm_no_normalize(r.c[2], r.c[2], t2);
    bench_lufp2_add_renorm_no_normalize(r.c[2], r.c[2], t1);
}

void mul_fp6_pack_lufp_renorm_no_normalize(Fp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    UFp6Pack8x1 wide{};
    mul_fp6_pack_lufp_renorm_no_normalize_wide(wide, a, b);
    UFp6Pack8x1::redc(r, wide);
}

void sqr_fp6_pack_lufp_renorm_no_normalize_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a)
{
    Fp2Pack8x1x1 t0{};
    UFp2Pack8x1x1 s0{}, m01{}, m12{}, s2{}, rx{};

    bench_lufp2_sqr_renorm_no_normalize(s0, a.c[0]);
    bench_lufp2_mul_renorm_no_normalize(m01, a.c[0], a.c[1]);
    bench_lufp2_add_renorm_no_normalize(m01, m01, m01);
    bench_lufp2_mul_renorm_no_normalize(m12, a.c[1], a.c[2]);
    bench_lufp2_add_renorm_no_normalize(m12, m12, m12);
    bench_lufp2_sqr_renorm_no_normalize(s2, a.c[2]);

    Fp2Pack8x1x1::add(t0, a.c[2], a.c[1]);
    Fp2Pack8x1x1::add(t0, t0, a.c[0]);
    bench_lufp2_mul_renorm_no_normalize(rx, t0, t0);
    bench_lufp2_sub_renorm_no_normalize(r.c[2], rx, s0);
    bench_lufp2_sub_renorm_no_normalize(r.c[2], r.c[2], s2);
    bench_lufp2_sub_renorm_no_normalize(r.c[2], r.c[2], m01);
    bench_lufp2_sub_renorm_no_normalize(r.c[2], r.c[2], m12);

    bench_lufp2_mul_by_u_plus_1_renorm_no_normalize(rx, m12);
    bench_lufp2_add_renorm_no_normalize(r.c[0], rx, s0);

    bench_lufp2_mul_by_u_plus_1_renorm_no_normalize(rx, s2);
    bench_lufp2_add_renorm_no_normalize(r.c[1], rx, m01);
}

void sqr_fp6_pack_lufp_renorm_no_normalize(Fp6Pack8x1& r, const Fp6Pack8x1& a)
{
    UFp6Pack8x1 wide{};
    sqr_fp6_pack_lufp_renorm_no_normalize_wide(wide, a);
    UFp6Pack8x1::redc(r, wide);
}

void mul_fp6_pack_baseline_no_repair_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    Fp2Pack8x1x1 aa{}, bb{};
    UFp2Pack8x1x1 t0{}, t1{}, t2{}, t3{};

    bench_ufp2_mul_no_repair(t0, a.c[0], b.c[0]);
    bench_ufp2_mul_no_repair(t1, a.c[1], b.c[1]);
    bench_ufp2_mul_no_repair(t2, a.c[2], b.c[2]);

    Fp2Pack8x1x1::add(aa, a.c[1], a.c[2]);
    Fp2Pack8x1x1::add(bb, b.c[1], b.c[2]);
    bench_ufp2_mul_no_repair(r.c[0], aa, bb);
    bench_ufp2_sub_no_repair(r.c[0], r.c[0], t1);
    bench_ufp2_sub_no_repair(r.c[0], r.c[0], t2);
    bench_ufp2_mul_by_u_plus_1_no_repair(t3, r.c[0]);
    bench_ufp2_add_no_repair(r.c[0], t3, t0);

    Fp2Pack8x1x1::add(aa, a.c[0], a.c[1]);
    Fp2Pack8x1x1::add(bb, b.c[0], b.c[1]);
    bench_ufp2_mul_no_repair(r.c[1], aa, bb);
    bench_ufp2_sub_no_repair(r.c[1], r.c[1], t0);
    bench_ufp2_sub_no_repair(r.c[1], r.c[1], t1);
    bench_ufp2_mul_by_u_plus_1_no_repair(t3, t2);
    bench_ufp2_add_no_repair(r.c[1], r.c[1], t3);

    Fp2Pack8x1x1::add(aa, a.c[0], a.c[2]);
    Fp2Pack8x1x1::add(bb, b.c[0], b.c[2]);
    bench_ufp2_mul_no_repair(r.c[2], aa, bb);
    bench_ufp2_sub_no_repair(r.c[2], r.c[2], t0);
    bench_ufp2_sub_no_repair(r.c[2], r.c[2], t2);
    bench_ufp2_add_no_repair(r.c[2], r.c[2], t1);
}

void mul_fp6_pack_baseline_no_repair(Fp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    UFp6Pack8x1 wide{};
    mul_fp6_pack_baseline_no_repair_wide(wide, a, b);
    UFp6Pack8x1::redc(r, wide);
}

void sqr_fp6_pack_baseline_no_repair_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a)
{
    Fp2Pack8x1x1 t0{};
    UFp2Pack8x1x1 s0{}, m01{}, m12{}, s2{}, rx{};

    bench_ufp2_sqr_no_repair(s0, a.c[0]);
    bench_ufp2_mul_no_repair(m01, a.c[0], a.c[1]);
    bench_ufp2_add_no_repair(m01, m01, m01);
    bench_ufp2_mul_no_repair(m12, a.c[1], a.c[2]);
    bench_ufp2_add_no_repair(m12, m12, m12);
    bench_ufp2_sqr_no_repair(s2, a.c[2]);

    Fp2Pack8x1x1::add(t0, a.c[2], a.c[1]);
    Fp2Pack8x1x1::add(t0, t0, a.c[0]);
    bench_ufp2_mul_no_repair(rx, t0, t0);
    bench_ufp2_sub_no_repair(r.c[2], rx, s0);
    bench_ufp2_sub_no_repair(r.c[2], r.c[2], s2);
    bench_ufp2_sub_no_repair(r.c[2], r.c[2], m01);
    bench_ufp2_sub_no_repair(r.c[2], r.c[2], m12);

    bench_ufp2_mul_by_u_plus_1_no_repair(rx, m12);
    bench_ufp2_add_no_repair(r.c[0], rx, s0);

    bench_ufp2_mul_by_u_plus_1_no_repair(rx, s2);
    bench_ufp2_add_no_repair(r.c[1], rx, m01);
}

void sqr_fp6_pack_baseline_no_repair(Fp6Pack8x1& r, const Fp6Pack8x1& a)
{
    UFp6Pack8x1 wide{};
    sqr_fp6_pack_baseline_no_repair_wide(wide, a);
    UFp6Pack8x1::redc(r, wide);
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

void cyclotomic_sqr_pack_fused(Fp12Pack8x1& r, const Fp12Pack8x1& a)
{
    Fp2Pack8x1x1 t0_0{};
    Fp2Pack8x1x1 t0_1{};
    Fp2Pack8x1x1 tx_0{};
    Fp2Pack8x1x1 tx_1{};
    Fp2Pack8x1x1 tmp{};

    // t0 = sqr_fp4(a00, a11)
    sqr_fp4_pack_bench_v2(t0_0, t0_1, a.c[0].c[0], a.c[1].c[1]);

    // ret00 = 2*(t0_0 - a00) + t0_0
    Fp2Pack8x1x1::sub(r.c[0].c[0], t0_0, a.c[0].c[0]);
    Fp2Pack8x1x1::add(r.c[0].c[0], r.c[0].c[0], r.c[0].c[0]);
    Fp2Pack8x1x1::add(r.c[0].c[0], r.c[0].c[0], t0_0);

    // ret11 = 2*(t0_1 + a11) + t0_1
    Fp2Pack8x1x1::add(r.c[1].c[1], t0_1, a.c[1].c[1]);
    Fp2Pack8x1x1::add(r.c[1].c[1], r.c[1].c[1], r.c[1].c[1]);
    Fp2Pack8x1x1::add(r.c[1].c[1], r.c[1].c[1], t0_1);

    // t1 = sqr_fp4(a10, a02)
    sqr_fp4_pack_bench_v2(tx_0, tx_1, a.c[1].c[0], a.c[0].c[2]);

    // ret01 = 2*(t1_0 - a01) + t1_0
    Fp2Pack8x1x1::sub(r.c[0].c[1], tx_0, a.c[0].c[1]);
    Fp2Pack8x1x1::add(r.c[0].c[1], r.c[0].c[1], r.c[0].c[1]);
    Fp2Pack8x1x1::add(r.c[0].c[1], r.c[0].c[1], tx_0);

    // ret12 = 2*(t1_1 + a12) + t1_1
    Fp2Pack8x1x1::add(r.c[1].c[2], tx_1, a.c[1].c[2]);
    Fp2Pack8x1x1::add(r.c[1].c[2], r.c[1].c[2], r.c[1].c[2]);
    Fp2Pack8x1x1::add(r.c[1].c[2], r.c[1].c[2], tx_1);

    // t2 = sqr_fp4(a01, a12)
    sqr_fp4_pack_bench_v2(tx_0, tx_1, a.c[0].c[1], a.c[1].c[2]);

    // ret02 = 2*(t2_0 - a02) + t2_0
    Fp2Pack8x1x1::sub(r.c[0].c[2], tx_0, a.c[0].c[2]);
    Fp2Pack8x1x1::add(r.c[0].c[2], r.c[0].c[2], r.c[0].c[2]);
    Fp2Pack8x1x1::add(r.c[0].c[2], r.c[0].c[2], tx_0);

    // ret10 = 2*((u+1)*t2_1 + a10) + (u+1)*t2_1
    Fp2Pack8x1x1::mul_by_u_plus_1(tmp, tx_1);
    Fp2Pack8x1x1::add(r.c[1].c[0], tmp, a.c[1].c[0]);
    Fp2Pack8x1x1::add(r.c[1].c[0], r.c[1].c[0], r.c[1].c[0]);
    Fp2Pack8x1x1::add(r.c[1].c[0], r.c[1].c[0], tmp);
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

void validate_fp6_pack_lufp()
{
    static const bool validated = []() {
        Fp6BenchData data(8);
        std::array<Fp6, 8> got{};
        std::array<Fp6, 8> ref{};

        for (std::size_t p = 0; p < data.n_packs; ++p) {
            Fp6Pack8x1 got_mul{};
            Fp6Pack8x1 ref_mul{};
            Fp6Pack8x1 got_sqr{};
            Fp6Pack8x1 ref_sqr{};

            mul_fp6_pack_lufp(got_mul, data.a_pack[p], data.b_pack[p]);
            Fp6Pack8x1::mul(ref_mul, data.a_pack[p], data.b_pack[p]);
            Fp6Pack8x1::to_aos(got.data(), got_mul);
            Fp6Pack8x1::to_aos(ref.data(), ref_mul);
            for (std::size_t lane = 0; lane < got.size(); ++lane) {
                if (!fp6_equal(got[lane], ref[lane])) {
                    std::fprintf(stderr, "Fp6 LUFp mul validation failed at pack %zu lane %zu\n", p, lane);
                    std::abort();
                }
            }

            sqr_fp6_pack_lufp(got_sqr, data.a_pack[p]);
            Fp6Pack8x1::sqr(ref_sqr, data.a_pack[p]);
            Fp6Pack8x1::to_aos(got.data(), got_sqr);
            Fp6Pack8x1::to_aos(ref.data(), ref_sqr);
            for (std::size_t lane = 0; lane < got.size(); ++lane) {
                if (!fp6_equal(got[lane], ref[lane])) {
                    std::fprintf(stderr, "Fp6 LUFp sqr validation failed at pack %zu lane %zu\n", p, lane);
                    std::abort();
                }
            }
        }

        return true;
    }();

    benchmark::DoNotOptimize(validated);
}

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
        , wide_mul_pack_lufpv2(n_packs)
        , wide_mul_pack_lufpv2_fp6_safe(n_packs)
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
            blst::bench_lufp_v2::mul_fp12_pack_lufpv2_wide(wide_mul_pack_lufpv2[p], a_pack[p], b_pack[p]);
            blst::bench_lufp_v2::mul_fp12_pack_lufpv2_fp6_safe_wide(
                wide_mul_pack_lufpv2_fp6_safe[p], a_pack[p], b_pack[p]);
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
    std::vector<UFp12Pack8x1> wide_mul_pack_lufpv2;
    std::vector<UFp12Pack8x1> wide_mul_pack_lufpv2_fp6_safe;
    std::vector<UFp12Pack8x1> wide_sqr_pack;
    std::vector<CyclotomicSqrTerms8x1> cyclo_terms;
    std::vector<Fp12Pack8x1> out_pack;
};

void debug_fp12_pack_lufpv2_fp6_failure(
    const Fp12Pack8x1& a_pack,
    const Fp12Pack8x1& b_pack,
    const Fp12Pack8x1& got_pack,
    const Fp12Pack8x1& ref_pack,
    const std::size_t lane)
{
    UFp6Pack8x1 lufp_t0{};
    UFp6Pack8x1 lufp_t1{};
    UFp6Pack8x1 lufp_m{};
    Fp6Pack8x1 lufp_t0_redc{};
    Fp6Pack8x1 lufp_t1_redc{};
    Fp6Pack8x1 lufp_m_redc{};
    Fp6Pack8x1 ref_t0{};
    Fp6Pack8x1 ref_t1{};
    Fp6Pack8x1 ref_m{};
    Fp12Pack8x1 reduced_products_recompose{};
    Fp6Pack8x1 t3{};
    std::array<Fp12, 8> a{};
    std::array<Fp12, 8> b{};
    std::array<Fp12, 8> got{};
    std::array<Fp12, 8> ref{};
    std::array<Fp6, 8> t0_lufp{};
    std::array<Fp6, 8> t1_lufp{};
    std::array<Fp6, 8> m_lufp{};
    std::array<Fp6, 8> t0_ref{};
    std::array<Fp6, 8> t1_ref{};
    std::array<Fp6, 8> m_ref{};
    std::array<Fp12, 8> reduced_products_recompose_aos{};

    blst::bench_lufp_v2::mul_fp6_pack_lufpv2_wide(lufp_t0, a_pack.c[0], b_pack.c[0]);
    blst::bench_lufp_v2::mul_fp6_pack_lufpv2_wide(lufp_t1, a_pack.c[1], b_pack.c[1]);
    Fp6Pack8x1::add(t3, a_pack.c[0], a_pack.c[1]);
    Fp6Pack8x1 tmp{};
    Fp6Pack8x1::add(tmp, b_pack.c[0], b_pack.c[1]);
    blst::bench_lufp_v2::mul_fp6_pack_lufpv2_wide(lufp_m, t3, tmp);

    UFp6Pack8x1::redc(lufp_t0_redc, lufp_t0);
    UFp6Pack8x1::redc(lufp_t1_redc, lufp_t1);
    UFp6Pack8x1::redc(lufp_m_redc, lufp_m);

    Fp6Pack8x1::mul(ref_t0, a_pack.c[0], b_pack.c[0]);
    Fp6Pack8x1::mul(ref_t1, a_pack.c[1], b_pack.c[1]);
    Fp6Pack8x1::mul(ref_m, t3, tmp);

    reduced_products_recompose.c[1] = lufp_m_redc;
    Fp6Pack8x1::sub(reduced_products_recompose.c[1], reduced_products_recompose.c[1], lufp_t0_redc);
    Fp6Pack8x1::sub(reduced_products_recompose.c[1], reduced_products_recompose.c[1], lufp_t1_redc);
    Fp6Pack8x1::mul_by_nonresidue(reduced_products_recompose.c[0], lufp_t1_redc);
    Fp6Pack8x1::add(reduced_products_recompose.c[0], reduced_products_recompose.c[0], lufp_t0_redc);

    Fp12Pack8x1::to_aos(a.data(), a_pack);
    Fp12Pack8x1::to_aos(b.data(), b_pack);
    Fp12Pack8x1::to_aos(got.data(), got_pack);
    Fp12Pack8x1::to_aos(ref.data(), ref_pack);
    Fp6Pack8x1::to_aos(t0_lufp.data(), lufp_t0_redc);
    Fp6Pack8x1::to_aos(t1_lufp.data(), lufp_t1_redc);
    Fp6Pack8x1::to_aos(m_lufp.data(), lufp_m_redc);
    Fp6Pack8x1::to_aos(t0_ref.data(), ref_t0);
    Fp6Pack8x1::to_aos(t1_ref.data(), ref_t1);
    Fp6Pack8x1::to_aos(m_ref.data(), ref_m);
    Fp12Pack8x1::to_aos(reduced_products_recompose_aos.data(), reduced_products_recompose);

    std::fprintf(stderr, "Fp12 LUFp v2 Fp6 baseline failure diagnostics for lane %zu\n", lane);
    dump_fp12("a", a[lane]);
    dump_fp12("b", b[lane]);
    dump_fp12("got", got[lane]);
    dump_fp12("ref", ref[lane]);
    dump_fp6("t0_lufp_redc", t0_lufp[lane]);
    dump_fp6("t0_ref", t0_ref[lane]);
    dump_fp6("t1_lufp_redc", t1_lufp[lane]);
    dump_fp6("t1_ref", t1_ref[lane]);
    dump_fp6("m_lufp_redc", m_lufp[lane]);
    dump_fp6("m_ref", m_ref[lane]);
    dump_fp12("recompose_from_reduced_lufp_products", reduced_products_recompose_aos[lane]);
}

void validate_fp12_pack_lufpv2_fp6_safe()
{
    static const bool validated = []() {
        Fp12BenchData data(8);
        std::array<Fp12, 8> got{};
        std::array<Fp12, 8> ref{};

        for (std::size_t p = 0; p < data.n_packs; ++p) {
            UFp12Pack8x1 wide;
            Fp12Pack8x1 got_mul;
            blst::bench_lufp_v2::mul_fp12_pack_lufpv2_fp6_safe_wide(wide, data.a_pack[p], data.b_pack[p]);
            UFp12Pack8x1::redc(got_mul, wide);
            Fp12Pack8x1::mul(data.out_pack[p], data.a_pack[p], data.b_pack[p]);

            Fp12Pack8x1::to_aos(got.data(), got_mul);
            Fp12Pack8x1::to_aos(ref.data(), data.out_pack[p]);
            for (std::size_t lane = 0; lane < 8; ++lane) {
                if (!fp12_equal(got[lane], ref[lane])) {
                    std::fprintf(stderr, "Fp12 LUFp v2 Fp6 safe baseline validation failed at pack %zu lane %zu\n", p, lane);
                    debug_fp12_pack_lufpv2_fp6_failure(data.a_pack[p], data.b_pack[p], got_mul, data.out_pack[p], lane);
                    std::fflush(stderr);
                    std::abort();
                }
            }
        }

        return true;
    }();

    benchmark::DoNotOptimize(validated);
}

void validate_fp12_pack_lufpv2()
{
    static const bool validated = []() {
        Fp12BenchData data(8);
        std::array<Fp12, 8> got{};
        std::array<Fp12, 8> ref{};

        for (std::size_t p = 0; p < data.n_packs; ++p) {
            UFp12Pack8x1 wide;
            Fp12Pack8x1 got_mul;
            Fp12Pack8x1 ref_mul;

            blst::bench_lufp_v2::mul_fp12_pack_lufpv2_correct(wide, data.a_pack[p], data.b_pack[p]);
            UFp12Pack8x1::redc(got_mul, wide);
            Fp12Pack8x1::mul(ref_mul, data.a_pack[p], data.b_pack[p]);

            Fp12Pack8x1::to_aos(got.data(), got_mul);
            Fp12Pack8x1::to_aos(ref.data(), ref_mul);
            for (std::size_t lane = 0; lane < 8; ++lane) {
                if (!fp12_equal(got[lane], ref[lane])) {
                    std::fprintf(stderr, "Fp12 LUFp v2 validation failed at pack %zu lane %zu\n", p, lane);
                    dump_fp12("got", got[lane]);
                    dump_fp12("ref", ref[lane]);
                    std::fflush(stderr);
                    std::abort();
                }
            }
        }

        return true;
    }();

    benchmark::DoNotOptimize(validated);
}

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

// remove? legacy Fp6 sqr baseline remains here because the new accounting suite
// does not cover sqr yet.
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

// remove? legacy LUFp sqr experiment; still unique until sqr enters the new suite.
void BM_Fp6Sqr_Pack8x1_LUFp(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp6BenchData data(packs);
    validate_fp6_pack_lufp();

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            sqr_fp6_pack_lufp(data.out_pack[p], data.a_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp6Sqr_Pack8x1_LUFp)->Arg(256);

// remove? legacy no-repair sqr ceiling experiment; still unique until sqr enters the new suite.
void BM_Fp6Sqr_Pack8x1_BaselineNoRepair(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp6BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            sqr_fp6_pack_baseline_no_repair(data.out_pack[p], data.a_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp6Sqr_Pack8x1_BaselineNoRepair)->Arg(256);

// remove? legacy no-normalize sqr experiment; still unique until sqr enters the new suite.
void BM_Fp6Sqr_Pack8x1_LUFp_NoNormalize(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp6BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            sqr_fp6_pack_lufp_no_normalize(data.out_pack[p], data.a_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp6Sqr_Pack8x1_LUFp_NoNormalize)->Arg(256);

// remove? legacy renorm-no-normalize sqr experiment; still unique until sqr enters the new suite.
void BM_Fp6Sqr_Pack8x1_LUFp_RenormNoNormalize(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp6BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            sqr_fp6_pack_lufp_renorm_no_normalize(data.out_pack[p], data.a_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp6Sqr_Pack8x1_LUFp_RenormNoNormalize)->Arg(256);

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

void BM_Fp12Mul_Pack8x1_LUFpV2(benchmark::State& state)
{
    validate_fp12_pack_lufpv2();

    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            UFp12Pack8x1 wide;
            blst::bench_lufp_v2::mul_fp12_pack_lufpv2_wide(wide, data.a_pack[p], data.b_pack[p]);
            UFp12Pack8x1::redc(data.out_pack[p], wide);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Mul_Pack8x1_LUFpV2)->Arg(128);

void BM_Fp12Mul_Pack8x1_LUFpV2Fp6Safe(benchmark::State& state)
{
    validate_fp12_pack_lufpv2_fp6_safe();

    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            UFp12Pack8x1 wide;
            blst::bench_lufp_v2::mul_fp12_pack_lufpv2_fp6_safe_wide(wide, data.a_pack[p], data.b_pack[p]);
            UFp12Pack8x1::redc(data.out_pack[p], wide);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Mul_Pack8x1_LUFpV2Fp6Safe)->Arg(128);

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

void BM_Fp12Mul_Pack8x1_LUFpV2_WideOnly(benchmark::State& state)
{
    validate_fp12_pack_lufpv2();

    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            blst::bench_lufp_v2::mul_fp12_pack_lufpv2_wide(data.wide_mul_pack_lufpv2[p], data.a_pack[p], data.b_pack[p]);
        }
        benchmark::DoNotOptimize(data.wide_mul_pack_lufpv2.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Mul_Pack8x1_LUFpV2_WideOnly)->Arg(128);

void BM_Fp12Mul_Pack8x1_LUFpV2Fp6Safe_WideOnly(benchmark::State& state)
{
    validate_fp12_pack_lufpv2_fp6_safe();

    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            blst::bench_lufp_v2::mul_fp12_pack_lufpv2_fp6_safe_wide(
                data.wide_mul_pack_lufpv2_fp6_safe[p], data.a_pack[p], data.b_pack[p]);
        }
        benchmark::DoNotOptimize(data.wide_mul_pack_lufpv2_fp6_safe.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Mul_Pack8x1_LUFpV2Fp6Safe_WideOnly)->Arg(128);

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

void BM_Fp12Mul_Pack8x1_LUFpV2_RedcOnly(benchmark::State& state)
{
    validate_fp12_pack_lufpv2();

    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            UFp12Pack8x1::redc(data.out_pack[p], data.wide_mul_pack_lufpv2[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Mul_Pack8x1_LUFpV2_RedcOnly)->Arg(128);

void BM_Fp12Mul_Pack8x1_LUFpV2Fp6Safe_RedcOnly(benchmark::State& state)
{
    validate_fp12_pack_lufpv2_fp6_safe();

    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            UFp12Pack8x1::redc(data.out_pack[p], data.wide_mul_pack_lufpv2_fp6_safe[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12Mul_Pack8x1_LUFpV2Fp6Safe_RedcOnly)->Arg(128);

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

void BM_Fp12CyclotomicSqr_Pack8x1_Fused(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    Fp12BenchData data(packs);

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            cyclotomic_sqr_pack_fused(data.out_pack[p], data.cyclo_pack[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_Fp12CyclotomicSqr_Pack8x1_Fused)->Arg(128);

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

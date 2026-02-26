#include "avxbls_bench.h"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "blst/common/defines.h"

#if defined(_MSC_VER) || defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include "blst/curve/G2_avx.h"
#include "blst/curve/GT.h"
#include "blst/curve/pairing_avx.h"
#include "blst/curve/point.h"
#include "blst/fields/Fp12_avx.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fp2_avx.h"
#include "blst/fields/Fp4_avx.h"
#include "blst/fields/Fp6_avx.h"
#include "blst/fields/Fp_avx.h"

using namespace blst;

namespace
{
constexpr u64 kFpWarmup = 10000;
constexpr u64 kFpIters = 100000;
constexpr u64 kPairingWarmup = 1000;
constexpr u64 kPairingIters = 10000;

constexpr std::array<u64, SWORDS> TV_A = {
    0x0123456789ABCDEFULL, 0x0123456789ABCDEFULL, 0x0123456789ABCDEFULL,
    0x0123456789ABCDEFULL, 0x0123456789ABCDEFULL, 0x1123456789ABCDEFULL,
};

constexpr std::array<u64, SWORDS> TV_B = {
    0xABCDEF0123456789ULL, 0xABCDEF0123456789ULL, 0xABCDEF0123456789ULL,
    0xABCDEF0123456789ULL, 0xABCDEF0123456789ULL, 0x19CDEF0123456789ULL,
};

struct BenchConfig
{
    u64 fp_warmup = kFpWarmup;
    u64 fp_iters = kFpIters;
    u64 pairing_warmup = kPairingWarmup;
    u64 pairing_iters = kPairingIters;
};

struct LineEval
{
    Fp2 c0;
    Fp2 c1;
    Fp2 c2;
};

struct JacobianPoint
{
    Fp2 X;
    Fp2 Y;
    Fp2 Z;
};

struct AffineInputPoint
{
    Fp2 X;
    Fp2 Y;
};

inline u64 read_tsc_cycles()
{
#if defined(_MSC_VER) || defined(__x86_64__) || defined(__i386__)
    _mm_lfence();
    const u64 t = __rdtsc();
    _mm_lfence();
    return t;
#else
    return 0;
#endif
}

template <class Fn>
u64 measure_cycles(const char* label, u64 warmup, u64 iters, Fn&& fn)
{
    for (u64 i = 0; i < warmup; ++i) {
        fn();
    }

    const u64 begin = read_tsc_cycles();
    for (u64 i = 0; i < iters; ++i) {
        fn();
    }
    const u64 end = read_tsc_cycles();
    const u64 per_iter = iters == 0 ? 0 : ((end - begin) / iters);

    std::printf("%s#cycle = %llu\n", label, static_cast<unsigned long long>(per_iter));
    return per_iter;
}

bool parse_u64(const char* text, u64* out)
{
    if (text == nullptr || *text == '\0' || out == nullptr) {
        return false;
    }
    char* end = nullptr;
    const unsigned long long value = std::strtoull(text, &end, 10);
    if (*end != '\0') {
        return false;
    }
    *out = static_cast<u64>(value);
    return true;
}

Fp make_fp_from_words(const std::array<u64, SWORDS>& words)
{
    Fp out{};
    for (u64 i = 0; i < SWORDS; ++i) {
        out[i] = words[i];
    }
    return out;
}

void fill_scalar_fp(Fp& out, u64 seed)
{
    for (u64 i = 0; i < SWORDS; ++i) {
        out[i] = (seed + i + 1) & BMASK;
    }
}

void fill_scalar_fp2(Fp2& out, u64 seed)
{
    fill_scalar_fp(out[0], seed + 0);
    fill_scalar_fp(out[1], seed + 17);
}

void fill_scalar_fp6(Fp6& out, u64 seed)
{
    fill_scalar_fp2(out[0], seed + 0);
    fill_scalar_fp2(out[1], seed + 41);
    fill_scalar_fp2(out[2], seed + 83);
}

void fill_scalar_fp12(Fp12& out, u64 seed)
{
    fill_scalar_fp6(out[0], seed + 0);
    fill_scalar_fp6(out[1], seed + 121);
}

void fill_vec_words(cc512i* out, u64 count, u64 seed)
{
    for (u64 i = 0; i < count; ++i) {
        out[i] = VSET1(static_cast<i64>((seed + i) & BMASK));
    }
}

template<class PackT>
void fill_vec_words(PackT& out, u64 count, u64 seed)
{
    fill_vec_words(out.data(), count, seed);
}

void init_fp_inputs(FpPack8x1& a_8x1w, FpPack8x1& b_8x1w, FpPack4x2& a_4x2w, FpPack4x2& b_4x2w)
{
    std::array<u64, NWORDS> a48 = {};
    std::array<u64, NWORDS> b48 = {};
    conv_64to48_mpi(a48.data(), TV_A.data(), static_cast<int>(NWORDS), static_cast<int>(SWORDS));
    conv_64to48_mpi(b48.data(), TV_B.data(), static_cast<int>(NWORDS), static_cast<int>(SWORDS));

    for (u64 i = 0; i < NWORDS; ++i) {
        a_8x1w[i] = VSET1(static_cast<i64>(a48[i]));
        b_8x1w[i] = VSET1(static_cast<i64>(b48[i]));
    }

    for (u64 i = 0; i < VWORDS; ++i) {
        a_4x2w[i] = VSET(0, 0, 0, 0, 0, 0, static_cast<i64>(a48[i + VWORDS]), static_cast<i64>(a48[i]));
        b_4x2w[i] = VSET(0, 0, 0, 0, 0, 0, static_cast<i64>(b48[i + VWORDS]), static_cast<i64>(b48[i]));
    }
}

void init_fp2_inputs(
    Fp2Pack8x1x1& a_8x1x1w,
    Fp2Pack8x1x1& b_8x1x1w,
    Fp2Pack4x2x1& a_4x2x1w,
    Fp2Pack4x2x1& b_4x2x1w,
    Fp2Pack2x2x2& a_2x2x2w,
    Fp2Pack2x2x2& b_2x2x2w,
    Fp2Pack2x4x1& a_2x4x1w,
    Fp2Pack2x4x1& b_2x4x1w,
    Fp2Pack1x4x2& a_1x4x2w,
    Fp2Pack1x4x2& b_1x4x2w)
{
    std::array<u64, NWORDS> a48 = {};
    std::array<u64, NWORDS> b48 = {};
    conv_64to48_mpi(a48.data(), TV_A.data(), static_cast<int>(NWORDS), static_cast<int>(SWORDS));
    conv_64to48_mpi(b48.data(), TV_B.data(), static_cast<int>(NWORDS), static_cast<int>(SWORDS));

    for (u64 i = 0; i < NWORDS; ++i) {
        a_4x2x1w[i] = VSET1(static_cast<i64>(a48[i]));
        b_4x2x1w[i] = VSET1(static_cast<i64>(b48[i]));
        a_2x4x1w[i] = VSET1(static_cast<i64>((a48[i] + 17) & BMASK));
        b_2x4x1w[i] = VSET1(static_cast<i64>((b48[i] + 19) & BMASK));
        a_8x1x1w[0][i] = VSET1(static_cast<i64>((a48[i] + 5) & BMASK));
        a_8x1x1w[1][i] = VSET1(static_cast<i64>((b48[i] + 7) & BMASK));
        b_8x1x1w[0][i] = VSET1(static_cast<i64>((b48[i] + 11) & BMASK));
        b_8x1x1w[1][i] = VSET1(static_cast<i64>((a48[i] + 13) & BMASK));
    }

    for (u64 i = 0; i < VWORDS; ++i) {
        a_2x2x2w[i] = VSET(
            0,
            0,
            static_cast<i64>(a48[i + VWORDS]),
            static_cast<i64>(a48[i]),
            0,
            0,
            static_cast<i64>(a48[i + VWORDS]),
            static_cast<i64>(a48[i]));
        b_2x2x2w[i] = VSET(
            0,
            0,
            static_cast<i64>(b48[i + VWORDS]),
            static_cast<i64>(b48[i]),
            0,
            0,
            static_cast<i64>(b48[i + VWORDS]),
            static_cast<i64>(b48[i]));
        a_1x4x2w[i] = VSET1(static_cast<i64>((a48[i] + 23) & BMASK));
        b_1x4x2w[i] = VSET1(static_cast<i64>((b48[i] + 29) & BMASK));
    }
}

void sqr_fp4_scalar(Fp4& ret, const Fp2& a0, const Fp2& a1)
{
    Fp2 t0{};
    Fp2 t1{};
    Fp2 t2{};

    Fp2::sqr(t0, a0);
    Fp2::sqr(t1, a1);
    Fp2::mul(t2, a0, a1);

    Fp2::mul_by_u_plus_1(t1, t1);
    Fp2::add(ret[0], t0, t1);
    Fp2::add(ret[1], t2, t2);
}

#if defined(BLST_ASM)
struct UFp6Bench
{
    UFp2 data[3];

    UFp2& operator[](std::size_t idx) { return data[idx]; }
    const UFp2& operator[](std::size_t idx) const { return data[idx]; }
};

void add_fp2x2_scalar(UFp2& ret, const UFp2& a, const UFp2& b)
{
    add_mod_384x384(&ret[0], &a[0], &b[0], &get_BLST_P());
    add_mod_384x384(&ret[1], &a[1], &b[1], &get_BLST_P());
}

void sub_fp2x2_scalar(UFp2& ret, const UFp2& a, const UFp2& b)
{
    sub_mod_384x384(&ret[0], &a[0], &b[0], &get_BLST_P());
    sub_mod_384x384(&ret[1], &a[1], &b[1], &get_BLST_P());
}

void mul_by_u_plus_1_fp2x2_scalar(UFp2& ret, const UFp2& a)
{
    sub_mod_384x384(&ret[0], &a[0], &a[1], &get_BLST_P());
    add_mod_384x384(&ret[1], &a[0], &a[1], &get_BLST_P());
}

void redc_fp2x2_scalar(Fp2& ret, const UFp2& a)
{
    Fp::redc(ret[0], a[0]);
    Fp::redc(ret[1], a[1]);
}

void mul_fp2x2_scalar(UFp2& ret, const Fp2& a, const Fp2& b)
{
    mul_v382x(ret, a, b, get_BLST_P());
}

void sqr_fp2x2_scalar(UFp2& ret, const Fp2& a)
{
    sqr_v382x(ret, a, get_BLST_P());
}

void mul_fp6x2_scalar(UFp6Bench& ret, const Fp6& a, const Fp6& b)
{
    UFp2 t0{};
    UFp2 t1{};
    UFp2 t2{};
    Fp2 aa{};
    Fp2 bb{};
    const Fp2& a0 = a[0];

    mul_fp2x2_scalar(t0, a0, b[0]);
    mul_fp2x2_scalar(t1, a[1], b[1]);
    mul_fp2x2_scalar(t2, a[2], b[2]);

    Fp2::add(aa, a[1], a[2]);
    Fp2::add(bb, b[1], b[2]);
    mul_fp2x2_scalar(ret[0], aa, bb);
    sub_fp2x2_scalar(ret[0], ret[0], t1);
    sub_fp2x2_scalar(ret[0], ret[0], t2);
    mul_by_u_plus_1_fp2x2_scalar(ret[1], ret[0]);
    add_fp2x2_scalar(ret[0], ret[1], t0);

    Fp2::add(aa, a[0], a[1]);
    Fp2::add(bb, b[0], b[1]);
    mul_fp2x2_scalar(ret[1], aa, bb);
    sub_fp2x2_scalar(ret[1], ret[1], t0);
    sub_fp2x2_scalar(ret[1], ret[1], t1);
    mul_by_u_plus_1_fp2x2_scalar(ret[2], t2);
    add_fp2x2_scalar(ret[1], ret[1], ret[2]);

    Fp2::add(aa, a[0], a[2]);
    Fp2::add(bb, b[0], b[2]);
    mul_fp2x2_scalar(ret[2], aa, bb);
    sub_fp2x2_scalar(ret[2], ret[2], t0);
    sub_fp2x2_scalar(ret[2], ret[2], t2);
    add_fp2x2_scalar(ret[2], ret[2], t1);
}

void redc_fp6x2_scalar(Fp6& ret, const UFp6Bench& a)
{
    redc_fp2x2_scalar(ret[0], a[0]);
    redc_fp2x2_scalar(ret[1], a[1]);
    redc_fp2x2_scalar(ret[2], a[2]);
}
#endif

void mul_fp6_scalar(Fp6& ret, const Fp6& a, const Fp6& b)
{
#if defined(BLST_ASM)
    UFp6Bench wide{};
    mul_fp6x2_scalar(wide, a, b);
    redc_fp6x2_scalar(ret, wide);
#else
    Fp2 a0b0{};
    Fp2 a1b2{};
    Fp2 a2b1{};
    Fp2 a0b1{};
    Fp2 a1b0{};
    Fp2 a2b2{};
    Fp2 a0b2{};
    Fp2 a2b0{};
    Fp2 a1b1{};
    Fp2 t{};

    Fp2::mul(a0b0, a[0], b[0]);
    Fp2::mul(a1b2, a[1], b[2]);
    Fp2::mul(a2b1, a[2], b[1]);
    Fp2::add(t, a1b2, a2b1);
    Fp2::mul_by_u_plus_1(t, t);
    Fp2::add(ret[0], a0b0, t);

    Fp2::mul(a0b1, a[0], b[1]);
    Fp2::mul(a1b0, a[1], b[0]);
    Fp2::mul(a2b2, a[2], b[2]);
    Fp2::mul_by_u_plus_1(a2b2, a2b2);
    Fp2::add(t, a0b1, a1b0);
    Fp2::add(ret[1], t, a2b2);

    Fp2::mul(a0b2, a[0], b[2]);
    Fp2::mul(a2b0, a[2], b[0]);
    Fp2::mul(a1b1, a[1], b[1]);
    Fp2::add(t, a0b2, a2b0);
    Fp2::add(ret[2], t, a1b1);
#endif
}

void line_dbl_scalar(LineEval& line, JacobianPoint& T, const JacobianPoint& Q)
{
    Fp2 ZZ{};
    Fp2 A{};
    Fp2 B{};
    Fp2 C{};
    Fp2 D{};
    Fp2 E{};
    Fp2 F{};

    Fp2::sqr(A, Q.X);
    Fp2::sqr(B, Q.Y);
    Fp2::sqr(ZZ, Q.Z);
    Fp2::sqr(C, B);

    Fp2::add(D, Q.X, B);
    Fp2::sqr(D, D);
    Fp2::sub(D, D, A);
    Fp2::sub(D, D, C);
    Fp2::add(D, D, D);

    Fp2::mul_by_3(E, A);
    Fp2::sqr(F, E);

    Fp2::add(line.c0, E, Q.X);

    Fp2::sub(T.X, F, D);
    Fp2::sub(T.X, T.X, D);

    Fp2::add(T.Z, Q.Y, Q.Z);
    Fp2::sqr(T.Z, T.Z);
    Fp2::sub(T.Z, T.Z, B);
    Fp2::sub(T.Z, T.Z, ZZ);

    Fp2::mul_by_8(C, C);
    Fp2::sub(T.Y, D, T.X);
    Fp2::mul(T.Y, T.Y, E);
    Fp2::sub(T.Y, T.Y, C);

    Fp2::sqr(line.c0, line.c0);
    Fp2::sub(line.c0, line.c0, A);
    Fp2::sub(line.c0, line.c0, F);
    Fp2::lshift(B, B, 2);
    Fp2::sub(line.c0, line.c0, B);

    Fp2::mul(line.c1, E, ZZ);
    Fp2::mul(line.c2, T.Z, ZZ);
}

void line_add_scalar(LineEval& line, JacobianPoint& T, const JacobianPoint& R, const AffineInputPoint& Q)
{
    Fp2 Z1Z1{};
    Fp2 U2{};
    Fp2 S2{};
    Fp2 H{};
    Fp2 HH{};
    Fp2 I{};
    Fp2 J{};
    Fp2 V{};
    Fp2 r{};

    Fp2::sqr(Z1Z1, R.Z);
    Fp2::mul(U2, Q.X, Z1Z1);

    Fp2::mul(S2, Q.Y, R.Z);
    Fp2::mul(S2, S2, Z1Z1);

    Fp2::sub(H, U2, R.X);

    Fp2::sqr(HH, H);
    Fp2::add(I, HH, HH);
    Fp2::add(I, I, I);

    Fp2::mul(J, H, I);

    Fp2::sub(r, S2, R.Y);
    Fp2::add(r, r, r);

    Fp2::mul(V, R.X, I);

    Fp2::sqr(T.X, r);
    Fp2::sub(T.X, T.X, J);
    Fp2::sub(T.X, T.X, V);
    Fp2::sub(T.X, T.X, V);

    Fp2::mul(J, J, R.Y);
    Fp2::sub(T.Y, V, T.X);
    Fp2::mul(T.Y, T.Y, r);
    Fp2::sub(T.Y, T.Y, J);
    Fp2::sub(T.Y, T.Y, J);

    Fp2::add(T.Z, R.Z, H);
    Fp2::sqr(T.Z, T.Z);
    Fp2::sub(T.Z, T.Z, Z1Z1);
    Fp2::sub(T.Z, T.Z, HH);

    Fp2::mul(I, r, Q.X);
    Fp2::mul(J, Q.Y, T.Z);
    Fp2::sub(I, I, J);
    Fp2::add(line.c0, I, I);
    line.c1 = r;
    line.c2 = T.Z;
}

void init_line_inputs(
    JacobianPoint& q_scalar,
    AffineInputPoint& r_scalar,
    Fp2Pack2x2x2& X1Y1,
    Fp2Pack2x2x2& Z1,
    Fp2Pack2x2x2& Z1Y2,
    Fp2Pack2x2x2& X2,
    Fp2Pack4x2x1& X1Y1Z1)
{
    q_scalar.X[0][0] = 1;
    q_scalar.X[1][0] = 2;
    q_scalar.Y[0][0] = 3;
    q_scalar.Y[1][0] = 4;
    q_scalar.Z[0][0] = 5;
    q_scalar.Z[1][0] = 6;

    r_scalar.X[0][0] = 7;
    r_scalar.X[1][0] = 8;
    r_scalar.Y[0][0] = 9;
    r_scalar.Y[1][0] = 10;

    for (u64 i = 1; i < SWORDS; ++i) {
        q_scalar.X[0][i] = 0;
        q_scalar.X[1][i] = 0;
        q_scalar.Y[0][i] = 0;
        q_scalar.Y[1][i] = 0;
        q_scalar.Z[0][i] = 0;
        q_scalar.Z[1][i] = 0;

        r_scalar.X[0][i] = 0;
        r_scalar.X[1][i] = 0;
        r_scalar.Y[0][i] = 0;
        r_scalar.Y[1][i] = 0;
    }

    X1Y1[0] = VSET(0, 4, 0, 3, 0, 2, 0, 1);
    Z1[0] = VSET(0, 6, 0, 5, 0, 6, 0, 5);
    Z1Y2[0] = VSET(0, 6, 0, 5, 0, 10, 0, 9);
    X2[0] = VSET(0, 8, 0, 7, 0, 8, 0, 7);
    X1Y1Z1[0] = VSET(4, 3, 2, 1, 6, 5, 4, 3);

    for (u64 i = 1; i < VWORDS; ++i) {
        X1Y1[i] = VZERO();
        Z1[i] = VZERO();
        Z1Y2[i] = VZERO();
        X2[i] = VZERO();
    }

    for (u64 i = 1; i < NWORDS; ++i) {
        X1Y1Z1[i] = VZERO();
    }
}

void bench_pairing(const BenchConfig& cfg)
{
    std::puts("\n=============================================================\n");
    std::puts("TIMING - PAIRING\n");

    G1AffinePoint p{};
    G2AffinePoint q{};
    G2Point q_j = get_BLST_G2();
    G2Point::doubling(q_j, q_j);
    G2Point::to_affine(q, q_j);
    G1Point::to_affine(p, get_BLST_G1());

    GT f{};
    GT e1{};

    measure_cycles("- miller_loop:        ", cfg.pairing_warmup, cfg.pairing_iters, [&]() {
        avx::miller_loop_n(f, &q, &p, 1);
    });

    avx::miller_loop_n(f, &q, &p, 1);
    measure_cycles("- final_exp:          ", cfg.pairing_warmup, cfg.pairing_iters, [&]() {
        avx::final_exp(e1, f);
    });

    measure_cycles("- optimal_ate_pairing:", cfg.pairing_warmup, cfg.pairing_iters, [&]() {
        avx::optimal_ate_pairing(e1, &q, &p, 1);
    });
}

void bench_fp(const BenchConfig& cfg)
{
    FpPack8x1 a_8x1w{};
    FpPack8x1 b_8x1w{};
    FpPack8x1 r_8x1w{};
    UFpPack8x1 z_8x1w{};
    FpPack4x2 a_4x2w{};
    FpPack4x2 b_4x2w{};
    FpPack4x2 r_4x2w{};
    UFpPack4x2 z_4x2w{};
    init_fp_inputs(a_8x1w, b_8x1w, a_4x2w, b_4x2w);

    std::puts("\n=============================================================\n");
    std::puts("TIMING - FP\n");

    measure_cycles("- FpPack8x1::mul_v2:     ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        FpPack8x1::mul_v2(z_8x1w, a_8x1w, b_8x1w);
    });

    measure_cycles("- FpPack4x2::mul_v1:     ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        FpPack4x2::mul_v1(z_4x2w, a_4x2w, b_4x2w);
    });

    measure_cycles("- UFpPack8x1::redc:       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        UFpPack8x1::redc(r_8x1w, z_8x1w);
    });

    measure_cycles("- UFpPack4x2::redc:       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        UFpPack4x2::redc(r_4x2w, z_4x2w);
    });
}

void bench_fp2(const BenchConfig& cfg)
{
    Fp2Pack8x1x1 a_8x1x1w{};
    Fp2Pack8x1x1 b_8x1x1w{};
    UFp2Pack8x1x1 z_8x1x1w{};
    Fp2Pack4x2x1 a_4x2x1w{};
    Fp2Pack4x2x1 b_4x2x1w{};
    Fp2Pack4x2x1 r_4x2x1w{};
    UFp2Pack4x2x1 z_4x2x1w{};
    Fp2Pack2x2x2 a_2x2x2w{};
    Fp2Pack2x2x2 b_2x2x2w{};
    Fp2Pack2x2x2 r_2x2x2w{};
    UFp2Pack2x2x2 z_2x2x2w{};
    Fp2Pack2x4x1 a_2x4x1w{};
    Fp2Pack2x4x1 b_2x4x1w{};
    UFp2Pack2x4x1 z_2x4x1w{};
    Fp2Pack1x4x2 a_1x4x2w{};
    Fp2Pack1x4x2 b_1x4x2w{};
    UFp2Pack1x4x2 z_1x4x2w{};
    Fp2 a_scalar{};
    Fp2 b_scalar{};
    Fp2 r_scalar{};
#if defined(BLST_ASM)
    UFp2 z_scalar{};
#endif

    init_fp2_inputs(
        a_8x1x1w,
        b_8x1x1w,
        a_4x2x1w,
        b_4x2x1w,
        a_2x2x2w,
        b_2x2x2w,
        a_2x4x1w,
        b_2x4x1w,
        a_1x4x2w,
        b_1x4x2w);
    a_scalar[0] = make_fp_from_words(TV_A);
    a_scalar[1] = make_fp_from_words(TV_B);
    b_scalar[0] = make_fp_from_words(TV_B);
    b_scalar[1] = make_fp_from_words(TV_A);

    std::puts("\n=============================================================\n");
    std::puts("TIMING - FP2\n");

    measure_cycles("- mul_fp2_scalar:       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        Fp2::mul(r_scalar, a_scalar, b_scalar);
    });

    measure_cycles("- Fp2Pack4x2x1::mul:       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        Fp2Pack4x2x1::mul(r_4x2x1w, a_4x2x1w, b_4x2x1w);
    });

    measure_cycles("- mul_Fp2Pack2x2x2:       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        Fp2Pack2x2x2::mul(r_2x2x2w, a_2x2x2w, b_2x2x2w);
    });

#if defined(BLST_ASM)
    measure_cycles("- mul_fp2x2_scalar:     ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_fp2x2_scalar(z_scalar, a_scalar, b_scalar);
    });
#endif

    measure_cycles("- UFp2Pack8x1x1::mul:     ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        UFp2Pack8x1x1::mul(z_8x1x1w, a_8x1x1w, b_8x1x1w);
    });

    measure_cycles("- UFp2Pack4x2x1::mul:     ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        UFp2Pack4x2x1::mul(z_4x2x1w, a_4x2x1w, b_4x2x1w);
    });

    measure_cycles("- mul_UFp2Pack2x4x1:     ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        UFp2Pack2x4x1::mul(z_2x4x1w, a_2x4x1w, b_2x4x1w);
    });

    measure_cycles("- mul_UFp2Pack2x2x2:     ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        UFp2Pack2x2x2::mul(z_2x2x2w, a_2x2x2w, b_2x2x2w);
    });

    measure_cycles("- mul_UFp2Pack1x4x2:     ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        UFp2Pack1x4x2::mul(z_1x4x2w, a_1x4x2w, b_1x4x2w);
    });

    measure_cycles("- sqr_fp2_scalar:       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        Fp2::sqr(r_scalar, a_scalar);
    });

    measure_cycles("- Fp2Pack4x2x1::sqr:       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        Fp2Pack4x2x1::sqr(r_4x2x1w, a_4x2x1w);
    });

    measure_cycles("- sqr_Fp2Pack2x2x2:       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        Fp2Pack2x2x2::sqr(r_2x2x2w, a_2x2x2w);
    });

#if defined(BLST_ASM)
    measure_cycles("- sqr_fp2x2_scalar:     ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        sqr_fp2x2_scalar(z_scalar, a_scalar);
    });
#endif

    measure_cycles("- UFp2Pack4x2x1::sqr:     ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        UFp2Pack4x2x1::sqr(z_4x2x1w, a_4x2x1w);
    });

    measure_cycles("- sqr_UFp2Pack2x2x2:     ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        UFp2Pack2x2x2::sqr(z_2x2x2w, a_2x2x2w);
    });
}

void bench_fp4(const BenchConfig& cfg)
{
    Fp4Pack2x2x2x1 a_2x2x2x1w{};
    Fp4Pack2x2x2x1 r_2x2x2x1w{};
    Fp4Pack1x2x2x2 a_1x2x2x2w{};
    Fp4Pack1x2x2x2 r_1x2x2x2w{};
    Fp4 r_scalar{};
    Fp2 a0_scalar{};
    Fp2 a1_scalar{};

    fill_vec_words(a_2x2x2x1w, NWORDS, 101);
    fill_vec_words(a_1x2x2x2w, VWORDS, 211);
    fill_scalar_fp2(a0_scalar, 7);
    fill_scalar_fp2(a1_scalar, 37);

    std::puts("\n=============================================================\n");
    std::puts("TIMING - FP4\n");

    measure_cycles("- sqr_fp4_scalar:       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        sqr_fp4_scalar(r_scalar, a0_scalar, a1_scalar);
    });

    measure_cycles("- Fp4Pack2x2x2x1::sqr_v1:  ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        Fp4Pack2x2x2x1::sqr_v1(r_2x2x2x1w, a_2x2x2x1w);
    });

    measure_cycles("- Fp4Pack1x2x2x2::sqr_v1:  ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        Fp4Pack1x2x2x2::sqr_v1(r_1x2x2x2w, a_1x2x2x2w);
    });
}

void bench_fp6(const BenchConfig& cfg)
{
    Fp2Pack8x1x1 a0_8x1x1w{};
    Fp2Pack8x1x1 a1_8x1x1w{};
    Fp2Pack8x1x1 a2_8x1x1w{};
    Fp2Pack8x1x1 a_8x1x1w{};
    Fp2Pack8x1x1 b_8x1x1w{};
    Fp2Pack4x2x1 a0_4x2x1w{};
    Fp2Pack4x2x1 a1_4x2x1w{};
    Fp2Pack4x2x1 a2_4x2x1w{};
    Fp2Pack4x2x1 a01_4x2x1w{};
    Fp2Pack4x2x1 a2_pair_4x2x1w{};
    Fp2Pack4x2x1 b01_4x2x1w{};
    Fp2Pack4x2x1 b4_4x2x1w{};
    Fp2Pack2x2x2 a0_2x2x2w{};
    Fp2Pack2x2x2 a1_2x2x2w{};
    Fp2Pack2x2x2 a2_2x2x2w{};

    UFp2Pack8x1x1 r01_8x1x1w{};
    UFp2Pack8x1x1 r23_8x1x1w{};
    UFp2Pack8x1x1 r45_8x1x1w{};
    UFp2Pack4x2x1 z01_4x2x1w{};
    UFp2Pack4x2x1 z2_4x2x1w{};
    UFp2Pack4x2x1 z3_4x2x1w{};
    UFp2Pack4x2x1 z45_4x2x1w{};
    UFp2Pack2x2x2 z01_2x2x2w{};
    UFp2Pack2x2x2 z2_2x2x2w{};

    Fp6 a_scalar{};
    Fp6 b_scalar{};
    Fp6 r_scalar{};

    fill_vec_words(a0_8x1x1w[0], NWORDS, 301);
    fill_vec_words(a0_8x1x1w[1], NWORDS, 331);
    fill_vec_words(a1_8x1x1w[0], NWORDS, 361);
    fill_vec_words(a1_8x1x1w[1], NWORDS, 391);
    fill_vec_words(a2_8x1x1w[0], NWORDS, 421);
    fill_vec_words(a2_8x1x1w[1], NWORDS, 451);
    fill_vec_words(a_8x1x1w[0], NWORDS, 481);
    fill_vec_words(a_8x1x1w[1], NWORDS, 511);
    fill_vec_words(b_8x1x1w[0], NWORDS, 541);
    fill_vec_words(b_8x1x1w[1], NWORDS, 571);

    fill_vec_words(a0_4x2x1w, NWORDS, 601);
    fill_vec_words(a1_4x2x1w, NWORDS, 631);
    fill_vec_words(a2_4x2x1w, NWORDS, 661);
    fill_vec_words(a01_4x2x1w, NWORDS, 691);
    fill_vec_words(a2_pair_4x2x1w, NWORDS, 721);
    fill_vec_words(b01_4x2x1w, NWORDS, 751);
    fill_vec_words(b4_4x2x1w, NWORDS, 781);

    fill_vec_words(a0_2x2x2w, VWORDS, 811);
    fill_vec_words(a1_2x2x2w, VWORDS, 841);
    fill_vec_words(a2_2x2x2w, VWORDS, 871);

    fill_scalar_fp6(a_scalar, 17);
    fill_scalar_fp6(b_scalar, 117);

    std::puts("\n=============================================================\n");
    std::puts("TIMING - FP6\n");

    measure_cycles("- mul_fp6x2_scalar:     ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_fp6_scalar(r_scalar, a_scalar, b_scalar);
    });

    measure_cycles("- mul_fp6x2_2x2x2x1w:   ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_fp6x2_2x2x2x1w(z01_4x2x1w, z2_4x2x1w, a0_4x2x1w, a1_4x2x1w, a2_4x2x1w);
    });

    measure_cycles("- mul_fp6x2_1x2x2x2w:   ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_fp6x2_1x2x2x2w(z01_2x2x2w, z2_2x2x2w, a0_2x2x2w, a1_2x2x2w, a2_2x2x2w);
    });

    measure_cycles("- mul_by_xy00z0_fp6x2_2x2x2x1w: ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_by_xy00z0_fp6x2_2x2x2x1w(
            z01_4x2x1w,
            z2_4x2x1w,
            z3_4x2x1w,
            z45_4x2x1w,
            a01_4x2x1w,
            a2_pair_4x2x1w,
            b01_4x2x1w,
            b4_4x2x1w);
    });

    measure_cycles("- mul_by_xy00z0_fp6x2_2x4x1x1w: ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_by_xy00z0_fp6x2_2x4x1x1w(r01_8x1x1w, r23_8x1x1w, r45_8x1x1w, a_8x1x1w, b_8x1x1w);
    });
}

void bench_fp12(const BenchConfig& cfg)
{
    Fp12 a_scalar{};
    Fp12 r_scalar{};
    Fp6 b_sparse{};

    Fp4Pack1x2x2x2 ra_1x2x2x2w{};
    Fp4Pack1x2x2x2 a_1x2x2x2w{};
    Fp4Pack2x2x2x1 rbc_2x2x2x1w{};
    Fp4Pack2x2x2x1 bc_2x2x2x1w{};
    Fp2Pack8x1x1 a0_8x1x1w{};
    Fp2Pack8x1x1 a1_8x1x1w{};
    Fp2Pack8x1x1 a2_8x1x1w{};
    Fp2Pack8x1x1 a_8x1x1w{};
    Fp2Pack8x1x1 b_8x1x1w{};
    Fp2Pack4x2x1 r01_4x2x1w{};
    Fp2Pack4x2x1 r2_4x2x1w{};
    Fp2Pack4x2x1 a0_4x2x1w{};
    Fp2Pack4x2x1 a1_4x2x1w{};
    Fp2Pack4x2x1 a2_4x2x1w{};
    Fp2Pack2x2x2 r001_2x2x2w{};
    Fp2Pack2x2x2 r101_2x2x2w{};
    Fp2Pack2x2x2 r2_2x2x2w{};
    Fp2Pack2x2x2 r12_2x2x2w{};
    Fp2Pack4x2x1 r0_4x2x1w{};
    Fp2Pack4x2x1 r1_4x2x1w{};
    Fp2Pack4x2x1 a01_4x2x1w{};
    Fp2Pack4x2x1 b01_4x2x1w{};
    Fp2Pack4x2x1 b4_4x2x1w{};
    Fp2Pack2x2x2 r1_2x2x2w{};

    fill_scalar_fp12(a_scalar, 19);
    fill_scalar_fp6(b_sparse, 311);

    fill_vec_words(ra_1x2x2x2w, VWORDS, 901);
    fill_vec_words(a_1x2x2x2w, VWORDS, 931);
    fill_vec_words(rbc_2x2x2x1w, NWORDS, 961);
    fill_vec_words(bc_2x2x2x1w, NWORDS, 991);
    fill_vec_words(a0_8x1x1w[0], NWORDS, 1021);
    fill_vec_words(a0_8x1x1w[1], NWORDS, 1051);
    fill_vec_words(a1_8x1x1w[0], NWORDS, 1081);
    fill_vec_words(a1_8x1x1w[1], NWORDS, 1111);
    fill_vec_words(a2_8x1x1w[0], NWORDS, 1141);
    fill_vec_words(a2_8x1x1w[1], NWORDS, 1171);
    fill_vec_words(a_8x1x1w[0], NWORDS, 1201);
    fill_vec_words(a_8x1x1w[1], NWORDS, 1231);
    fill_vec_words(b_8x1x1w[0], NWORDS, 1261);
    fill_vec_words(b_8x1x1w[1], NWORDS, 1291);

    fill_vec_words(r01_4x2x1w, NWORDS, 1321);
    fill_vec_words(r2_4x2x1w, NWORDS, 1351);
    fill_vec_words(a0_4x2x1w, NWORDS, 1381);
    fill_vec_words(a1_4x2x1w, NWORDS, 1411);
    fill_vec_words(a2_4x2x1w, NWORDS, 1441);
    fill_vec_words(r001_2x2x2w, VWORDS, 1471);
    fill_vec_words(r101_2x2x2w, VWORDS, 1501);
    fill_vec_words(r2_2x2x2w, VWORDS, 1531);
    fill_vec_words(r12_2x2x2w, VWORDS, 1561);
    fill_vec_words(r0_4x2x1w, NWORDS, 1591);
    fill_vec_words(r1_4x2x1w, NWORDS, 1621);
    fill_vec_words(a01_4x2x1w, NWORDS, 1651);
    fill_vec_words(b01_4x2x1w, NWORDS, 1681);
    fill_vec_words(b4_4x2x1w, NWORDS, 1711);
    fill_vec_words(r1_2x2x2w, VWORDS, 1741);

    std::puts("\n=============================================================\n");
    std::puts("TIMING - FP12\n");

    measure_cycles("- cyclotomic_sqr_fp12_scalar:            ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        Fp12::cyclotomic_sqr(r_scalar, a_scalar);
    });

    measure_cycles("- cyclotomic_sqr_fp12_vec_v1:            ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        cyclotomic_sqr_fp12_vec_v1(ra_1x2x2x2w, rbc_2x2x2x1w, a_1x2x2x2w, bc_2x2x2x1w);
    });

    measure_cycles("- cyclotomic_sqr_fp12_vec_v2:            ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        cyclotomic_sqr_fp12_vec_v2(ra_1x2x2x2w, rbc_2x2x2x1w, a_1x2x2x2w, bc_2x2x2x1w);
    });

    measure_cycles("- mul_fp12_scalar:                       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        Fp12::mul(r_scalar, a_scalar, a_scalar);
    });

    measure_cycles("- mul_fp12_vec_v1:                       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_fp12_vec_v1(r01_4x2x1w, r2_4x2x1w, a0_8x1x1w, a1_8x1x1w, a2_8x1x1w);
    });

    measure_cycles("- mul_fp12_vec_v2:                       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_fp12_vec_v2(r01_4x2x1w, r2_2x2x2w, a0_8x1x1w, a1_8x1x1w, a2_8x1x1w);
    });

    measure_cycles("- mul_fp12_vec_v3:                       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_fp12_vec_v3(r001_2x2x2w, r101_2x2x2w, r2_2x2x2w, a0_4x2x1w, a1_4x2x1w, a2_4x2x1w);
    });

    measure_cycles("- mul_fp12_vec_v4:                       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_fp12_vec_v4(r0_4x2x1w, r101_2x2x2w, r12_2x2x2w, a0_4x2x1w, a1_4x2x1w, a2_4x2x1w);
    });

    measure_cycles("- mul_by_xy00z0_fp12_scalar:             ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        Fp12::mul_by_xy00z0(r_scalar, a_scalar, b_sparse);
    });

    measure_cycles("- mul_by_xy00z0_fp12_vec_v1:             ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_by_xy00z0_fp12_vec_v1(r0_4x2x1w, r1_2x2x2w, a01_4x2x1w, a2_4x2x1w, b01_4x2x1w, b4_4x2x1w);
    });

    measure_cycles("- mul_by_xy00z0_fp12_vec_v2:             ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_by_xy00z0_fp12_vec_v2(r0_4x2x1w, r1_2x2x2w, a_8x1x1w, b_8x1x1w);
    });

    measure_cycles("- sqr_fp12_scalar:                       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        Fp12::sqr(r_scalar, a_scalar);
    });

    measure_cycles("- sqr_fp12_vec_v1:                       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        sqr_fp12_vec_v1(r0_4x2x1w, r1_4x2x1w, a0_4x2x1w, a1_4x2x1w);
    });

    measure_cycles("- mul_fp12_vector:                       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_fp12_vector(r_scalar, a_scalar, a_scalar);
    });

    measure_cycles("- mul_by_xy00z0_fp12_vector:             ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        mul_by_xy00z0_fp12_vector(r_scalar, a_scalar, b_sparse);
    });

    measure_cycles("- sqr_fp12_vector:                       ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        sqr_fp12_vector(r_scalar, a_scalar);
    });

    measure_cycles("- cyclotomic_sqr_fp12_vector:            ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        cyclotomic_sqr_fp12_vector(r_scalar, a_scalar);
    });
}

void bench_line(const BenchConfig& cfg)
{
    JacobianPoint q_scalar{};
    AffineInputPoint r_scalar{};
    JacobianPoint t_scalar{};
    LineEval line_scalar{};

    Fp2Pack2x2x2 X1Y1{};
    Fp2Pack2x2x2 Z1{};
    Fp2Pack2x2x2 Z1Y2{};
    Fp2Pack2x2x2 X2{};
    Fp2Pack2x2x2 l0Y3{};
    Fp2Pack2x2x2 l1{};
    Fp2Pack2x2x2 l2{};
    Fp2Pack2x2x2 X3{};
    Fp2Pack2x2x2 Z3{};

    Fp2Pack4x2x1 X1Y1Z1{};
    Fp2Pack4x2x1 l0{};
    Fp2Pack4x2x1 l12{};
    Fp2Pack4x2x1 _X3{};
    Fp2Pack4x2x1 _Y3{};
    Fp2Pack4x2x1 _Z3{};

    init_line_inputs(q_scalar, r_scalar, X1Y1, Z1, Z1Y2, X2, X1Y1Z1);

    std::puts("\n=============================================================\n");
    std::puts("TIMING - LINE\n");

    measure_cycles("- line_dbl_scalar: ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        line_dbl_scalar(line_scalar, t_scalar, q_scalar);
    });

    measure_cycles("- line_dbl_vec_v1: ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        line_dbl_vec_v1(l0Y3, l1, l2, X3, Z3, X1Y1, Z1);
    });

    measure_cycles("- line_dbl_vec_v2: ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        line_dbl_vec_v2(l0, l12, _X3, _Y3, _Z3, X1Y1Z1);
    });

    measure_cycles("- line_add_scalar: ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        line_add_scalar(line_scalar, t_scalar, q_scalar, r_scalar);
    });

    measure_cycles("- line_add_vec_v1: ", cfg.fp_warmup, cfg.fp_iters, [&]() {
        line_add_vec_v1(l0Y3, l1, X3, Z3, X1Y1, Z1Y2, X2);
    });
}

bool parse_config(int argc, char** argv, BenchConfig& cfg)
{
    for (int i = 3; i < argc; ++i) {
        if (std::strcmp(argv[i], "--iters") == 0) {
            if (i + 1 >= argc || !parse_u64(argv[i + 1], &cfg.fp_iters)) {
                std::fprintf(stderr, "invalid --iters value\n");
                return false;
            }
            cfg.pairing_iters = cfg.fp_iters < 10 ? 1 : (cfg.fp_iters / 10);
            ++i;
            continue;
        }

        if (std::strcmp(argv[i], "--warmup") == 0) {
            if (i + 1 >= argc || !parse_u64(argv[i + 1], &cfg.fp_warmup)) {
                std::fprintf(stderr, "invalid --warmup value\n");
                return false;
            }
            cfg.pairing_warmup = cfg.fp_warmup < 10 ? 1 : (cfg.fp_warmup / 10);
            ++i;
            continue;
        }

        if (std::strcmp(argv[i], "--help") == 0) {
            std::puts("usage: blst_benchmarks --user --avxbls [--iters N] [--warmup N]");
            std::puts("defaults: fp/group warmup=10000 iters=100000, pairing warmup=1000 iters=10000");
            return false;
        }

        std::fprintf(stderr, "unknown flag for --avxbls mode: %s\n", argv[i]);
        return false;
    }
    return true;
}
} // namespace

int run_avxbls_cycle_bench(int argc, char** argv)
{
    (void)argc;
    (void)argv;

    BenchConfig cfg{};
    if (!parse_config(argc, argv, cfg)) {
        return 1;
    }

    bench_pairing(cfg);
    bench_fp(cfg);
    bench_fp2(cfg);
    bench_fp4(cfg);
    bench_fp6(cfg);
    bench_fp12(cfg);
    bench_line(cfg);
    return 0;
}

#else

int run_avxbls_cycle_bench(int argc, char** argv)
{
    (void)argc;
    (void)argv;
    std::puts("AVX benchmark mode requires BLST_AVX512 or BLST_AVX512_PORTABLE.");
    return 1;
}

#endif

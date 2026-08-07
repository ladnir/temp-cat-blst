// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs

#include "Matrixnok_fft_preview.h"

#include <benchmark/benchmark.h>

#include "blst/backend/backend.h"
#include "blst/avx/curve/matrixnok_g1.h"
#include "blst/avx/curve/point.h"
#include "blst/avx/fft/gt_context.h"
#include "blst/common/bytes.h"
#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/primitveRoot.h"
#include "blst/scalar/curve/GT.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace
{
    using AvxG1 = blst::avx::G1Point;
    using AvxGt = blst::avx::AvxGt52;
    namespace matrixnok = blst::avx::matrixnok;

    std::vector<blst::Fr> make_roots(const std::size_t n, const bool inverse)
    {
        blst::Fr root = blst::primRootOfUnity<blst::Fr>(n);
        if (inverse) {
            blst::Fr::reciprocal(root, root);
        }
        std::vector<blst::Fr> roots(n);
        blst::Fr current = blst::Fr::one();
        for (std::size_t i = 0; i < n; ++i) {
            roots[i] = current;
            blst::Fr::mul(current, current, root);
        }
        return roots;
    }

    std::array<unsigned char, 32> scalar_bytes(const blst::Fr& value)
    {
        blst::Fr canonical = {};
        blst::Fr::from_mont(canonical, value);
        std::array<unsigned char, 32> out = {};
        blst::le_bytes_from_limbs(out.data(), canonical.l, out.size());
        return out;
    }

    struct AvxG1FftCtx
    {
        static AvxG1 plus(const AvxG1& left, const AvxG1& right)
        {
            AvxG1 out;
            AvxG1::dadd(out, left, right);
            return out;
        }

        static AvxG1 minus(const AvxG1& left, const AvxG1& right)
        {
            AvxG1 negative;
            AvxG1::neg(negative, right);
            return plus(left, negative);
        }

        static AvxG1 mult(const AvxG1& point, const blst::Fr& scalar)
        {
            const auto bytes = scalar_bytes(scalar);
            AvxG1 out;
            AvxG1::mult_vartime(out, point, bytes.data(), 256);
            return out;
        }

        static blst::Fr inverse(const blst::Fr& value)
        {
            blst::Fr out;
            blst::Fr::reciprocal(out, value);
            return out;
        }

        template<class>
        static blst::Fr from_u64(const uint64_t value)
        {
            return blst::Fr::from_uint64(value);
        }
    };

    struct MatrixNoKG1FftCtx
    {
        inline static const matrixnok::FpRingCandidate ring = {};
        inline static const matrixnok::G1CurveCandidate curve = {};
        inline static const matrixnok::G1GLVCandidate glv = {curve, ring};

        static matrixnok::G1ProjectiveCandidate plus(
            const matrixnok::G1ProjectiveCandidate& left,
            const matrixnok::G1ProjectiveCandidate& right)
        {
            return curve.add_point(left, right, ring);
        }

        static matrixnok::G1ProjectiveCandidate minus(
            const matrixnok::G1ProjectiveCandidate& left,
            const matrixnok::G1ProjectiveCandidate& right)
        {
            return plus(left, curve.negate(right, ring));
        }

        static matrixnok::G1ProjectiveCandidate mult(
            const matrixnok::G1ProjectiveCandidate& point,
            const blst::Fr& scalar)
        {
            return glv.multiply_vartime(point, scalar_bytes(scalar), ring);
        }

        static blst::Fr inverse(const blst::Fr& value)
        {
            blst::Fr out;
            blst::Fr::reciprocal(out, value);
            return out;
        }

        template<class>
        static blst::Fr from_u64(const uint64_t value)
        {
            return blst::Fr::from_uint64(value);
        }
    };

    std::vector<AvxG1> make_g1_data(const std::size_t n)
    {
        std::vector<AvxG1> out(n);
        const AvxG1 generator = AvxG1::generator();
        AvxG1 current = generator;
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = current;
            current = AvxG1FftCtx::plus(current, generator);
        }
        return out;
    }

    matrixnok::G1ProjectiveCandidate make_matrixnok_g1_generator()
    {
        const auto affine = blst::scalar::G1AffinePoint::generator();
        return matrixnok::G1ProjectiveCandidate(
            MatrixNoKG1FftCtx::ring.from_scalar(affine.x),
            MatrixNoKG1FftCtx::ring.from_scalar(affine.y),
            MatrixNoKG1FftCtx::ring.one());
    }

    std::vector<matrixnok::G1ProjectiveCandidate> make_matrixnok_g1_data(
        const std::size_t n)
    {
        std::vector<matrixnok::G1ProjectiveCandidate> out(n);
        const auto generator = make_matrixnok_g1_generator();
        auto current = generator;
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = current;
            current = MatrixNoKG1FftCtx::plus(current, generator);
        }
        return out;
    }

    void BM_CatBLST_MatrixNoKCandidate_G1Add(benchmark::State& state)
    {
        const auto generator = make_matrixnok_g1_generator();
        auto current = generator;
        for (auto _ : state) {
            current = MatrixNoKG1FftCtx::curve.add_point(
                current, generator, MatrixNoKG1FftCtx::ring);
            benchmark::DoNotOptimize(current);
        }
    }

    void BM_CatBLST_MatrixNoKCandidate_G1Double(benchmark::State& state)
    {
        auto current = make_matrixnok_g1_generator();
        for (auto _ : state) {
            current = MatrixNoKG1FftCtx::curve.double_point(
                current, MatrixNoKG1FftCtx::ring);
            benchmark::DoNotOptimize(current);
        }
    }

    std::vector<AvxGt> make_gt_data(const std::size_t n)
    {
        std::vector<AvxGt> out(n);
        const AvxGt identity = AvxGt::identity();
        const AvxGt generator = AvxGt::from_scalar(blst::scalar::GT::generator());
        out[0] = identity;
        AvxGt current = generator;
        for (std::size_t i = 1; i < n; ++i) {
            out[i] = current;
            current = blst::avx::AvxGt52CoeffCtx::plus(current, generator);
        }
        return out;
    }

    void BM_CatBLST_Incumbent_G1Forward(benchmark::State& state)
    {
        const std::size_t n = static_cast<std::size_t>(state.range(0));
        const auto roots = make_roots(n, false);
        const auto original = make_g1_data(n);
        std::vector<AvxG1> data(n);

        for (auto _ : state) {
            state.PauseTiming();
            data = original;
            state.ResumeTiming();
            const auto start = std::chrono::steady_clock::now();
            fft<AvxG1, blst::Fr, AvxG1FftCtx>(data, roots, false);
            const auto finish = std::chrono::steady_clock::now();
            state.SetIterationTime(std::chrono::duration<double>(finish - start).count());
            benchmark::DoNotOptimize(data.data());
        }
    }

    void BM_CatBLST_MatrixNoKCandidate_G1Forward(benchmark::State& state)
    {
        const std::size_t n = static_cast<std::size_t>(state.range(0));
        const auto roots = make_roots(n, false);
        const auto original = make_matrixnok_g1_data(n);
        std::vector<matrixnok::G1ProjectiveCandidate> data(n);

        for (auto _ : state) {
            state.PauseTiming();
            data = original;
            state.ResumeTiming();
            const auto start = std::chrono::steady_clock::now();
            fft<matrixnok::G1ProjectiveCandidate, blst::Fr, MatrixNoKG1FftCtx>(
                data, roots, false);
            const auto finish = std::chrono::steady_clock::now();
            state.SetIterationTime(std::chrono::duration<double>(finish - start).count());
            benchmark::DoNotOptimize(data.data());
        }
    }

    void BM_CatBLST_Incumbent_GTInverse(benchmark::State& state)
    {
        const std::size_t n = static_cast<std::size_t>(state.range(0));
        const auto roots = make_roots(n, false);
        const auto inverse_roots = make_roots(n, true);
        auto frequency = make_gt_data(n);
        {
            blst::avx::AvxGt52PreparedRoots prepared(roots);
            blst::avx::AvxGt52CoeffCtx::ScopedRoots scope(prepared);
            fft<AvxGt, blst::Fr, blst::avx::AvxGt52CoeffCtx>(frequency, roots, false);
        }
        blst::avx::AvxGt52PreparedRoots prepared_inverse(inverse_roots);
        prepared_inverse.set_scale_from_size(n);
        blst::avx::AvxGt52CoeffCtx::ScopedRoots scope(prepared_inverse);
        std::vector<AvxGt> data(n);

        for (auto _ : state) {
            state.PauseTiming();
            data = frequency;
            state.ResumeTiming();
            const auto start = std::chrono::steady_clock::now();
            fft<AvxGt, blst::Fr, blst::avx::AvxGt52CoeffCtx>(
                data, inverse_roots, true);
            const auto finish = std::chrono::steady_clock::now();
            state.SetIterationTime(std::chrono::duration<double>(finish - start).count());
            benchmark::DoNotOptimize(data.data());
        }
    }

    BENCHMARK(BM_CatBLST_Incumbent_G1Forward)
        ->Arg(256)->Arg(1024)->Arg(4096)->Arg(16384)
        ->UseManualTime()->Unit(benchmark::kMillisecond);
    BENCHMARK(BM_CatBLST_MatrixNoKCandidate_G1Forward)
        ->Arg(256)->Arg(1024)->Arg(4096)->Arg(16384)
        ->UseManualTime()->Unit(benchmark::kMillisecond);
    BENCHMARK(BM_CatBLST_MatrixNoKCandidate_G1Add);
    BENCHMARK(BM_CatBLST_MatrixNoKCandidate_G1Double);
    BENCHMARK(BM_CatBLST_Incumbent_GTInverse)
        ->Arg(256)->Arg(1024)->Arg(4096)->Arg(16384)
        ->UseManualTime()->Unit(benchmark::kMillisecond);
}

extern "C" bool catblst_matrixnok_fft_roots(
    unsigned char* roots_le,
    MatrixNoKFftGlv4Schedule* schedules,
    const std::size_t n,
    const bool inverse,
    unsigned char scale_le[32],
    MatrixNoKFftGlv4Schedule* scale_schedule)
{
    if (roots_le == nullptr || schedules == nullptr || scale_le == nullptr ||
        scale_schedule == nullptr || n == 0 || (n & (n - 1)) != 0) {
        return false;
    }
    const auto roots = make_roots(n, inverse);
    for (std::size_t i = 0; i < n; ++i) {
        const auto bytes = scalar_bytes(roots[i]);
        std::copy(bytes.begin(), bytes.end(), roots_le + i * bytes.size());
        std::array<uint64_t, 4> limbs = {};
        blst::limbs_from_le_bytes(limbs.data(), bytes.data(), bytes.size());
        const blst::avx::AvxGLV4Exp compiled = blst::avx::avx_glv4_compile(limbs.data());
        schedules[i].sign_mask = compiled.sign_mask;
        schedules[i].sched_len = compiled.sched_len;
        schedules[i].schedule = compiled.schedule;
    }

    blst::Fr scale = blst::Fr::from_uint64(n);
    blst::Fr::reciprocal(scale, scale);
    const auto scale_bytes = scalar_bytes(scale);
    std::copy(scale_bytes.begin(), scale_bytes.end(), scale_le);
    std::array<uint64_t, 4> scale_limbs = {};
    blst::limbs_from_le_bytes(scale_limbs.data(), scale_bytes.data(), scale_bytes.size());
    const blst::avx::AvxGLV4Exp compiled_scale =
        blst::avx::avx_glv4_compile(scale_limbs.data());
    scale_schedule->sign_mask = compiled_scale.sign_mask;
    scale_schedule->sched_len = compiled_scale.sched_len;
    scale_schedule->schedule = compiled_scale.schedule;
    return true;
}

extern "C" void catblst_matrixnok_gt_generator_coefficients(
    unsigned char coefficients_be[12 * 48])
{
    const blst::Fp12& value = blst::scalar::GT::generator().data;
    const std::array<const blst::Fp*, 12> coefficients = {
        &value[0][0][0], &value[0][0][1],
        &value[0][1][0], &value[0][1][1],
        &value[0][2][0], &value[0][2][1],
        &value[1][0][0], &value[1][0][1],
        &value[1][1][0], &value[1][1][1],
        &value[1][2][0], &value[1][2][1],
    };
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        blst::Fp canonical = {};
        blst::Fp::from_mont(canonical, *coefficients[i]);
        blst::be_bytes_from_limbs(coefficients_be + i * 48, canonical.l, 48);
    }
}

extern "C" void catblst_matrixnok_g1_generator_coordinates(
    unsigned char coordinates_be[2 * 48])
{
    blst::scalar::G1AffinePoint::serialize_BE(
        coordinates_be, blst::scalar::G1AffinePoint::generator());
}

// Pinned VROOM still names the inherited BLST division helper as a C symbol.
// Keep that compatibility shim confined to this benchmark executable.
extern "C" void blst_div_by_zz(std::uint64_t val[])
{
    blst::div_by_zz(val);
}

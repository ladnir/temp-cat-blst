#pragma once

#include <algorithm>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "blst/common/defines.h"
#include "blst/common/timer.h"
#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/primitveRoot.h"
#include "blst/fields/Fr.h"

namespace blst {

enum class MiddleProductAlgorithm {
    BasicFFT,
    SpecializedFFT,
};

struct MiddleProductPrecomputedB {
    MiddleProductAlgorithm algorithm;
    std::size_t n;
    std::size_t fft_size;
    std::vector<G2Point> first;
    std::vector<G2Point> second;
};

inline std::size_t middle_product_fft_size(std::size_t n, MiddleProductAlgorithm algorithm)
{
    BLST_ASSERT(n > 0);

    std::size_t target = 0;
    if (algorithm == MiddleProductAlgorithm::BasicFFT) {
        target = (3 * n) - 1;
    }
    else {
        target = 2 * n;
    }

    std::size_t fft_size = 1;
    while (fft_size < target) {
        fft_size <<= 1;
    }
    return fft_size;
}

inline std::size_t fft_admissible_size(std::size_t x, std::size_t fft_size)
{
    BLST_ASSERT(fft_size > 0);
    BLST_ASSERT((fft_size & (fft_size - 1)) == 0);

    constexpr std::size_t kRoundPow2 = 16;
    if (x == 0) {
        x = 1;
    }

    x = ((x + (kRoundPow2 - 1)) / kRoundPow2) * kRoundPow2;
    if (x > fft_size) {
        x = fft_size;
    }

    // Match NTL's near-full truncation heuristic: for small sizes round up
    // earlier (n - n/8), for larger sizes later (n - n/16).
    const std::size_t near_full =
        (fft_size >= 1024) ? (fft_size - (fft_size >> 4)) : (fft_size - (fft_size >> 3));
    if (x > near_full) {
        x = fft_size;
    }

    return x;
}

inline std::vector<Fr> make_fft_roots(std::size_t n)
{
    Fr prim = primRootOfUnity<Fr>(n);
    std::vector<Fr> roots(n);
    Fr cur = Fr::one();
    for (std::size_t i = 0; i < n; ++i) {
        roots[i] = cur;
        Fr::mul(cur, cur, prim);
    }
    return roots;
}

inline std::vector<Fr> make_ifft_roots(std::span<const Fr> roots)
{
    std::vector<Fr> inv_roots(roots.size());
    for (std::size_t i = 0; i < roots.size(); ++i) {
        Fr::reciprocal(inv_roots[i], roots[i]);
    }
    return inv_roots;
}

inline MiddleProductPrecomputedB precompute_middle_product_b(
    std::span<const G2Point> b,
    std::size_t n,
    MiddleProductAlgorithm algorithm,
    TimingStats* timings = nullptr)
{
    ScopedTimer total_timer(timings, "precompute.total_ns");
    BLST_ASSERT(n > 0);
    BLST_ASSERT(b.size() == 2 * n);

    MiddleProductPrecomputedB precomputed;
    precomputed.algorithm = algorithm;
    precomputed.n = n;
    precomputed.fft_size = middle_product_fft_size(n, algorithm);
    precomputed.first.assign(precomputed.fft_size, G2Point::infinity());

    std::vector<Fr> roots;
    {
        ScopedTimer timer(timings, "precompute.roots_ns");
        roots = make_fft_roots(precomputed.fft_size);
    }

    if (algorithm == MiddleProductAlgorithm::BasicFFT) {
        const std::size_t xn_b = fft_admissible_size(2 * n, precomputed.fft_size);
        for (std::size_t i = 0; i < 2 * n; ++i) {
            precomputed.first[i] = b[i];
        }
        ScopedTimer timer(timings, "precompute.basic_fft_b_ns");
        fft_forward_truncated_prefix<G2Point, Fr, G2CoeffCtx>(
            precomputed.first, roots, xn_b, timings, "precompute.fft_b");
        return precomputed;
    }

    precomputed.second.assign(precomputed.fft_size, G2Point::infinity());
    for (std::size_t i = 0; i < n; ++i) {
        precomputed.second[i] = b[i];
        precomputed.first[i] = b[n + i];
    }
    for (std::size_t i = 0; i < n / 2; ++i) {
        std::swap(precomputed.second[i], precomputed.second[n - 1 - i]);
    }

    {
        ScopedTimer timer(timings, "precompute.spec_fft_b1_ns");
        fft<G2Point, Fr, G2CoeffCtx>(precomputed.first, roots, false);
    }
    {
        ScopedTimer timer(timings, "precompute.spec_fft_b0_ns");
        fft<G2Point, Fr, G2CoeffCtx>(precomputed.second, roots, false);
    }
    return precomputed;
}

inline void middle_product(
    std::span<const G1Point> a,
    std::span<const G2Point> b,
    std::span<GT> out,
    MiddleProductAlgorithm algorithm,
    TimingStats* timings = nullptr)
{
    ScopedTimer total_timer(
        timings,
        algorithm == MiddleProductAlgorithm::BasicFFT ? "basic.total_ns" : "middle.total_ns");
    const std::size_t n = a.size();
    BLST_ASSERT(n > 0);
    BLST_ASSERT(b.size() == 2 * n);
    BLST_ASSERT(out.size() >= n);

    const std::size_t fft_size = middle_product_fft_size(n, algorithm);
    std::vector<Fr> roots;
    std::vector<Fr> inv_roots;
    {
        ScopedTimer timer(
            timings,
            algorithm == MiddleProductAlgorithm::BasicFFT ? "basic.roots_ns" : "middle.roots_ns");
        roots = make_fft_roots(fft_size);
    }
    {
        ScopedTimer timer(
            timings,
            algorithm == MiddleProductAlgorithm::BasicFFT ? "basic.inv_roots_ns" : "middle.inv_roots_ns");
        inv_roots = make_ifft_roots(roots);
    }

    std::fill(out.begin(), out.begin() + static_cast<std::ptrdiff_t>(n), GT::identity());

    if (algorithm == MiddleProductAlgorithm::BasicFFT) {
        const std::size_t xn_a = fft_admissible_size(n, fft_size);
        const std::size_t xn_b = fft_admissible_size(2 * n, fft_size);
        std::vector<G1Point> a_fft(fft_size, G1Point::infinity());
        std::vector<G2Point> b_fft(fft_size, G2Point::infinity());
        std::vector<GT> c_fft(fft_size, GT::identity());

        for (std::size_t i = 0; i < n; ++i) {
            a_fft[i] = a[i];
        }
        for (std::size_t i = 0; i < 2 * n; ++i) {
            b_fft[i] = b[i];
        }

        poly_mult_fft<G1Point, G2Point, GT, Fr, CoeffCtx>(
            a_fft, b_fft, false, c_fft, roots, inv_roots, timings, xn_a, xn_b);

        {
            ScopedTimer timer(timings, "basic.extract_ns");
            for (std::size_t i = 0; i < n; ++i) {
                out[i] = c_fft[n + i];
            }
        }
        return;
    }

    std::vector<G1Point> a_fft(fft_size, G1Point::infinity());
    std::vector<G2Point> b0_fft(fft_size, G2Point::infinity());
    std::vector<G2Point> b1_fft(fft_size, G2Point::infinity());
    std::vector<GT> scratch(fft_size, GT::identity());

    for (std::size_t i = 0; i < n; ++i) {
        a_fft[i] = a[i];
        b0_fft[i] = b[i];
        b1_fft[i] = b[n + i];
    }

    const Fr prim = primRootOfUnity<Fr>(fft_size);
    middle_product_fft<G1Point, G2Point, GT, Fr, CoeffCtx>(
        a_fft, b0_fft, b1_fft, false, n, fft_size, out, scratch, roots, inv_roots, prim,
        GT::identity(), G1Point::infinity(), G2Point::infinity(), timings);
}

inline void middle_product(
    std::span<const G1Point> a,
    const MiddleProductPrecomputedB& precomputed_b,
    std::span<GT> out,
    TimingStats* timings = nullptr)
{
    ScopedTimer total_timer(
        timings,
        precomputed_b.algorithm == MiddleProductAlgorithm::BasicFFT
            ? "basic_pre.total_ns"
            : "middle_pre.total_ns");
    const std::size_t n = a.size();
    BLST_ASSERT(n > 0);
    BLST_ASSERT(precomputed_b.n == n);
    BLST_ASSERT(out.size() >= n);

    const std::size_t fft_size = precomputed_b.fft_size;
    std::vector<Fr> roots;
    std::vector<Fr> inv_roots;
    {
        ScopedTimer timer(
            timings,
            precomputed_b.algorithm == MiddleProductAlgorithm::BasicFFT
                ? "basic_pre.roots_ns"
                : "middle_pre.roots_ns");
        roots = make_fft_roots(fft_size);
    }
    {
        ScopedTimer timer(
            timings,
            precomputed_b.algorithm == MiddleProductAlgorithm::BasicFFT
                ? "basic_pre.inv_roots_ns"
                : "middle_pre.inv_roots_ns");
        inv_roots = make_ifft_roots(roots);
    }

    std::fill(out.begin(), out.begin() + static_cast<std::ptrdiff_t>(n), GT::identity());

    if (precomputed_b.algorithm == MiddleProductAlgorithm::BasicFFT) {
        const std::size_t xn_a = fft_admissible_size(n, fft_size);
        const std::size_t xn_b = fft_admissible_size(2 * n, fft_size);
        BLST_ASSERT(precomputed_b.first.size() == fft_size);

        std::vector<G1Point> a_fft(fft_size, G1Point::infinity());
        std::vector<G2Point> b_fft = precomputed_b.first;
        std::vector<GT> c_fft(fft_size, GT::identity());

        for (std::size_t i = 0; i < n; ++i) {
            a_fft[i] = a[i];
        }

        poly_mult_fft<G1Point, G2Point, GT, Fr, CoeffCtx>(
            a_fft, b_fft, true, c_fft, roots, inv_roots, timings, xn_a, xn_b);

        {
            ScopedTimer timer(timings, "basic_pre.extract_ns");
            for (std::size_t i = 0; i < n; ++i) {
                out[i] = c_fft[n + i];
            }
        }
        return;
    }

    BLST_ASSERT(precomputed_b.first.size() == fft_size);
    BLST_ASSERT(precomputed_b.second.size() == fft_size);

    std::vector<G1Point> a_fft(fft_size, G1Point::infinity());
    std::vector<G2Point> b0_fft = precomputed_b.second;
    std::vector<G2Point> b1_fft = precomputed_b.first;
    std::vector<GT> scratch(fft_size, GT::identity());

    for (std::size_t i = 0; i < n; ++i) {
        a_fft[i] = a[i];
    }

    const Fr prim = primRootOfUnity<Fr>(fft_size);
    middle_product_fft<G1Point, G2Point, GT, Fr, CoeffCtx>(
        a_fft, b0_fft, b1_fft, true, n, fft_size, out, scratch, roots, inv_roots, prim,
        GT::identity(), G1Point::infinity(), G2Point::infinity(), timings);
}

// Standalone baseline variant for experimentation:
// Basic FFT path with full forward FFT and truncated inverse FFT over middle window.
inline void middle_product_basic_ifft_truncated(
    std::span<const G1Point> a,
    const MiddleProductPrecomputedB& precomputed_b,
    std::span<GT> out,
    TimingStats* timings = nullptr)
{
    ScopedTimer total_timer(timings, "basic_ifft_trunc.total_ns");
    const std::size_t n = a.size();
    BLST_ASSERT(n > 0);
    BLST_ASSERT(precomputed_b.algorithm == MiddleProductAlgorithm::BasicFFT);
    BLST_ASSERT(precomputed_b.n == n);
    BLST_ASSERT(out.size() >= n);

    const std::size_t fft_size = precomputed_b.fft_size;
    BLST_ASSERT(precomputed_b.first.size() == fft_size);

    std::vector<Fr> roots;
    std::vector<Fr> inv_roots;
    {
        ScopedTimer timer(timings, "basic_ifft_trunc.roots_ns");
        roots = make_fft_roots(fft_size);
    }
    {
        ScopedTimer timer(timings, "basic_ifft_trunc.inv_roots_ns");
        inv_roots = make_ifft_roots(roots);
    }

    std::vector<G1Point> a_fft(fft_size, G1Point::infinity());
    std::vector<G2Point> b_fft = precomputed_b.first;
    std::vector<GT> c_fft(fft_size, GT::identity());

    for (std::size_t i = 0; i < n; ++i) {
        a_fft[i] = a[i];
    }

    {
        ScopedTimer timer(timings, "basic_ifft_trunc.fft_a_ns");
        fft_forward_truncated_prefix<G1Point, Fr, CoeffCtx>(
            a_fft,
            roots,
            fft_size,
            timings,
            "basic_ifft_trunc.fft_a");
    }

    {
        ScopedTimer timer(timings, "basic_ifft_trunc.pointwise_ns");
        for (std::size_t i = 0; i < fft_size; ++i) {
            c_fft[i] = CoeffCtx::mult(a_fft[i], b_fft[i]);
        }
    }

    {
        ScopedTimer timer(timings, "basic_ifft_trunc.ifft_ns");
        fft_inverse_truncated_range<GT, Fr, CoeffCtx>(
            c_fft,
            inv_roots,
            n,
            n,
            timings,
            "basic_ifft_trunc.ifft");
    }

    {
        ScopedTimer timer(timings, "basic_ifft_trunc.extract_ns");
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = c_fft[n + i];
        }
    }
}

// TODO(BLST_CUDA): mirror middle_product() and precompute_middle_product_b() on GPU.

} // namespace blst

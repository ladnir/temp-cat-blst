#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
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
    OneIfftCyclic,
};

enum class OneIfftPairingMode {
    Direct,
    PrecomputedLines,
};

struct MiddleProductPrecomputedB {
    MiddleProductAlgorithm algorithm;
    std::size_t n;
    std::size_t fft_size;
    OneIfftPairingMode oneifft_pairing_mode = OneIfftPairingMode::Direct;
    std::vector<G2Point> first;
    std::vector<G2Point> second;
    std::vector<std::array<Fp6, 68>> oneifft_precomputed_lines;
    std::vector<uint8_t> oneifft_d_is_zero;
    std::vector<Fr> roots;
    std::vector<Fr> inv_roots;
    std::size_t xn_a = 0;
    FftTruncatedForwardPlan fft_a_plan;
    FftInverseTruncatedRangePlan ifft_middle_plan;
};

struct MiddleProductBasicPrunedWorkspace {
    std::vector<G1Point> a_fft;
    std::vector<GT> c_fft;
    std::vector<Fr> roots;
    std::vector<Fr> inv_roots;

    void ensure_size(std::size_t m)
    {
        if (a_fft.size() != m) {
            a_fft.assign(m, G1Point::infinity());
        }
        if (c_fft.size() != m) {
            c_fft.assign(m, GT::identity());
        }
    }
};

struct MiddleProductOneIfftWorkspace {
    std::vector<G1Point> a_fft;
    std::vector<GT> c_fft;
    std::vector<Fr> roots;
    std::vector<Fr> inv_roots;

    void ensure_size(std::size_t m)
    {
        if (a_fft.size() != m) {
            a_fft.assign(m, G1Point::infinity());
        }
        if (c_fft.size() != m) {
            c_fft.assign(m, GT::identity());
        }
    }
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
    TimingStats* timings = nullptr,
    OneIfftPairingMode oneifft_pairing_mode = OneIfftPairingMode::Direct)
{
    ScopedTimer total_timer(timings, "precompute.total_ns");
    BLST_ASSERT(n > 0);
    BLST_ASSERT(b.size() == 2 * n);

    MiddleProductPrecomputedB precomputed;
    precomputed.algorithm = algorithm;
    precomputed.n = n;
    precomputed.fft_size = middle_product_fft_size(n, algorithm);
    precomputed.oneifft_pairing_mode = oneifft_pairing_mode;
    precomputed.first.assign(precomputed.fft_size, G2Point::infinity());

    std::vector<Fr> roots;
    {
        ScopedTimer timer(timings, "precompute.roots_ns");
        roots = make_fft_roots(precomputed.fft_size);
    }

    if (algorithm == MiddleProductAlgorithm::BasicFFT) {
        precomputed.roots = roots;
        precomputed.inv_roots = make_ifft_roots(roots);
        precomputed.xn_a = fft_admissible_size(n, precomputed.fft_size);
        precomputed.fft_a_plan =
            make_fft_truncated_forward_plan(precomputed.fft_size, precomputed.xn_a);
        precomputed.ifft_middle_plan =
            make_fft_inverse_truncated_range_plan(precomputed.fft_size, n, n);

        const std::size_t xn_b = fft_admissible_size(2 * n, precomputed.fft_size);
        for (std::size_t i = 0; i < 2 * n; ++i) {
            precomputed.first[i] = b[i];
        }
        ScopedTimer timer(timings, "precompute.basic_fft_b_ns");
        fft_forward_truncated_prefix<G2Point, Fr, G2CoeffCtx>(
            precomputed.first, roots, xn_b, timings, "precompute.fft_b");
        return precomputed;
    }

    if (algorithm == MiddleProductAlgorithm::OneIfftCyclic) {
        precomputed.roots = roots;
        precomputed.inv_roots = make_ifft_roots(roots);
        Fr m_inv = Fr::from_uint64(static_cast<uint64_t>(precomputed.fft_size));
        Fr::reciprocal(m_inv, m_inv);

        std::vector<G2Point> b0_fft(precomputed.fft_size, G2Point::infinity());
        std::vector<G2Point> b1_fft(precomputed.fft_size, G2Point::infinity());
        for (std::size_t i = 0; i < n; ++i) {
            b0_fft[i] = b[i];
            b1_fft[i] = b[n + i];
        }

        const std::size_t xn_b = fft_admissible_size(n, precomputed.fft_size);
        {
            ScopedTimer timer(timings, "precompute.oneifft_fft_b0_ns");
            fft_forward_truncated_prefix<G2Point, Fr, G2CoeffCtx>(
                b0_fft, roots, xn_b, timings, "precompute.oneifft.fft_b0");
        }
        {
            ScopedTimer timer(timings, "precompute.oneifft_fft_b1_ns");
            fft_forward_truncated_prefix<G2Point, Fr, G2CoeffCtx>(
                b1_fft, roots, xn_b, timings, "precompute.oneifft.fft_b1");
        }
        {
            ScopedTimer timer(timings, "precompute.oneifft_fold_ns");
            for (std::size_t j = 0; j < precomputed.fft_size; ++j) {
                const std::size_t idx = (j * (n % precomputed.fft_size)) % precomputed.fft_size;
                const Fr zeta = roots[idx];
                const G2Point zeta_b1 = G2CoeffCtx::mult(b1_fft[j], zeta);
                const G2Point folded = G2CoeffCtx::plus(b0_fft[j], zeta_b1);
                // Move 1/m iFFT normalization into fixed-B precompute.
                precomputed.first[j] = G2CoeffCtx::mult(folded, m_inv);
            }
        }

        if (oneifft_pairing_mode == OneIfftPairingMode::PrecomputedLines) {
            ScopedTimer timer(timings, "precompute.oneifft_lines_ns");
            precomputed.oneifft_precomputed_lines.resize(precomputed.fft_size);
            precomputed.oneifft_d_is_zero.resize(precomputed.fft_size);
            for (std::size_t j = 0; j < precomputed.fft_size; ++j) {
                const bool d_is_zero = (G2Point::is_inf(precomputed.first[j]) != 0);
                precomputed.oneifft_d_is_zero[j] = static_cast<uint8_t>(d_is_zero ? 1 : 0);
                if (d_is_zero) {
                    continue;
                }
                G2AffinePoint q_aff;
                G2Point::to_affine(q_aff, precomputed.first[j]);
                GT::precompute_lines(precomputed.oneifft_precomputed_lines[j].data(), q_aff);
            }
        }
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
    const char* total_metric = "middle.total_ns";
    if (algorithm == MiddleProductAlgorithm::BasicFFT) {
        total_metric = "basic.total_ns";
    }
    else if (algorithm == MiddleProductAlgorithm::OneIfftCyclic) {
        total_metric = "oneifft.total_ns";
    }
    ScopedTimer total_timer(timings, total_metric);
    const std::size_t n = a.size();
    BLST_ASSERT(n > 0);
    BLST_ASSERT(b.size() == 2 * n);
    BLST_ASSERT(out.size() >= n);

    const std::size_t fft_size = middle_product_fft_size(n, algorithm);
    std::vector<Fr> roots;
    std::vector<Fr> inv_roots;
    const char* roots_metric = "middle.roots_ns";
    const char* inv_roots_metric = "middle.inv_roots_ns";
    if (algorithm == MiddleProductAlgorithm::BasicFFT) {
        roots_metric = "basic.roots_ns";
        inv_roots_metric = "basic.inv_roots_ns";
    }
    else if (algorithm == MiddleProductAlgorithm::OneIfftCyclic) {
        roots_metric = "oneifft.roots_ns";
        inv_roots_metric = "oneifft.inv_roots_ns";
    }
    {
        ScopedTimer timer(timings, roots_metric);
        roots = make_fft_roots(fft_size);
    }
    {
        ScopedTimer timer(timings, inv_roots_metric);
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

    if (algorithm == MiddleProductAlgorithm::OneIfftCyclic) {
        std::vector<G2Point> d_fft(fft_size, G2Point::infinity());
        std::vector<G1Point> a_fft(fft_size, G1Point::infinity());
        std::vector<GT> c_fft(fft_size, GT::identity());
        Fr m_inv = Fr::from_uint64(static_cast<uint64_t>(fft_size));
        Fr::reciprocal(m_inv, m_inv);
        const Fr no_scale = Fr::one();

        for (std::size_t i = 0; i < n; ++i) {
            a_fft[i] = a[i];
            d_fft[i] = b[i];
        }
        std::vector<G2Point> b1_fft(fft_size, G2Point::infinity());
        for (std::size_t i = 0; i < n; ++i) {
            b1_fft[i] = b[n + i];
        }

        const std::size_t xn_a = fft_admissible_size(n, fft_size);
        const std::size_t xn_b = fft_admissible_size(n, fft_size);
        {
            ScopedTimer timer(timings, "oneifft.fft_b0_ns");
            fft_forward_truncated_prefix<G2Point, Fr, G2CoeffCtx>(
                d_fft, roots, xn_b, timings, "oneifft.fft_b0");
        }
        {
            ScopedTimer timer(timings, "oneifft.fft_b1_ns");
            fft_forward_truncated_prefix<G2Point, Fr, G2CoeffCtx>(
                b1_fft, roots, xn_b, timings, "oneifft.fft_b1");
        }
        {
            ScopedTimer timer(timings, "oneifft.fold_ns");
            for (std::size_t j = 0; j < fft_size; ++j) {
                const std::size_t idx = (j * (n % fft_size)) % fft_size;
                const Fr zeta = roots[idx];
                const G2Point zeta_b1 = G2CoeffCtx::mult(b1_fft[j], zeta);
                const G2Point folded = G2CoeffCtx::plus(d_fft[j], zeta_b1);
                d_fft[j] = G2CoeffCtx::mult(folded, m_inv);
            }
        }

        {
            ScopedTimer timer(timings, "oneifft.fft_a_ns");
            fft_forward_truncated_prefix<G1Point, Fr, G1CoeffCtx>(
                a_fft, roots, xn_a, timings, "oneifft.fft_a");
        }
        {
            ScopedTimer timer(timings, "oneifft.pointwise_ns");
            for (std::size_t i = 0; i < fft_size; ++i) {
                c_fft[i] = CoeffCtx::mult(a_fft[i], d_fft[i]);
            }
        }
        {
            ScopedTimer timer(timings, "oneifft.ifft_ns");
            fft_inverse_full<GT, Fr, GTCoeffCtx>(
                std::span<GT>(c_fft), fft_size, std::span<const Fr>(inv_roots), no_scale);
        }
        {
            ScopedTimer timer(timings, "oneifft.extract_ns");
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

inline void middle_product_oneifft_cyclic(
    std::span<const G1Point> a,
    const MiddleProductPrecomputedB& precomputed_b,
    std::span<GT> out,
    MiddleProductOneIfftWorkspace& workspace,
    TimingStats* timings = nullptr)
{
    ScopedTimer total_timer(timings, "oneifft_pre.total_ns");
    const std::size_t n = a.size();
    BLST_ASSERT(n > 0);
    BLST_ASSERT(precomputed_b.algorithm == MiddleProductAlgorithm::OneIfftCyclic);
    BLST_ASSERT(precomputed_b.n == n);
    BLST_ASSERT(out.size() >= n);

    const std::size_t fft_size = precomputed_b.fft_size;
    BLST_ASSERT(precomputed_b.first.size() == fft_size);
    const Fr no_scale = Fr::one();

    std::span<const Fr> roots;
    std::span<const Fr> inv_roots;
    if (precomputed_b.roots.size() == fft_size &&
        precomputed_b.inv_roots.size() == fft_size) {
        roots = std::span<const Fr>(precomputed_b.roots);
        inv_roots = std::span<const Fr>(precomputed_b.inv_roots);
    }
    else {
        ScopedTimer timer(timings, "oneifft_pre.roots_ns");
        workspace.roots = make_fft_roots(fft_size);
        workspace.inv_roots = make_ifft_roots(workspace.roots);
        roots = std::span<const Fr>(workspace.roots);
        inv_roots = std::span<const Fr>(workspace.inv_roots);
    }

    workspace.ensure_size(fft_size);
    std::span<G1Point> a_fft(workspace.a_fft);
    std::span<GT> c_fft(workspace.c_fft);
    {
        ScopedTimer timer(timings, "oneifft_pre.fill_a_fft_ns");
        std::fill(a_fft.begin(), a_fft.end(), G1Point::infinity());
    }
    {
        ScopedTimer timer(timings, "oneifft_pre.fill_c_fft_ns");
        std::fill(c_fft.begin(), c_fft.end(), GT::identity());
    }

    for (std::size_t i = 0; i < n; ++i) {
        a_fft[i] = a[i];
    }

    const std::size_t xn_a = fft_admissible_size(n, fft_size);
    {
        ScopedTimer timer(timings, "oneifft_pre.fft_a_ns");
        fft_forward_truncated_prefix<G1Point, Fr, G1CoeffCtx>(
            a_fft, roots, xn_a, timings, "oneifft_pre.fft_a");
    }
    {
        ScopedTimer timer(timings, "oneifft_pre.pointwise_ns");
        if (precomputed_b.oneifft_pairing_mode == OneIfftPairingMode::PrecomputedLines &&
            precomputed_b.oneifft_precomputed_lines.size() == fft_size) {
            const bool has_d_zero_mask = (precomputed_b.oneifft_d_is_zero.size() == fft_size);
            for (std::size_t i = 0; i < fft_size; ++i) {
                if (G1Point::is_inf(a_fft[i]) != 0 ||
                    (has_d_zero_mask && precomputed_b.oneifft_d_is_zero[i] != 0)) {
                    c_fft[i] = GT::identity();
                    continue;
                }
                G1AffinePoint p_aff;
                G1Point::to_affine(p_aff, a_fft[i]);
                GT::miller_loop_lines(c_fft[i], precomputed_b.oneifft_precomputed_lines[i].data(), p_aff);
                GT::final_mult(c_fft[i], c_fft[i]);
            }
        }
        else {
            for (std::size_t i = 0; i < fft_size; ++i) {
                c_fft[i] = CoeffCtx::mult(a_fft[i], precomputed_b.first[i]);
            }
        }
    }
    {
        ScopedTimer timer(timings, "oneifft_pre.ifft_ns");
        fft_inverse_full<GT, Fr, GTCoeffCtx>(
            c_fft, fft_size, inv_roots, no_scale);
    }
    {
        ScopedTimer timer(timings, "oneifft_pre.extract_ns");
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = c_fft[n + i];
        }
    }
}

inline void middle_product_oneifft_cyclic(
    std::span<const G1Point> a,
    const MiddleProductPrecomputedB& precomputed_b,
    std::span<GT> out,
    TimingStats* timings = nullptr)
{
    MiddleProductOneIfftWorkspace workspace;
    middle_product_oneifft_cyclic(a, precomputed_b, out, workspace, timings);
}

inline void middle_product(
    std::span<const G1Point> a,
    const MiddleProductPrecomputedB& precomputed_b,
    std::span<GT> out,
    TimingStats* timings = nullptr)
{
    const char* total_metric = "middle_pre.total_ns";
    if (precomputed_b.algorithm == MiddleProductAlgorithm::BasicFFT) {
        total_metric = "basic_pre.total_ns";
    }
    else if (precomputed_b.algorithm == MiddleProductAlgorithm::OneIfftCyclic) {
        total_metric = "oneifft_pre.total_ns";
    }
    ScopedTimer total_timer(timings, total_metric);
    const std::size_t n = a.size();
    BLST_ASSERT(n > 0);
    BLST_ASSERT(precomputed_b.n == n);
    BLST_ASSERT(out.size() >= n);

    const std::size_t fft_size = precomputed_b.fft_size;
    if (precomputed_b.algorithm == MiddleProductAlgorithm::OneIfftCyclic) {
        MiddleProductOneIfftWorkspace workspace;
        middle_product_oneifft_cyclic(a, precomputed_b, out, workspace, timings);
        return;
    }

    std::vector<Fr> roots;
    std::vector<Fr> inv_roots;
    const char* roots_metric = "middle_pre.roots_ns";
    const char* inv_roots_metric = "middle_pre.inv_roots_ns";
    if (precomputed_b.algorithm == MiddleProductAlgorithm::BasicFFT) {
        roots_metric = "basic_pre.roots_ns";
        inv_roots_metric = "basic_pre.inv_roots_ns";
    }
    else if (precomputed_b.algorithm == MiddleProductAlgorithm::OneIfftCyclic) {
        roots_metric = "oneifft_pre.roots_ns";
        inv_roots_metric = "oneifft_pre.inv_roots_ns";
    }
    if (precomputed_b.roots.size() == fft_size && precomputed_b.inv_roots.size() == fft_size) {
        roots = precomputed_b.roots;
        inv_roots = precomputed_b.inv_roots;
    }
    else {
        {
            ScopedTimer timer(timings, roots_metric);
            roots = make_fft_roots(fft_size);
        }
        {
            ScopedTimer timer(timings, inv_roots_metric);
            inv_roots = make_ifft_roots(roots);
        }
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

// Dedicated online path for |A|=n, |B|=2n using:
// - truncated forward FFT on A
// - precomputed FFT(B)
// - middle-window truncated inverse FFT on C
inline void middle_product_basic_pruned(
    std::span<const G1Point> a,
    const MiddleProductPrecomputedB& precomputed_b,
    std::span<GT> out,
    MiddleProductBasicPrunedWorkspace& workspace,
    TimingStats* timings = nullptr)
{
    ScopedTimer total_timer(timings, "basic_pruned.total_ns");
    const std::size_t n = a.size();
    BLST_ASSERT(n > 0);
    BLST_ASSERT(precomputed_b.algorithm == MiddleProductAlgorithm::BasicFFT);
    BLST_ASSERT(precomputed_b.n == n);
    BLST_ASSERT(out.size() >= n);

    const std::size_t fft_size = precomputed_b.fft_size;
    BLST_ASSERT(precomputed_b.first.size() == fft_size);

    std::span<const Fr> roots;
    std::span<const Fr> inv_roots;
    if (precomputed_b.roots.size() == fft_size &&
        precomputed_b.inv_roots.size() == fft_size) {
        roots = std::span<const Fr>(precomputed_b.roots);
        inv_roots = std::span<const Fr>(precomputed_b.inv_roots);
    }
    else {
        ScopedTimer timer(timings, "basic_pruned.roots_ns");
        workspace.roots = make_fft_roots(fft_size);
        workspace.inv_roots = make_ifft_roots(workspace.roots);
        roots = std::span<const Fr>(workspace.roots);
        inv_roots = std::span<const Fr>(workspace.inv_roots);
    }

    const std::size_t xn_a =
        (precomputed_b.xn_a > 0) ? precomputed_b.xn_a : fft_admissible_size(n, fft_size);
    const FftTruncatedForwardPlan* fft_a_plan = nullptr;
    const FftInverseTruncatedRangePlan* ifft_plan = nullptr;
    FftTruncatedForwardPlan local_fft_a_plan;
    FftInverseTruncatedRangePlan local_ifft_plan;
    if (precomputed_b.fft_a_plan.n == fft_size && precomputed_b.fft_a_plan.xn == xn_a) {
        fft_a_plan = &precomputed_b.fft_a_plan;
    }
    else {
        local_fft_a_plan = make_fft_truncated_forward_plan(fft_size, xn_a);
        fft_a_plan = &local_fft_a_plan;
    }
    if (precomputed_b.ifft_middle_plan.n == fft_size &&
        precomputed_b.ifft_middle_plan.start == n &&
        precomputed_b.ifft_middle_plan.length == n) {
        ifft_plan = &precomputed_b.ifft_middle_plan;
    }
    else {
        local_ifft_plan = make_fft_inverse_truncated_range_plan(fft_size, n, n);
        ifft_plan = &local_ifft_plan;
    }

    workspace.ensure_size(fft_size);
    std::span<G1Point> a_fft(workspace.a_fft);
    std::span<GT> c_fft(workspace.c_fft);
    std::span<const G2Point> b_fft(precomputed_b.first);
    std::fill(a_fft.begin(), a_fft.end(), G1Point::infinity());
    std::fill(c_fft.begin(), c_fft.end(), GT::identity());

    for (std::size_t i = 0; i < n; ++i) {
        a_fft[i] = a[i];
    }

    {
        ScopedTimer timer(timings, "basic_pruned.fft_a_ns");
        fft_forward_truncated_prefix<G1Point, Fr, CoeffCtx>(
            a_fft,
            roots,
            xn_a,
            *fft_a_plan,
            timings,
            "basic_pruned.fft_a");
    }

    {
        ScopedTimer timer(timings, "basic_pruned.pointwise_ns");
        for (std::size_t i = 0; i < fft_size; ++i) {
            c_fft[i] = CoeffCtx::mult(a_fft[i], b_fft[i]);
        }
    }

    {
        ScopedTimer timer(timings, "basic_pruned.ifft_ns");
        fft_inverse_truncated_range<GT, Fr, CoeffCtx>(
            c_fft,
            inv_roots,
            n,
            n,
            *ifft_plan,
            timings,
            "basic_pruned.ifft");
    }

    {
        ScopedTimer timer(timings, "basic_pruned.extract_ns");
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = c_fft[n + i];
        }
    }
}

inline void middle_product_basic_pruned(
    std::span<const G1Point> a,
    const MiddleProductPrecomputedB& precomputed_b,
    std::span<GT> out,
    TimingStats* timings = nullptr)
{
    MiddleProductBasicPrunedWorkspace workspace;
    middle_product_basic_pruned(a, precomputed_b, out, workspace, timings);
}

// Standalone specialized-flow variant with truncation-enabled substeps:
// - truncated forward FFT on A
// - truncated inverse FFT windows for both phase outputs
inline void middle_product_specialized_truncated(
    std::span<const G1Point> a,
    const MiddleProductPrecomputedB& precomputed_b,
    std::span<GT> out,
    TimingStats* timings = nullptr)
{
    ScopedTimer total_timer(timings, "middle_spec_trunc.total_ns");
    const std::size_t n = a.size();
    BLST_ASSERT(n > 0);
    BLST_ASSERT(precomputed_b.algorithm == MiddleProductAlgorithm::SpecializedFFT);
    BLST_ASSERT(precomputed_b.n == n);
    BLST_ASSERT(out.size() >= n);

    const std::size_t m = precomputed_b.fft_size;
    BLST_ASSERT(precomputed_b.first.size() == m);
    BLST_ASSERT(precomputed_b.second.size() == m);

    std::vector<Fr> roots;
    std::vector<Fr> inv_roots;
    {
        ScopedTimer timer(timings, "middle_spec_trunc.roots_ns");
        roots = make_fft_roots(m);
    }
    {
        ScopedTimer timer(timings, "middle_spec_trunc.inv_roots_ns");
        inv_roots = make_ifft_roots(roots);
    }

    std::vector<G1Point> a_fft(m, G1Point::infinity());
    std::vector<G2Point> b0_fft = precomputed_b.second;
    std::vector<G2Point> b1_fft = precomputed_b.first;
    std::vector<GT> scratch(m, GT::identity());
    std::vector<G1Point> a_rev_fft(m, G1Point::infinity());

    for (std::size_t i = 0; i < n; ++i) {
        a_fft[i] = a[i];
    }
    std::fill(out.begin(), out.begin() + static_cast<std::ptrdiff_t>(n), GT::identity());

    const std::size_t xn_a = fft_admissible_size(n, m);
    const FftTruncatedForwardPlan a_fwd_plan = make_fft_truncated_forward_plan(m, xn_a);
    const FftInverseTruncatedRangePlan phase1_ifft_plan =
        make_fft_inverse_truncated_range_plan(m, 0, n);
    const FftInverseTruncatedRangePlan phase2_ifft_plan =
        (n > 1) ? make_fft_inverse_truncated_range_plan(m, 0, n - 1)
                : make_fft_inverse_truncated_range_plan(m, 0, 1);

    // Phase 1: p = low_n(A * B1)
    {
        ScopedTimer timer(timings, "middle_spec_trunc.phase1_fft_a_ns");
        fft_forward_truncated_prefix<G1Point, Fr, CoeffCtx>(
            a_fft, roots, xn_a, a_fwd_plan, timings, "middle_spec_trunc.phase1_fft_a");
    }
    {
        ScopedTimer timer(timings, "middle_spec_trunc.phase1_pointwise_ns");
        for (std::size_t i = 0; i < m; ++i) {
            scratch[i] = CoeffCtx::mult(a_fft[i], b1_fft[i]);
        }
    }
    {
        ScopedTimer timer(timings, "middle_spec_trunc.phase1_ifft_ns");
        fft_inverse_truncated_range<GT, Fr, CoeffCtx>(
            scratch,
            inv_roots,
            0,
            n,
            phase1_ifft_plan,
            timings,
            "middle_spec_trunc.phase1_ifft");
    }
    {
        ScopedTimer timer(timings, "middle_spec_trunc.phase1_copy_ns");
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = scratch[i];
        }
    }

    // Phase 2: q = high_n(A * B0) via reversal trick.
    const Fr prim = primRootOfUnity<Fr>(m);
    {
        ScopedTimer timer(timings, "middle_spec_trunc.phase2_derive_ntt_ns");
        derive_ntt_of_reversed<G1Point, Fr, CoeffCtx>(
            std::span<const G1Point>(a_fft),
            std::span<G1Point>(a_rev_fft),
            n,
            m,
            prim,
            G1Point::infinity());
    }
    {
        ScopedTimer timer(timings, "middle_spec_trunc.phase2_pointwise_ns");
        for (std::size_t i = 0; i < m; ++i) {
            scratch[i] = CoeffCtx::mult(a_rev_fft[i], b0_fft[i]);
        }
    }
    if (n > 1) {
        ScopedTimer timer(timings, "middle_spec_trunc.phase2_ifft_ns");
        fft_inverse_truncated_range<GT, Fr, CoeffCtx>(
            scratch,
            inv_roots,
            0,
            n - 1,
            phase2_ifft_plan,
            timings,
            "middle_spec_trunc.phase2_ifft");
    }
    {
        ScopedTimer timer(timings, "middle_spec_trunc.phase2_combine_ns");
        for (std::size_t k = 0; k + 1 < n; ++k) {
            out[k] = CoeffCtx::plus(out[k], scratch[n - 2 - k]);
        }
    }
}

// TODO(BLST_CUDA): mirror middle_product() and precompute_middle_product_b() on GPU.

} // namespace blst

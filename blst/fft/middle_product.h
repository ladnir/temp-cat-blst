#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
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
        OneIfftCyclic,
    };

    enum class OneIfftPairingMode {
        Direct,
        PrecomputedLines,
    };

    struct alignas(64) OneIfftPreparedLineTable {
        Fp6 lines[68];
    };

    struct MiddleProductPrecomputedB {
        MiddleProductAlgorithm algorithm;
        std::size_t n;
        std::size_t fft_size;
        OneIfftPairingMode oneifft_pairing_mode = OneIfftPairingMode::Direct;
        std::vector<G2Point> first;
        std::vector<OneIfftPreparedLineTable> oneifft_precomputed_lines;
        std::vector<uint8_t> oneifft_d_is_zero;
        std::vector<Fr> roots;
        std::vector<Fr> inv_roots;
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
                    precomputed.first[j] = folded;
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
                    GT::precompute_lines(precomputed.oneifft_precomputed_lines[j].lines, q_aff);
                }
            }
            return precomputed;
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
                    d_fft[j] = G2CoeffCtx::plus(d_fft[j], zeta_b1);
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
                fft_inverse_truncated_range<GT, Fr, GTCoeffCtx>(
                    std::span<GT>(c_fft), std::span<const Fr>(inv_roots), n, n);
            }
            {
                ScopedTimer timer(timings, "oneifft.extract_ns");
                std::memcpy(out.data(), c_fft.data() + n, n * sizeof(GT));
            }
            return;
        }

        BLST_ASSERT(false);
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
            std::memset(a_fft.data(), 0, a_fft.size_bytes());
        }
        {
            ScopedTimer timer(timings, "oneifft_pre.fill_c_fft_ns");
            std::memset(c_fft.data(), 0, c_fft.size_bytes());
        }

        {

            ScopedTimer timer(timings, "oneifft_pre.copy_a_ns");
            memcpy(a_fft.data(), a.data(), n * sizeof(G1Point));
            //for (std::size_t i = 0; i < n; ++i) {
            //    a_fft[i] = a[i];
            //}
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
                    GT::miller_loop_lines(c_fft[i], precomputed_b.oneifft_precomputed_lines[i].lines, p_aff);
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
            fft_inverse_truncated_range<GT, Fr, GTCoeffCtx>(c_fft, inv_roots, n, n);
        }
        {
            ScopedTimer timer(timings, "oneifft_pre.extract_ns");
            std::memcpy(out.data(), c_fft.data() + n, n * sizeof(GT));
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
        BLST_ASSERT(false);
    }

    // TODO(BLST_CUDA): mirror middle_product() and precompute_middle_product_b() on GPU.

} // namespace blst

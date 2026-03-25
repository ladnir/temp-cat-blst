#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
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

    // Middle product for G1 x G2 -> GT:
    // - input A has n coefficients in G1
    // - input B has 2n coefficients in G2
    // - output is the middle n coefficients of the full convolution
    //
    // This file intentionally exposes only two production paths:
    // 1. middle_product_basic: a simple reference implementation that recomputes
    //    everything each call.
    // 2. middle_product: the reduced-iFFT implementation used for the faster path.
    //
    // The reduced-iFFT path also supports a precomputed-B form for repeated calls
    // with the same G2 input.

    enum class MiddleProductPairingMode {
        Direct,
        PrecomputedLines,
    };

    // Stores precomputed Miller-loop line data for one G2 point in the reduced
    // middle-product path.
    struct alignas(64) MiddleProductPreparedLineTable {
        Fp6 lines[68];
    };

    // Holds the fixed-B precomputation for the reduced middle-product algorithm.
    //
    // Invariants:
    // - n > 0
    // - fft_size is a power of two and fft_size >= 2n
    // - b_fft.size() == fft_size
    // - roots.size() == fft_size
    // - inv_roots.size() == fft_size
    // - when pairing_mode == PrecomputedLines, precomputed_lines.size() is
    //   either 0 or fft_size
    struct MiddleProductPrecomputedB {
        std::size_t n = 0;
        std::size_t fft_size = 0;
        MiddleProductPairingMode pairing_mode = MiddleProductPairingMode::Direct;
        std::vector<G2Point> b_fft;
        std::vector<MiddleProductPreparedLineTable> precomputed_lines;
        std::vector<uint8_t> b_is_zero;
        std::vector<Fr> roots;
        std::vector<Fr> inv_roots;
    };

    // Scratch buffers reused across calls to the reduced middle-product path.
    //
    // Invariants after ensure_size(m):
    // - a_fft.size() == m
    // - c_fft.size() == m
    struct MiddleProductWorkspace {
        std::vector<G1Point> a_fft;
        std::vector<GT> c_fft;

        // Resize scratch storage to exactly m coefficients.
        //
        // @param m FFT domain size. Must be the fft_size used by the associated
        //   middle-product precomputation.
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

    // Return the FFT size used by the simple reference implementation.
    //
    // The basic path computes the full convolution and extracts the middle n
    // coefficients, so it needs an FFT size large enough for 3n - 1 terms.
    //
    // @param n Number of G1 coefficients in A. Must satisfy n > 0.
    // @return Smallest power of two >= 3n - 1.
    //
    // Invariants:
    // - return value is a power of two
    // - return value >= 3n - 1
    inline std::size_t middle_product_basic_fft_size(std::size_t n)
    {
        BLST_ASSERT(n > 0);

        std::size_t fft_size = 1;
        const std::size_t target = (3 * n) - 1;
        while (fft_size < target) {
            fft_size <<= 1;
        }
        return fft_size;
    }

    // Return the FFT size used by the reduced middle-product algorithm.
    //
    // The reduced path works in a 2n-sized cyclic domain and only computes the
    // output window that is actually needed.
    //
    // @param n Number of G1 coefficients in A. Must satisfy n > 0.
    // @return Smallest power of two >= 2n.
    //
    // Invariants:
    // - return value is a power of two
    // - return value >= 2n
    inline std::size_t middle_product_fft_size(std::size_t n)
    {
        BLST_ASSERT(n > 0);

        std::size_t fft_size = 1;
        const std::size_t target = 2 * n;
        while (fft_size < target) {
            fft_size <<= 1;
        }
        return fft_size;
    }

    // Round an active prefix length up to the truncated-FFT schedule used by the
    // implementation.
    //
    // @param x Requested active prefix length.
    // @param fft_size FFT domain size. Must be a non-zero power of two.
    // @return Rounded prefix length, clamped to [1, fft_size].
    //
    // Invariants:
    // - return value >= 1
    // - return value <= fft_size
    // - return value is a multiple of 16 unless it is exactly fft_size
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

    // Construct the forward roots of unity for an n-point FFT domain.
    //
    // @param n FFT domain size. Must be a supported power of two in Fr.
    // @return Vector of size n containing successive powers of a primitive
    //   n-th root of unity, starting at 1.
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

    // Construct inverse roots corresponding to a forward root table.
    //
    // @param roots Forward roots table.
    // @return Vector of size roots.size() where each entry is the multiplicative
    //   inverse of the corresponding forward root.
    inline std::vector<Fr> make_ifft_roots(std::span<const Fr> roots)
    {
        std::vector<Fr> inv_roots(roots.size());
        for (std::size_t i = 0; i < roots.size(); ++i) {
            Fr::reciprocal(inv_roots[i], roots[i]);
        }
        return inv_roots;
    }

    // Precompute the fixed 2n-sized G2 input for the reduced-iFFT middle-product
    // path. This is intended for repeated calls with the same B input.
    //
    // @param b Fixed G2 coefficients. Must contain exactly 2n points.
    // @param n Number of output coefficients, and the number of G1 coefficients
    //   expected by the corresponding middle_product(...) calls.
    // @param timings Optional timing sink. May be null.
    // @param pairing_mode Selects direct pairings or precomputed Miller lines
    //   during the pointwise multiplication step.
    // @return Precomputed state for repeated reduced middle-product evaluations.
    //
    // Invariants:
    // - n > 0
    // - b.size() == 2 * n
    // - returned state satisfies the invariants documented on
    //   MiddleProductPrecomputedB
    inline MiddleProductPrecomputedB precompute_middle_product_b(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr,
        MiddleProductPairingMode pairing_mode = MiddleProductPairingMode::Direct)
    {
        ScopedTimer total_timer(timings, "middle_product_precompute.total_ns");
        BLST_ASSERT(n > 0);
        BLST_ASSERT(b.size() == 2 * n);

        MiddleProductPrecomputedB precomputed;
        precomputed.n = n;
        precomputed.fft_size = middle_product_fft_size(n);
        precomputed.pairing_mode = pairing_mode;
        precomputed.b_fft.assign(precomputed.fft_size, G2Point::infinity());

        {
            ScopedTimer timer(timings, "middle_product_precompute.roots_ns");
            precomputed.roots = make_fft_roots(precomputed.fft_size);
        }
        precomputed.inv_roots = make_ifft_roots(precomputed.roots);

        std::vector<G2Point> b0_fft(precomputed.fft_size, G2Point::infinity());
        std::vector<G2Point> b1_fft(precomputed.fft_size, G2Point::infinity());
        for (std::size_t i = 0; i < n; ++i) {
            b0_fft[i] = b[i];
            b1_fft[i] = b[n + i];
        }

        const std::size_t xn_b = fft_admissible_size(n, precomputed.fft_size);
        {
            ScopedTimer timer(timings, "middle_product_precompute.fft_b0_ns");
            fft_forward_truncated_prefix<G2Point, Fr, G2CoeffCtx>(
                b0_fft,
                std::span<const Fr>(precomputed.roots),
                xn_b,
                timings,
                "middle_product_precompute.fft_b0");
        }
        {
            ScopedTimer timer(timings, "middle_product_precompute.fft_b1_ns");
            fft_forward_truncated_prefix<G2Point, Fr, G2CoeffCtx>(
                b1_fft,
                std::span<const Fr>(precomputed.roots),
                xn_b,
                timings,
                "middle_product_precompute.fft_b1");
        }
        {
            ScopedTimer timer(timings, "middle_product_precompute.fold_ns");
            for (std::size_t j = 0; j < precomputed.fft_size; ++j) {
                const std::size_t idx = (j * (n % precomputed.fft_size)) % precomputed.fft_size;
                const Fr zeta = precomputed.roots[idx];
                const G2Point zeta_b1 = G2CoeffCtx::mult(b1_fft[j], zeta);
                precomputed.b_fft[j] = G2CoeffCtx::plus(b0_fft[j], zeta_b1);
            }
        }

        if (pairing_mode == MiddleProductPairingMode::PrecomputedLines) {
            ScopedTimer timer(timings, "middle_product_precompute.lines_ns");
            precomputed.precomputed_lines.resize(precomputed.fft_size);
            precomputed.b_is_zero.resize(precomputed.fft_size);
            for (std::size_t j = 0; j < precomputed.fft_size; ++j) {
                const bool is_zero = (G2Point::is_inf(precomputed.b_fft[j]) != 0);
                precomputed.b_is_zero[j] = static_cast<uint8_t>(is_zero ? 1 : 0);
                if (is_zero) {
                    continue;
                }
                G2AffinePoint q_aff;
                G2Point::to_affine(q_aff, precomputed.b_fft[j]);
                GT::precompute_lines(precomputed.precomputed_lines[j].lines, q_aff);
            }
        }

        return precomputed;
    }

    // Reference implementation. This keeps the full basic FFT flow in one place
    // and intentionally does not offer a precomputed-B variant.
    //
    // @param a G1 coefficients of length n.
    // @param b G2 coefficients of length 2n.
    // @param out Output buffer for the middle n GT coefficients.
    // @param timings Optional timing sink. May be null.
    //
    // Invariants:
    // - a.size() > 0
    // - b.size() == 2 * a.size()
    // - out.size() >= a.size()
    inline void middle_product_basic(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        ScopedTimer total_timer(timings, "middle_product_basic.total_ns");
        const std::size_t n = a.size();
        BLST_ASSERT(n > 0);
        BLST_ASSERT(b.size() == 2 * n);
        BLST_ASSERT(out.size() >= n);

        const std::size_t fft_size = middle_product_basic_fft_size(n);
        std::vector<Fr> roots;
        std::vector<Fr> inv_roots;
        {
            ScopedTimer timer(timings, "middle_product_basic.roots_ns");
            roots = make_fft_roots(fft_size);
        }
        {
            ScopedTimer timer(timings, "middle_product_basic.inv_roots_ns");
            inv_roots = make_ifft_roots(roots);
        }

        std::vector<G1Point> a_fft(fft_size, G1Point::infinity());
        std::vector<G2Point> b_fft(fft_size, G2Point::infinity());
        std::vector<GT> c_fft(fft_size, GT::identity());

        for (std::size_t i = 0; i < n; ++i) {
            a_fft[i] = a[i];
        }
        for (std::size_t i = 0; i < 2 * n; ++i) {
            b_fft[i] = b[i];
        }

        const std::size_t xn_a = fft_admissible_size(n, fft_size);
        const std::size_t xn_b = fft_admissible_size(2 * n, fft_size);
        poly_mult_fft<G1Point, G2Point, GT, Fr, CoeffCtx>(
            a_fft,
            b_fft,
            false,
            c_fft,
            std::span<const Fr>(roots),
            std::span<const Fr>(inv_roots),
            timings,
            xn_a,
            xn_b);

        {
            ScopedTimer timer(timings, "middle_product_basic.extract_ns");
            for (std::size_t i = 0; i < n; ++i) {
                out[i] = c_fft[n + i];
            }
        }
    }

    // Reduced middle-product overload for a fixed, precomputed B input with an
    // explicit reusable workspace.
    //
    // @param a G1 coefficients of length n.
    // @param precomputed_b Fixed-B precomputation produced by
    //   precompute_middle_product_b(...).
    // @param out Output buffer for the middle n GT coefficients.
    // @param workspace Scratch storage reused across calls.
    // @param timings Optional timing sink. May be null.
    //
    // Invariants:
    // - a.size() > 0
    // - precomputed_b.n == a.size()
    // - out.size() >= a.size()
    // - precomputed_b satisfies the invariants documented on
    //   MiddleProductPrecomputedB
    inline void middle_product(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedB& precomputed_b,
        std::span<GT> out,
        MiddleProductWorkspace& workspace,
        TimingStats* timings = nullptr);

    // Reduced middle-product overload that accepts the raw 2n-sized G2 input.
    //
    // This is the production fast path without explicit precompute reuse. It
    // internally precomputes the fixed-B state and then runs the reduced-iFFT
    // algorithm.
    //
    // @param a G1 coefficients of length n.
    // @param b G2 coefficients of length 2n.
    // @param out Output buffer for the middle n GT coefficients.
    // @param timings Optional timing sink. May be null.
    //
    // Invariants:
    // - a.size() > 0
    // - b.size() == 2 * a.size()
    // - out.size() >= a.size()
    inline void middle_product(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        ScopedTimer total_timer(timings, "middle_product.total_ns");
        const std::size_t n = a.size();
        BLST_ASSERT(n > 0);
        BLST_ASSERT(b.size() == 2 * n);
        BLST_ASSERT(out.size() >= n);

        const MiddleProductPrecomputedB precomputed = precompute_middle_product_b(b, n, timings);
        MiddleProductWorkspace workspace;
        middle_product(a, precomputed, out, workspace, timings);
    }

    // Reduced middle-product overload for a fixed, precomputed B input with an
    // explicit reusable workspace.
    //
    // This is the lowest-level production entry point for repeated calls with
    // the same G2 input.
    //
    // @param a G1 coefficients of length n.
    // @param precomputed_b Fixed-B precomputation produced by
    //   precompute_middle_product_b(...).
    // @param out Output buffer for the middle n GT coefficients.
    // @param workspace Scratch storage reused across calls.
    // @param timings Optional timing sink. May be null.
    //
    // Invariants:
    // - a.size() > 0
    // - precomputed_b.n == a.size()
    // - out.size() >= a.size()
    // - precomputed_b.b_fft.size() == precomputed_b.fft_size
    // - precomputed_b.roots.size() == precomputed_b.fft_size
    // - precomputed_b.inv_roots.size() == precomputed_b.fft_size
    inline void middle_product(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedB& precomputed_b,
        std::span<GT> out,
        MiddleProductWorkspace& workspace,
        TimingStats* timings)
    {
        ScopedTimer total_timer(timings, "middle_product_pre.total_ns");
        const std::size_t n = a.size();
        BLST_ASSERT(n > 0);
        BLST_ASSERT(precomputed_b.n == n);
        BLST_ASSERT(out.size() >= n);
        BLST_ASSERT(precomputed_b.b_fft.size() == precomputed_b.fft_size);
        BLST_ASSERT(precomputed_b.roots.size() == precomputed_b.fft_size);
        BLST_ASSERT(precomputed_b.inv_roots.size() == precomputed_b.fft_size);

        workspace.ensure_size(precomputed_b.fft_size);
        std::span<G1Point> a_fft(workspace.a_fft);
        std::span<GT> c_fft(workspace.c_fft);

        {
            ScopedTimer timer(timings, "middle_product_pre.fill_a_fft_ns");
            std::memset(a_fft.data(), 0, a_fft.size_bytes());
        }
        {
            ScopedTimer timer(timings, "middle_product_pre.fill_c_fft_ns");
            std::memset(c_fft.data(), 0, c_fft.size_bytes());
        }
        {
            ScopedTimer timer(timings, "middle_product_pre.copy_a_ns");
            std::memcpy(a_fft.data(), a.data(), n * sizeof(G1Point));
        }

        const std::size_t xn_a = fft_admissible_size(n, precomputed_b.fft_size);
        {
            ScopedTimer timer(timings, "middle_product_pre.fft_a_ns");
            fft_forward_truncated_prefix<G1Point, Fr, G1CoeffCtx>(
                a_fft,
                std::span<const Fr>(precomputed_b.roots),
                xn_a,
                timings,
                "middle_product_pre.fft_a");
        }
        {
            ScopedTimer timer(timings, "middle_product_pre.pointwise_ns");
            if (precomputed_b.pairing_mode == MiddleProductPairingMode::PrecomputedLines &&
                precomputed_b.precomputed_lines.size() == precomputed_b.fft_size) {
                const bool has_zero_mask = (precomputed_b.b_is_zero.size() == precomputed_b.fft_size);
                for (std::size_t i = 0; i < precomputed_b.fft_size; ++i) {
                    if (G1Point::is_inf(a_fft[i]) != 0 ||
                        (has_zero_mask && precomputed_b.b_is_zero[i] != 0)) {
                        c_fft[i] = GT::identity();
                        continue;
                    }
                    G1AffinePoint p_aff;
                    G1Point::to_affine(p_aff, a_fft[i]);
                    GT::miller_loop_lines(c_fft[i], precomputed_b.precomputed_lines[i].lines, p_aff);
                    GT::final_mult(c_fft[i], c_fft[i]);
                }
            }
            else {
                for (std::size_t i = 0; i < precomputed_b.fft_size; ++i) {
                    c_fft[i] = CoeffCtx::mult(a_fft[i], precomputed_b.b_fft[i]);
                }
            }
        }
        {
            ScopedTimer timer(timings, "middle_product_pre.ifft_ns");
            fft_inverse_truncated_range<GT, Fr, GTCoeffCtx>(
                c_fft,
                std::span<const Fr>(precomputed_b.inv_roots),
                n,
                n);
        }
        {
            ScopedTimer timer(timings, "middle_product_pre.extract_ns");
            std::memcpy(out.data(), c_fft.data() + n, n * sizeof(GT));
        }
    }

    // Reduced middle-product overload for a fixed, precomputed B input using an
    // internal temporary workspace.
    //
    // @param a G1 coefficients of length n.
    // @param precomputed_b Fixed-B precomputation produced by
    //   precompute_middle_product_b(...).
    // @param out Output buffer for the middle n GT coefficients.
    // @param timings Optional timing sink. May be null.
    //
    // Invariants:
    // - a.size() > 0
    // - precomputed_b.n == a.size()
    // - out.size() >= a.size()
    inline void middle_product(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedB& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        MiddleProductWorkspace workspace;
        middle_product(a, precomputed_b, out, workspace, timings);
    }

    // TODO(BLST_CUDA): mirror middle_product_basic(), middle_product(), and
    // precompute_middle_product_b() on GPU.

} // namespace blst

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>

#include "blst/avx/curve/GT.h"
#include "blst/avx/curve/map_to_G.h"
#include "blst/avx/fft/gt_context.h"
#include "blst/common/defines.h"
#include "blst/common/timer.h"
#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/primitveRoot.h"
#include "blst/fields/Fr.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/point.h"

namespace blst {

    enum class MiddleProductPairingMode {
        Direct,
        PrecomputedLines,
    };

    struct alignas(64) MiddleProductPreparedScalarLineTable {
        Fp6 lines[68];
    };

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

    namespace detail {

        struct BasicMiddleProductCoeffCtx {
            template<typename G, typename F>
            static auto plus(const G& a, const F& b)
            {
                G r;
                G::dadd(r, a, b);
                return r;
            }

            template<typename G, typename F>
            static auto minus(const G& a, const F& b)
            {
                G r;
                F neg_b;
                F::neg(neg_b, b);
                G::dadd(r, a, neg_b);
                return r;
            }

            template<typename G, typename F>
            static auto mult(const G& a, const F& b)
            {
                if constexpr (std::is_same_v<G, F>) {
                    G r;
                    G::mul(r, a, b);
                    return r;
                } else {
                    G r;
                    F b_non_mont = {};
                    F::from_mont(b_non_mont, b);
                    G::mult(r, a, reinterpret_cast<const byte*>(&b_non_mont), 256);
                    return r;
                }
            }

            static auto mult(const G1Point& a, const G2Point& b)
            {
                G1AffinePoint p_aff = {};
                G1Point::to_affine(p_aff, a);
                G2AffinePoint q_aff = {};
                G2Point::to_affine(q_aff, b);

                GT pair = GT::identity();
                GT::miller_loop(pair, &q_aff, &p_aff);
                GT::final_mult(pair, pair);
                return pair;
            }

            static auto inverse(const Fr& a)
            {
                Fr r = {};
                Fr::reciprocal(r, a);
                return r;
            }

            template<typename F>
            static auto from_u64(uint64_t a)
            {
                Fr r = {};
                Fr::from_uint64(r, a);
                return r;
            }
        };

        template<class Group>
        struct GroupFftCoeffCtx {
            static auto plus(const Group& a, const Group& b)
            {
                Group r;
                Group::dadd(r, a, b);
                return r;
            }

            static auto minus(const Group& a, const Group& b)
            {
                Group neg_b;
                Group::neg(neg_b, b);
                Group r;
                Group::dadd(r, a, neg_b);
                return r;
            }

            static auto mult(const Group& a, const Fr& b)
            {
                Group r;
                Fr b_non_mont = {};
                Fr::from_mont(b_non_mont, b);
                Group::mult(r, a, reinterpret_cast<const byte*>(&b_non_mont), 256);
                return r;
            }

            static auto inverse(const Fr& a)
            {
                Fr r = {};
                Fr::reciprocal(r, a);
                return r;
            }

            template<typename F>
            static auto from_u64(uint64_t a)
            {
                Fr r = {};
                Fr::from_uint64(r, a);
                return r;
            }
        };

        struct ScalarMiddleProductContext {
            using G1Point = scalar::G1Point;
            using G2Point = scalar::G2Point;
            using G1AffinePoint = scalar::G1AffinePoint;
            using G2AffinePoint = scalar::G2AffinePoint;
            using G1CoeffCtx = GroupFftCoeffCtx<G1Point>;
            using G2CoeffCtx = GroupFftCoeffCtx<G2Point>;
            using GtFftElement = scalar::GT;
            using GtCoeffCtx = GroupFftCoeffCtx<GtFftElement>;
            using PreparedInverseRoots = std::vector<Fr>;
            using PreparedLineTable = MiddleProductPreparedScalarLineTable;

            static G1Point g1_infinity()
            {
                return G1Point::infinity();
            }

            static G2Point g2_infinity()
            {
                return G2Point::infinity();
            }

            static G1Point to_backend_g1(const blst::G1Point& point)
            {
#if BLST_DEFAULT_AVX
                return avx::avx_map_to_G_detail::to_scalar(point);
#else
                return point;
#endif
            }

            static G2Point to_backend_g2(const blst::G2Point& point)
            {
#if BLST_DEFAULT_AVX
                return avx::avx_map_to_G_detail::to_scalar(point);
#else
                return point;
#endif
            }

            static GT to_public_gt(const GtFftElement& value)
            {
#if BLST_DEFAULT_AVX
                return avx::AvxGt52::to_default(avx::AvxGt52::from_scalar(value));
#else
                return value;
#endif
            }

            static PreparedInverseRoots make_inverse_roots(std::span<const Fr> roots)
            {
                return make_ifft_roots(roots);
            }

            static std::size_t inverse_roots_size(const PreparedInverseRoots& roots)
            {
                return roots.size();
            }

            static void prepare_lines(PreparedLineTable& out, const G2Point& point)
            {
                G2AffinePoint point_affine = {};
                G2Point::to_affine(point_affine, point);
                scalar::GT::precompute_lines(out.lines, point_affine);
            }

            static GtFftElement pointwise_pair_direct(const G1Point& left, const G2Point& right)
            {
                G1AffinePoint p_aff = {};
                G1Point::to_affine(p_aff, left);
                G2AffinePoint q_aff = {};
                G2Point::to_affine(q_aff, right);

                GtFftElement pair = GtFftElement::identity();
                GtFftElement::miller_loop(pair, &q_aff, &p_aff);
                GtFftElement::final_mult(pair, pair);
                return pair;
            }

            static GtFftElement pointwise_pair_lines(
                const PreparedLineTable& lines,
                const G1Point& point)
            {
                G1AffinePoint point_affine = {};
                G1Point::to_affine(point_affine, point);

                GtFftElement pair = GtFftElement::identity();
                GtFftElement::miller_loop_lines(pair, lines.lines, point_affine);
                GtFftElement::final_mult(pair, pair);
                return pair;
            }

            static GtFftElement identity()
            {
                return GtFftElement::identity();
            }

            static void inverse_fft_middle_window(
                std::span<GtFftElement> values,
                const PreparedInverseRoots& inv_roots,
                std::size_t start,
                std::size_t length)
            {
                fft_inverse_truncated_range<GtFftElement, Fr, GtCoeffCtx>(
                    values,
                    std::span<const Fr>(inv_roots),
                    start,
                    length);
            }

            static void extract_middle_window(
                std::span<GT> out,
                std::span<const GtFftElement> values,
                std::size_t start,
                std::size_t length)
            {
                for (std::size_t i = 0; i < length; ++i) {
                    out[i] = to_public_gt(values[start + i]);
                }
            }
        };

        struct AvxMiddleProductContext {
            using G1Point = avx::G1Point;
            using G2Point = avx::G2Point;
            using G1AffinePoint = avx::G1AffinePoint;
            using G2AffinePoint = avx::G2AffinePoint;
            using G1CoeffCtx = GroupFftCoeffCtx<G1Point>;
            using G2CoeffCtx = GroupFftCoeffCtx<G2Point>;
            using GtFftElement = avx::AvxGt52;
            using PreparedInverseRoots = avx::AvxGt52PreparedRoots;
            using PreparedLineTable = avx::PreparedLineTable;

            static G1Point g1_infinity()
            {
                return G1Point::infinity();
            }

            static G2Point g2_infinity()
            {
                return G2Point::infinity();
            }

            static G1Point to_backend_g1(const blst::G1Point& point)
            {
#if BLST_DEFAULT_AVX
                return point;
#else
                return avx::avx_map_to_G_detail::from_scalar(point);
#endif
            }

            static G2Point to_backend_g2(const blst::G2Point& point)
            {
#if BLST_DEFAULT_AVX
                return point;
#else
                return avx::avx_map_to_G_detail::from_scalar(point);
#endif
            }

            static GT to_public_gt(const GtFftElement& value)
            {
                return avx::AvxGt52::to_default(value);
            }

            static PreparedInverseRoots make_inverse_roots(std::span<const Fr> roots)
            {
                const std::vector<Fr> inv_roots = make_ifft_roots(roots);
                PreparedInverseRoots prepared{std::span<const Fr>(inv_roots)};
                prepared.set_scale_from_size(roots.size());
                return prepared;
            }

            static std::size_t inverse_roots_size(const PreparedInverseRoots& roots)
            {
                return roots.size();
            }

            static void prepare_lines(PreparedLineTable& out, const G2Point& point)
            {
                G2AffinePoint point_affine = {};
                G2Point::to_affine(point_affine, point);
                avx::GT::precompute_lines(out.lines, point_affine);
            }

            static GtFftElement pointwise_pair_direct(const G1Point& left, const G2Point& right)
            {
                G1AffinePoint p_aff = {};
                G1Point::to_affine(p_aff, left);
                G2AffinePoint q_aff = {};
                G2Point::to_affine(q_aff, right);

                avx::GT pair = avx::GT::identity();
                avx::GT::miller_loop(pair, &q_aff, &p_aff);
                avx::GT::final_mult(pair, pair);
                return GtFftElement(pair.data);
            }

            static GtFftElement pointwise_pair_lines(
                const PreparedLineTable& lines,
                const G1Point& point)
            {
                G1AffinePoint point_affine = {};
                G1Point::to_affine(point_affine, point);

                avx::GT pair = avx::GT::identity();
                avx::GT::miller_loop_lines(pair, lines.lines, point_affine);
                avx::GT::final_mult(pair, pair);
                return GtFftElement(pair.data);
            }

            static GtFftElement identity()
            {
                return GtFftElement::identity();
            }

            static void inverse_fft_middle_window(
                std::span<GtFftElement> values,
                const PreparedInverseRoots& inv_roots,
                std::size_t start,
                std::size_t length)
            {
                avx::AvxGt52CoeffCtx::ScopedRoots scoped_roots(inv_roots);
                fft_inverse_truncated_range<GtFftElement, Fr, avx::AvxGt52CoeffCtx>(
                    values,
                    inv_roots.scalar_span(),
                    start,
                    length);
            }

            static void extract_middle_window(
                std::span<GT> out,
                std::span<const GtFftElement> values,
                std::size_t start,
                std::size_t length)
            {
                for (std::size_t i = 0; i < length; ++i) {
                    out[i] = to_public_gt(values[start + i]);
                }
            }
        };

#if BLST_DEFAULT_AVX
        using DefaultMiddleProductContext = AvxMiddleProductContext;
#else
        using DefaultMiddleProductContext = ScalarMiddleProductContext;
#endif

        template<class Context>
        struct MiddleProductPrecomputedBImpl {
            std::size_t n = 0;
            std::size_t fft_size = 0;
            MiddleProductPairingMode pairing_mode = MiddleProductPairingMode::Direct;
            std::vector<typename Context::G2Point> b_tail;
            std::vector<typename Context::G2Point> b_fft;
            std::vector<typename Context::PreparedLineTable> tail_precomputed_lines;
            std::vector<uint8_t> tail_b_is_zero;
            std::vector<typename Context::PreparedLineTable> precomputed_lines;
            std::vector<uint8_t> b_is_zero;
            std::vector<Fr> roots;
            typename Context::PreparedInverseRoots inv_roots;
        };

        template<class Context>
        struct MiddleProductWorkspaceImpl {
            std::vector<typename Context::G1Point> a_fft;
            std::vector<typename Context::GtFftElement> c_fft;

            void ensure_size(std::size_t m)
            {
                if (a_fft.size() != m) {
                    a_fft.assign(m, Context::g1_infinity());
                }
                if (c_fft.size() != m) {
                    c_fft.assign(m, Context::identity());
                }
            }
        };

        template<class Context>
        inline std::vector<typename Context::G1Point> convert_g1_inputs(
            std::span<const blst::G1Point> input)
        {
            std::vector<typename Context::G1Point> out;
            out.reserve(input.size());
            for (const blst::G1Point& point : input) {
                out.push_back(Context::to_backend_g1(point));
            }
            return out;
        }

        template<class Context>
        inline std::vector<typename Context::G2Point> convert_g2_inputs(
            std::span<const blst::G2Point> input)
        {
            std::vector<typename Context::G2Point> out;
            out.reserve(input.size());
            for (const blst::G2Point& point : input) {
                out.push_back(Context::to_backend_g2(point));
            }
            return out;
        }

        template<class Context>
        inline MiddleProductPrecomputedBImpl<Context> precompute_middle_product_b_impl(
            std::span<const typename Context::G2Point> b,
            std::size_t n,
            TimingStats* timings,
            MiddleProductPairingMode pairing_mode)
        {
            ScopedTimer total_timer(timings, "middle_product_precompute.total_ns");
            BLST_ASSERT(n > 0);
            BLST_ASSERT(b.size() == 2 * n);

            MiddleProductPrecomputedBImpl<Context> precomputed;
            precomputed.n = n;
            precomputed.fft_size = middle_product_fft_size(n);
            precomputed.pairing_mode = pairing_mode;
            precomputed.b_tail.assign(n, Context::g2_infinity());
            precomputed.b_fft.assign(precomputed.fft_size, Context::g2_infinity());

            {
                ScopedTimer timer(timings, "middle_product_precompute.roots_ns");
                precomputed.roots = make_fft_roots(precomputed.fft_size);
            }
            precomputed.inv_roots = Context::make_inverse_roots(precomputed.roots);

            std::vector<typename Context::G2Point> b0_fft(precomputed.fft_size, Context::g2_infinity());
            std::vector<typename Context::G2Point> b1_fft(precomputed.fft_size, Context::g2_infinity());
            for (std::size_t i = 0; i < n; ++i) {
                b0_fft[i] = b[i];
                precomputed.b_tail[i] = b[n + i];
                b1_fft[i] = precomputed.b_tail[i];
            }

            const std::size_t xn_b = fft_admissible_size(n, precomputed.fft_size);
            {
                ScopedTimer timer(timings, "middle_product_precompute.fft_b0_ns");
                fft_forward_truncated_prefix<typename Context::G2Point, Fr, typename Context::G2CoeffCtx>(
                    b0_fft,
                    std::span<const Fr>(precomputed.roots),
                    xn_b,
                    timings,
                    "middle_product_precompute.fft_b0");
            }
            {
                ScopedTimer timer(timings, "middle_product_precompute.fft_b1_ns");
                fft_forward_truncated_prefix<typename Context::G2Point, Fr, typename Context::G2CoeffCtx>(
                    b1_fft,
                    std::span<const Fr>(precomputed.roots),
                    xn_b,
                    timings,
                    "middle_product_precompute.fft_b1");
            }
            {
                ScopedTimer timer(timings, "middle_product_precompute.fold_ns");
                for (std::size_t j = 0; j < precomputed.fft_size; ++j) {
                    const std::size_t idx =
                        (j * (n % precomputed.fft_size)) % precomputed.fft_size;
                    const Fr& zeta = precomputed.roots[idx];
                    const auto zeta_b1 = Context::G2CoeffCtx::mult(b1_fft[j], zeta);
                    precomputed.b_fft[j] = Context::G2CoeffCtx::plus(b0_fft[j], zeta_b1);
                }
            }

            if (pairing_mode == MiddleProductPairingMode::PrecomputedLines) {
                ScopedTimer timer(timings, "middle_product_precompute.lines_ns");
                precomputed.tail_precomputed_lines.resize(n);
                precomputed.tail_b_is_zero.resize(n);
                for (std::size_t i = 0; i < n; ++i) {
                    const bool is_zero = (Context::G2Point::is_inf(precomputed.b_tail[i]) != 0);
                    precomputed.tail_b_is_zero[i] = static_cast<uint8_t>(is_zero ? 1 : 0);
                    if (is_zero) {
                        continue;
                    }
                    Context::prepare_lines(precomputed.tail_precomputed_lines[i], precomputed.b_tail[i]);
                }
                precomputed.precomputed_lines.resize(precomputed.fft_size);
                precomputed.b_is_zero.resize(precomputed.fft_size);
                for (std::size_t j = 0; j < precomputed.fft_size; ++j) {
                    const bool is_zero = (Context::G2Point::is_inf(precomputed.b_fft[j]) != 0);
                    precomputed.b_is_zero[j] = static_cast<uint8_t>(is_zero ? 1 : 0);
                    if (is_zero) {
                        continue;
                    }
                    Context::prepare_lines(precomputed.precomputed_lines[j], precomputed.b_fft[j]);
                }
            }

            return precomputed;
        }

        template<class Context>
        inline void middle_product_impl(
            std::span<const typename Context::G1Point> a,
            const MiddleProductPrecomputedBImpl<Context>& precomputed_b,
            std::span<GT> out,
            MiddleProductWorkspaceImpl<Context>& workspace,
            TimingStats* timings)
        {
            ScopedTimer total_timer(timings, "middle_product_pre.total_ns");
            const std::size_t n = a.size();
            BLST_ASSERT(n > 0);
            BLST_ASSERT(precomputed_b.n == n);
            BLST_ASSERT(out.size() >= n);
            BLST_ASSERT(precomputed_b.b_fft.size() == precomputed_b.fft_size);
            BLST_ASSERT(precomputed_b.roots.size() == precomputed_b.fft_size);
            BLST_ASSERT(Context::inverse_roots_size(precomputed_b.inv_roots) == precomputed_b.fft_size);

            workspace.ensure_size(precomputed_b.fft_size);
            std::span<typename Context::G1Point> a_fft(workspace.a_fft);
            std::span<typename Context::GtFftElement> c_fft(workspace.c_fft);

            {
                ScopedTimer timer(timings, "middle_product_pre.fill_a_fft_ns");
                std::fill(a_fft.begin(), a_fft.end(), Context::g1_infinity());
            }
            {
                ScopedTimer timer(timings, "middle_product_pre.copy_a_ns");
                std::copy(a.begin(), a.end(), a_fft.begin());
            }
            {
                ScopedTimer timer(timings, "middle_product_pre.fill_c_fft_ns");
            }

            const std::size_t xn_a = fft_admissible_size(n, precomputed_b.fft_size);
            {
                ScopedTimer timer(timings, "middle_product_pre.fft_a_ns");
                fft_forward_truncated_prefix<typename Context::G1Point, Fr, typename Context::G1CoeffCtx>(
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
                    const bool has_zero_mask =
                        (precomputed_b.b_is_zero.size() == precomputed_b.fft_size);
                    for (std::size_t i = 0; i < precomputed_b.fft_size; ++i) {
                        if (Context::G1Point::is_inf(a_fft[i]) != 0 ||
                            (has_zero_mask && precomputed_b.b_is_zero[i] != 0)) {
                            c_fft[i] = Context::identity();
                            continue;
                        }
                        c_fft[i] = Context::pointwise_pair_lines(
                            precomputed_b.precomputed_lines[i],
                            a_fft[i]);
                    }
                }
                else {
                    for (std::size_t i = 0; i < precomputed_b.fft_size; ++i) {
                        c_fft[i] = Context::pointwise_pair_direct(a_fft[i], precomputed_b.b_fft[i]);
                    }
                }
            }
            {
                ScopedTimer timer(timings, "middle_product_pre.ifft_ns");
                Context::inverse_fft_middle_window(c_fft, precomputed_b.inv_roots, n, n);
            }
            {
                ScopedTimer timer(timings, "middle_product_pre.extract_ns");
                Context::extract_middle_window(out, c_fft, n, n);
            }
        }

        template<class Context>
        inline void update_impl(
            const typename Context::G1Point& a0,
            const MiddleProductPrecomputedBImpl<Context>& precomputed_b,
            std::span<GT> out,
            TimingStats* timings)
        {
            ScopedTimer total_timer(timings, "middle_product_update.total_ns");
            const std::size_t n = precomputed_b.n;
            BLST_ASSERT(n > 0);
            BLST_ASSERT(precomputed_b.b_tail.size() == n);
            BLST_ASSERT(out.size() >= n);

            if (Context::G1Point::is_inf(a0) != 0) {
                std::fill(out.begin(), out.begin() + n, GT::identity());
                return;
            }

            {
                ScopedTimer timer(timings, "middle_product_update.pointwise_ns");
                if (precomputed_b.pairing_mode == MiddleProductPairingMode::PrecomputedLines &&
                    precomputed_b.tail_precomputed_lines.size() == n) {
                    const bool has_zero_mask = (precomputed_b.tail_b_is_zero.size() == n);
                    for (std::size_t i = 0; i < n; ++i) {
                        if (has_zero_mask && precomputed_b.tail_b_is_zero[i] != 0) {
                            out[i] = GT::identity();
                            continue;
                        }
                        out[i] = Context::to_public_gt(
                            Context::pointwise_pair_lines(precomputed_b.tail_precomputed_lines[i], a0));
                    }
                } else {
                    for (std::size_t i = 0; i < n; ++i) {
                        out[i] = Context::to_public_gt(
                            Context::pointwise_pair_direct(a0, precomputed_b.b_tail[i]));
                    }
                }
            }
        }

        template<class Context>
        inline MiddleProductPrecomputedBImpl<Context> precompute_middle_product_b_public(
            std::span<const blst::G2Point> b,
            std::size_t n,
            TimingStats* timings,
            MiddleProductPairingMode pairing_mode)
        {
            const std::vector<typename Context::G2Point> converted = convert_g2_inputs<Context>(b);
            return precompute_middle_product_b_impl<Context>(
                std::span<const typename Context::G2Point>(converted),
                n,
                timings,
                pairing_mode);
        }

        template<class Context>
        inline void middle_product_public(
            std::span<const blst::G1Point> a,
            const MiddleProductPrecomputedBImpl<Context>& precomputed_b,
            std::span<GT> out,
            MiddleProductWorkspaceImpl<Context>& workspace,
            TimingStats* timings)
        {
            const std::vector<typename Context::G1Point> converted_a = convert_g1_inputs<Context>(a);
            middle_product_impl<Context>(
                std::span<const typename Context::G1Point>(converted_a),
                precomputed_b,
                out,
                workspace,
                timings);
        }

        template<class Context>
        inline void update_public(
            const blst::G1Point& a0,
            const MiddleProductPrecomputedBImpl<Context>& precomputed_b,
            std::span<GT> out,
            TimingStats* timings)
        {
            const typename Context::G1Point converted_a0 = Context::to_backend_g1(a0);
            update_impl<Context>(converted_a0, precomputed_b, out, timings);
        }

        template<class Context>
        inline void middle_product_public(
            std::span<const blst::G1Point> a,
            std::span<const blst::G2Point> b,
            std::span<GT> out,
            TimingStats* timings)
        {
            ScopedTimer total_timer(timings, "middle_product.total_ns");
            const std::size_t n = a.size();
            BLST_ASSERT(n > 0);
            BLST_ASSERT(b.size() == 2 * n);
            BLST_ASSERT(out.size() >= n);

            const auto precomputed =
                precompute_middle_product_b_public<Context>(b, n, timings, MiddleProductPairingMode::Direct);
            MiddleProductWorkspaceImpl<Context> workspace;
            middle_product_public<Context>(a, precomputed, out, workspace, timings);
        }

    } // namespace detail

    using MiddleProductScalarContext = detail::ScalarMiddleProductContext;
    using MiddleProductAvxContext = detail::AvxMiddleProductContext;
    using MiddleProductDefaultContext = detail::DefaultMiddleProductContext;

    using MiddleProductPrecomputedBScalar =
        detail::MiddleProductPrecomputedBImpl<MiddleProductScalarContext>;
    using MiddleProductPrecomputedBAvx =
        detail::MiddleProductPrecomputedBImpl<MiddleProductAvxContext>;
    using MiddleProductPrecomputedB =
        detail::MiddleProductPrecomputedBImpl<MiddleProductDefaultContext>;

    using MiddleProductWorkspaceScalar =
        detail::MiddleProductWorkspaceImpl<MiddleProductScalarContext>;
    using MiddleProductWorkspaceAvx =
        detail::MiddleProductWorkspaceImpl<MiddleProductAvxContext>;
    using MiddleProductWorkspace =
        detail::MiddleProductWorkspaceImpl<MiddleProductDefaultContext>;

    inline MiddleProductPrecomputedBScalar precompute_middle_product_b_scalar(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr,
        MiddleProductPairingMode pairing_mode = MiddleProductPairingMode::Direct)
    {
        BLST_ASSERT(n > 0);
        BLST_ASSERT(b.size() == 2 * n);
        return detail::precompute_middle_product_b_public<MiddleProductScalarContext>(
            b,
            n,
            timings,
            pairing_mode);
    }

    inline MiddleProductPrecomputedBAvx precompute_middle_product_b_avx(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr,
        MiddleProductPairingMode pairing_mode = MiddleProductPairingMode::Direct)
    {
        BLST_ASSERT(n > 0);
        BLST_ASSERT(b.size() == 2 * n);
        return detail::precompute_middle_product_b_public<MiddleProductAvxContext>(
            b,
            n,
            timings,
            pairing_mode);
    }

    inline MiddleProductPrecomputedB precompute_middle_product_b(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr,
        MiddleProductPairingMode pairing_mode = MiddleProductPairingMode::Direct)
    {
#if BLST_DEFAULT_AVX
        return precompute_middle_product_b_avx(b, n, timings, pairing_mode);
#else
        return precompute_middle_product_b_scalar(b, n, timings, pairing_mode);
#endif
    }

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
        poly_mult_fft<G1Point, G2Point, GT, Fr, detail::BasicMiddleProductCoeffCtx>(
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

    inline void middle_product_scalar(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedBScalar& precomputed_b,
        std::span<GT> out,
        MiddleProductWorkspaceScalar& workspace,
        TimingStats* timings = nullptr)
    {
        BLST_ASSERT(precomputed_b.n == a.size());
        detail::middle_product_public<MiddleProductScalarContext>(
            a,
            precomputed_b,
            out,
            workspace,
            timings);
    }

    inline void middle_product_avx(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedBAvx& precomputed_b,
        std::span<GT> out,
        MiddleProductWorkspaceAvx& workspace,
        TimingStats* timings = nullptr)
    {
        BLST_ASSERT(precomputed_b.n == a.size());
        detail::middle_product_public<MiddleProductAvxContext>(
            a,
            precomputed_b,
            out,
            workspace,
            timings);
    }

    inline void middle_product(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedB& precomputed_b,
        std::span<GT> out,
        MiddleProductWorkspace& workspace,
        TimingStats* timings = nullptr)
    {
#if BLST_DEFAULT_AVX
        middle_product_avx(a, precomputed_b, out, workspace, timings);
#else
        middle_product_scalar(a, precomputed_b, out, workspace, timings);
#endif
    }

    inline void middle_product_scalar(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        detail::middle_product_public<MiddleProductScalarContext>(a, b, out, timings);
    }

    inline void middle_product_avx(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        detail::middle_product_public<MiddleProductAvxContext>(a, b, out, timings);
    }

    inline void middle_product(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
#if BLST_DEFAULT_AVX
        middle_product_avx(a, b, out, timings);
#else
        middle_product_scalar(a, b, out, timings);
#endif
    }

    inline void middle_product_scalar(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedBScalar& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        MiddleProductWorkspaceScalar workspace;
        middle_product_scalar(a, precomputed_b, out, workspace, timings);
    }

    inline void middle_product_avx(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedBAvx& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        MiddleProductWorkspaceAvx workspace;
        middle_product_avx(a, precomputed_b, out, workspace, timings);
    }

    inline void middle_product(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedB& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        MiddleProductWorkspace workspace;
        middle_product(a, precomputed_b, out, workspace, timings);
    }

    inline void update_scalar(
        const G1Point& a0,
        const MiddleProductPrecomputedBScalar& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        detail::update_public<MiddleProductScalarContext>(a0, precomputed_b, out, timings);
    }

    inline void update_avx(
        const G1Point& a0,
        const MiddleProductPrecomputedBAvx& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        detail::update_public<MiddleProductAvxContext>(a0, precomputed_b, out, timings);
    }

    inline void update(
        const G1Point& a0,
        const MiddleProductPrecomputedB& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
#if BLST_DEFAULT_AVX
        update_avx(a0, precomputed_b, out, timings);
#else
        update_scalar(a0, precomputed_b, out, timings);
#endif
    }

    inline void combine(
        std::span<const GT> mid,
        std::span<const GT> cnst,
        std::span<GT> out)
    {
        BLST_ASSERT(mid.size() == cnst.size());
        BLST_ASSERT(out.size() >= mid.size());

        for (std::size_t i = 0; i < mid.size(); ++i) {
            GT::add(out[i], mid[i], cnst[i]);
        }
    }

} // namespace blst

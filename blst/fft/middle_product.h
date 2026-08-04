#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>

#include "blst/common/defines.h"
#include "blst/common/timer.h"
#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/domain.h"
#include "blst/fft/primitveRoot.h"
#include "blst/fields/Fr.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/point.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/curve/GT.h"
#include "blst/avx/curve/map_to_G.h"
#include "blst/avx/fft/gt_context.h"
#endif

namespace blst {

    enum class MiddleProductAlgorithm {
        Fft,
        Multipair,
        Auto,
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
        return fft_detail::make_roots<Fr, FrCoeffCtx>(n, false);
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
					G::mult_vartime(r, a, reinterpret_cast<const byte*>(&b_non_mont), 256);
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
				Group::mult_vartime(r, a, reinterpret_cast<const byte*>(&b_non_mont), 256);
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

    #if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
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
    #endif

#if BLST_DEFAULT_AVX
        using DefaultMiddleProductContext = AvxMiddleProductContext;
#else
        using DefaultMiddleProductContext = ScalarMiddleProductContext;
#endif

        template<class Context>
        struct MiddleProductPrecomputedBImpl {
            std::size_t n = 0;
            std::size_t fft_size = 0;
            std::vector<typename Context::G2Point> b_tail;
            std::vector<typename Context::G2Point> b_fft;
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
            TimingStats* timings)
        {
            ScopedTimer total_timer(timings, "middle_product_precompute.total_ns");
            BLST_ASSERT(n > 0);
            BLST_ASSERT(b.size() == 2 * n);

            MiddleProductPrecomputedBImpl<Context> precomputed;
            precomputed.n = n;
            precomputed.fft_size = middle_product_fft_size(n);
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
                for (std::size_t i = 0; i < precomputed_b.fft_size; ++i) {
                    c_fft[i] = Context::pointwise_pair_direct(a_fft[i], precomputed_b.b_fft[i]);
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
                for (std::size_t i = 0; i < n; ++i) {
                    out[i] = Context::to_public_gt(
                        Context::pointwise_pair_direct(a0, precomputed_b.b_tail[i]));
                }
            }
        }

        template<class Context>
        inline MiddleProductPrecomputedBImpl<Context> precompute_middle_product_b_public(
            std::span<const blst::G2Point> b,
            std::size_t n,
            TimingStats* timings)
        {
            const std::vector<typename Context::G2Point> converted = convert_g2_inputs<Context>(b);
            return precompute_middle_product_b_impl<Context>(
                std::span<const typename Context::G2Point>(converted),
                n,
                timings);
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
                precompute_middle_product_b_public<Context>(b, n, timings);
            MiddleProductWorkspaceImpl<Context> workspace;
            middle_product_public<Context>(a, precomputed, out, workspace, timings);
        }

    } // namespace detail

    using MiddleProductScalarContext = detail::ScalarMiddleProductContext;
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    using MiddleProductAvxContext = detail::AvxMiddleProductContext;
#endif
    using MiddleProductDefaultContext = detail::DefaultMiddleProductContext;

    using MiddleProductPrecomputedBScalar =
        detail::MiddleProductPrecomputedBImpl<MiddleProductScalarContext>;
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    using MiddleProductPrecomputedBAvx =
        detail::MiddleProductPrecomputedBImpl<MiddleProductAvxContext>;
#endif
    using MiddleProductPrecomputedB =
        detail::MiddleProductPrecomputedBImpl<MiddleProductDefaultContext>;

    using MiddleProductWorkspaceScalar =
        detail::MiddleProductWorkspaceImpl<MiddleProductScalarContext>;
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    using MiddleProductWorkspaceAvx =
        detail::MiddleProductWorkspaceImpl<MiddleProductAvxContext>;
#endif
    using MiddleProductWorkspace =
        detail::MiddleProductWorkspaceImpl<MiddleProductDefaultContext>;

    template<class Context>
    struct MiddleProductPrecomputedBMultipairImpl {
        std::size_t n = 0;
        std::vector<typename Context::PreparedLineTable> prepared_lines;
        std::vector<uint8_t> b_is_zero;
    };

    template<class Context>
    struct MiddleProductWorkspaceMultipairImpl {
        std::vector<typename Context::G1AffinePoint> affine_a_nonzero;
        std::vector<std::size_t> a_indices_nonzero;
    };

    template<>
    struct MiddleProductWorkspaceMultipairImpl<MiddleProductScalarContext> {
        std::vector<MiddleProductScalarContext::G1AffinePoint> affine_a_nonzero;
        std::vector<std::size_t> a_indices_nonzero;
        std::vector<scalar::G1AffinePoint> a_line_points;
        std::vector<const Fp6*> line_batch_scalar;
        std::vector<scalar::G1AffinePoint> px2_batch_scalar;
    };

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    template<>
    struct MiddleProductWorkspaceMultipairImpl<MiddleProductAvxContext> {
        std::vector<MiddleProductAvxContext::G1AffinePoint> affine_a_nonzero;
        std::vector<std::size_t> a_indices_nonzero;
        std::vector<const MiddleProductAvxContext::G1AffinePoint*> p_batch;
        std::vector<const avx::PreparedLine*> line_batch_avx;
        std::vector<avx::GT> miller_batch;
        std::vector<avx::GT> final_batch;
    };
#endif

#if BLST_DEFAULT_AVX
    using MiddleProductDefaultMultipairContext = MiddleProductAvxContext;
#else
    using MiddleProductDefaultMultipairContext = MiddleProductScalarContext;
#endif

    using MiddleProductPrecomputedBMultipairScalar =
        MiddleProductPrecomputedBMultipairImpl<MiddleProductScalarContext>;
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    using MiddleProductPrecomputedBMultipairAvx =
        MiddleProductPrecomputedBMultipairImpl<MiddleProductAvxContext>;
#endif
    using MiddleProductPrecomputedBMultipair =
        MiddleProductPrecomputedBMultipairImpl<MiddleProductDefaultMultipairContext>;

    using MiddleProductWorkspaceMultipairScalar =
        MiddleProductWorkspaceMultipairImpl<MiddleProductScalarContext>;
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    using MiddleProductWorkspaceMultipairAvx =
        MiddleProductWorkspaceMultipairImpl<MiddleProductAvxContext>;
#endif
    using MiddleProductWorkspaceMultipair =
        MiddleProductWorkspaceMultipairImpl<MiddleProductDefaultMultipairContext>;

    inline bool middle_product_is_dense_mask(std::span<const uint8_t> is_zero)
    {
        return std::find(is_zero.begin(), is_zero.end(), static_cast<uint8_t>(1)) == is_zero.end();
    }

    template<class Context>
    inline std::vector<typename Context::G1Point> middle_product_convert_multipair_g1_inputs(
        std::span<const G1Point> input)
    {
        std::vector<typename Context::G1Point> out;
        out.reserve(input.size());
        for (const G1Point& point : input) {
            out.push_back(Context::to_backend_g1(point));
        }
        return out;
    }

    template<class Context>
    inline std::vector<typename Context::G2Point> middle_product_convert_multipair_g2_inputs(
        std::span<const G2Point> input)
    {
        std::vector<typename Context::G2Point> out;
        out.reserve(input.size());
        for (const G2Point& point : input) {
            out.push_back(Context::to_backend_g2(point));
        }
        return out;
    }

    template<class Context>
    inline void middle_product_collect_multipair_g1_affine(
        std::vector<typename Context::G1AffinePoint>& out_affine_nonzero,
        std::vector<std::size_t>& out_indices_nonzero,
        std::span<const typename Context::G1Point> input)
    {
        out_affine_nonzero.clear();
        out_indices_nonzero.clear();
        out_affine_nonzero.reserve(input.size());
        out_indices_nonzero.reserve(input.size());

        for (std::size_t i = 0; i < input.size(); ++i) {
            if (Context::G1Point::is_inf(input[i]) != 0) {
                continue;
            }
            typename Context::G1AffinePoint affine = {};
            Context::G1Point::to_affine(affine, input[i]);
            out_affine_nonzero.push_back(affine);
            out_indices_nonzero.push_back(i);
        }
    }

    inline void middle_product_prepare_scalar_line_eval_points(
        std::vector<scalar::G1AffinePoint>& out_px2,
        std::span<const scalar::G1AffinePoint> input)
    {
        out_px2.clear();
        out_px2.reserve(input.size());
        for (const scalar::G1AffinePoint& point : input) {
            scalar::G1AffinePoint px2 = {};
            scalar::Fp::add(px2.x, point.x, point.x);
            scalar::Fp::neg(px2.x, px2.x);
            scalar::Fp::add(px2.y, point.y, point.y);
            out_px2.push_back(px2);
        }
    }

    template<class Context>
    inline MiddleProductPrecomputedBMultipairImpl<Context> precompute_middle_product_b_multipair_impl(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr)
    {
        ScopedTimer total_timer(timings, "middle_product_multipair_precompute.total_ns");
        BLST_ASSERT(n > 0);
        BLST_ASSERT(b.size() == 2 * n);

        MiddleProductPrecomputedBMultipairImpl<Context> precomputed;
        precomputed.n = n;
        precomputed.prepared_lines.resize(2 * n);
        precomputed.b_is_zero.resize(2 * n);

        const std::vector<typename Context::G2Point> converted =
            middle_product_convert_multipair_g2_inputs<Context>(b);
        for (std::size_t j = 0; j < converted.size(); ++j) {
            const bool is_zero = (Context::G2Point::is_inf(converted[j]) != 0);
            precomputed.b_is_zero[j] = static_cast<uint8_t>(is_zero ? 1 : 0);
            if (is_zero) {
                continue;
            }
            Context::prepare_lines(precomputed.prepared_lines[j], converted[j]);
        }

        return precomputed;
    }

    inline GT middle_product_multipair_scalar_window_dense(
        std::span<const scalar::G1AffinePoint> a_line_points,
        const MiddleProductPrecomputedBMultipairScalar& precomputed_b,
        std::size_t output_index,
        std::vector<const Fp6*>& line_batch)
    {
        const std::size_t n = precomputed_b.n;
        line_batch.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            line_batch[i] = precomputed_b.prepared_lines[n + output_index - i].lines;
        }

        scalar::GT miller = scalar::GT::identity();
        scalar::GT::miller_loop_lines_n_preprocessed(
            miller,
            line_batch.data(),
            a_line_points.data(),
            line_batch.size());

        scalar::GT final_pair = scalar::GT::identity();
        scalar::GT::final_mult(final_pair, miller);
        return MiddleProductScalarContext::to_public_gt(final_pair);
    }

    inline GT middle_product_multipair_scalar_window_sparse(
        std::span<const scalar::G1AffinePoint> a_line_points,
        std::span<const std::size_t> a_indices_nonzero,
        const MiddleProductPrecomputedBMultipairScalar& precomputed_b,
        std::size_t output_index,
        std::vector<const Fp6*>& line_batch,
        std::vector<scalar::G1AffinePoint>& px2_batch)
    {
        line_batch.clear();
        px2_batch.clear();

        for (std::size_t nz = 0; nz < a_line_points.size(); ++nz) {
            const std::size_t i = a_indices_nonzero[nz];
            const std::size_t j = precomputed_b.n + output_index - i;
            if (precomputed_b.b_is_zero[j] != 0) {
                continue;
            }
            px2_batch.push_back(a_line_points[nz]);
            line_batch.push_back(precomputed_b.prepared_lines[j].lines);
        }

        if (line_batch.empty()) {
            return GT::identity();
        }

        scalar::GT miller = scalar::GT::identity();
        scalar::GT::miller_loop_lines_n_preprocessed(
            miller,
            line_batch.data(),
            px2_batch.data(),
            line_batch.size());

        scalar::GT final_pair = scalar::GT::identity();
        scalar::GT::final_mult(final_pair, miller);
        return MiddleProductScalarContext::to_public_gt(final_pair);
    }

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    inline avx::GT middle_product_multipair_avx_window_dense(
        std::span<const avx::G1AffinePoint> a_affine_nonzero,
        const MiddleProductPrecomputedBMultipairAvx& precomputed_b,
        std::size_t output_index,
        std::vector<const avx::PreparedLine*>& line_batch)
    {
        const std::size_t n = precomputed_b.n;
        line_batch.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            line_batch[i] = precomputed_b.prepared_lines[n + output_index - i].lines;
        }

        avx::GT miller = avx::GT::identity();
        avx::GT::miller_loop_lines_n_preprocessed(
            miller,
            line_batch.data(),
            a_affine_nonzero.data(),
            line_batch.size());
        return miller;
    }

    inline avx::GT middle_product_multipair_avx_window_sparse(
        std::span<const avx::G1AffinePoint> a_affine_nonzero,
        std::span<const std::size_t> a_indices_nonzero,
        const MiddleProductPrecomputedBMultipairAvx& precomputed_b,
        std::size_t output_index,
        std::vector<const avx::G1AffinePoint*>& p_batch,
        std::vector<const avx::PreparedLine*>& line_batch)
    {
        p_batch.clear();
        line_batch.clear();

        for (std::size_t nz = 0; nz < a_affine_nonzero.size(); ++nz) {
            const std::size_t i = a_indices_nonzero[nz];
            const std::size_t j = precomputed_b.n + output_index - i;
            if (precomputed_b.b_is_zero[j] != 0) {
                continue;
            }
            p_batch.push_back(&a_affine_nonzero[nz]);
            line_batch.push_back(precomputed_b.prepared_lines[j].lines);
        }

        if (line_batch.empty()) {
            return avx::GT::identity();
        }

        avx::GT miller = avx::GT::identity();
        avx::GT::miller_loop_lines_n(
            miller,
            line_batch.data(),
            p_batch.data(),
            line_batch.size());
        return miller;
    }

    inline void middle_product_finalize_avx_millers(
        std::span<GT> out,
        std::span<const avx::GT> millers,
        std::vector<avx::GT>& final_pairs)
    {
        BLST_ASSERT(out.size() >= millers.size());
        final_pairs.resize(millers.size());
        avx::GT::final_mult_n(final_pairs.data(), millers.data(), millers.size());
        for (std::size_t i = 0; i < millers.size(); ++i) {
            out[i] = MiddleProductAvxContext::to_public_gt(avx::AvxGt52(final_pairs[i].data));
        }
    }
#endif

    template<class Context>
    inline void middle_product_multipair_impl(
        std::span<const typename Context::G1Point> a,
        const MiddleProductPrecomputedBMultipairImpl<Context>& precomputed_b,
        std::span<GT> out,
        MiddleProductWorkspaceMultipairImpl<Context>& workspace,
        TimingStats* timings = nullptr)
    {
        ScopedTimer total_timer(timings, "middle_product_multipair.total_ns");
        const std::size_t n = a.size();
        BLST_ASSERT(n > 0);
        BLST_ASSERT(precomputed_b.n == n);
        BLST_ASSERT(precomputed_b.prepared_lines.size() == 2 * n);
        BLST_ASSERT(precomputed_b.b_is_zero.size() == 2 * n);
        BLST_ASSERT(out.size() >= n);

        middle_product_collect_multipair_g1_affine<Context>(
            workspace.affine_a_nonzero,
            workspace.a_indices_nonzero,
            a);

        {
            ScopedTimer timer(timings, "middle_product_multipair.pointwise_ns");
            if constexpr (std::is_same_v<Context, MiddleProductScalarContext>) {
                middle_product_prepare_scalar_line_eval_points(
                    workspace.a_line_points,
                    std::span<const scalar::G1AffinePoint>(workspace.affine_a_nonzero));

                const bool dense_a = (workspace.a_indices_nonzero.size() == n);
                const bool dense_b = middle_product_is_dense_mask(std::span<const uint8_t>(precomputed_b.b_is_zero));

                if (dense_a && dense_b) {
                    for (std::size_t k = 0; k < n; ++k) {
                        out[k] = middle_product_multipair_scalar_window_dense(
                            std::span<const scalar::G1AffinePoint>(workspace.a_line_points),
                            precomputed_b,
                            k,
                            workspace.line_batch_scalar);
                    }
                } else {
                    for (std::size_t k = 0; k < n; ++k) {
                        out[k] = middle_product_multipair_scalar_window_sparse(
                            std::span<const scalar::G1AffinePoint>(workspace.a_line_points),
                            std::span<const std::size_t>(workspace.a_indices_nonzero),
                            precomputed_b,
                            k,
                            workspace.line_batch_scalar,
                            workspace.px2_batch_scalar);
                    }
                }
            } else {
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
                const bool dense_a = (workspace.a_indices_nonzero.size() == n);
                const bool dense_b = middle_product_is_dense_mask(std::span<const uint8_t>(precomputed_b.b_is_zero));
                workspace.miller_batch.resize(n);

                if (dense_a && dense_b) {
                    for (std::size_t k = 0; k < n; ++k) {
                        workspace.miller_batch[k] = middle_product_multipair_avx_window_dense(
                            std::span<const avx::G1AffinePoint>(workspace.affine_a_nonzero),
                            precomputed_b,
                            k,
                            workspace.line_batch_avx);
                    }
                } else {
                    workspace.p_batch.reserve(workspace.affine_a_nonzero.size());
                    workspace.line_batch_avx.reserve(workspace.affine_a_nonzero.size());
                    for (std::size_t k = 0; k < n; ++k) {
                        workspace.miller_batch[k] = middle_product_multipair_avx_window_sparse(
                            std::span<const avx::G1AffinePoint>(workspace.affine_a_nonzero),
                            std::span<const std::size_t>(workspace.a_indices_nonzero),
                            precomputed_b,
                            k,
                            workspace.p_batch,
                            workspace.line_batch_avx);
                    }
                }
                middle_product_finalize_avx_millers(
                    out.first(n),
                    std::span<const avx::GT>(workspace.miller_batch),
                    workspace.final_batch);
#endif
            }
        }
    }

    template<class Context>
    inline void middle_product_multipair_public(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedBMultipairImpl<Context>& precomputed_b,
        std::span<GT> out,
        MiddleProductWorkspaceMultipairImpl<Context>& workspace,
        TimingStats* timings = nullptr)
    {
        const std::vector<typename Context::G1Point> converted_a =
            middle_product_convert_multipair_g1_inputs<Context>(a);
        middle_product_multipair_impl<Context>(
            std::span<const typename Context::G1Point>(converted_a),
            precomputed_b,
            out,
            workspace,
            timings);
    }

    template<class Context>
    inline void update_multipair_impl(
        const typename Context::G1Point& a0,
        const MiddleProductPrecomputedBMultipairImpl<Context>& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        ScopedTimer total_timer(timings, "middle_product_multipair_update.total_ns");
        const std::size_t n = precomputed_b.n;
        BLST_ASSERT(n > 0);
        BLST_ASSERT(precomputed_b.prepared_lines.size() == 2 * n);
        BLST_ASSERT(precomputed_b.b_is_zero.size() == 2 * n);
        BLST_ASSERT(out.size() >= n);

        if (Context::G1Point::is_inf(a0) != 0) {
            std::fill(out.begin(), out.begin() + n, GT::identity());
            return;
        }

        typename Context::G1AffinePoint affine_a0 = {};
        Context::G1Point::to_affine(affine_a0, a0);

        {
            ScopedTimer timer(timings, "middle_product_multipair_update.pointwise_ns");
            if constexpr (std::is_same_v<Context, MiddleProductScalarContext>) {
                for (std::size_t k = 0; k < n; ++k) {
                    if (precomputed_b.b_is_zero[n + k] != 0) {
                        out[k] = GT::identity();
                        continue;
                    }
                    scalar::GT miller = scalar::GT::identity();
                    scalar::GT::miller_loop_lines(
                        miller,
                        precomputed_b.prepared_lines[n + k].lines,
                        affine_a0);
                    scalar::GT final_pair = scalar::GT::identity();
                    scalar::GT::final_mult(final_pair, miller);
                    out[k] = MiddleProductScalarContext::to_public_gt(final_pair);
                }
            } else {
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
                std::vector<avx::GT> millers(n, avx::GT::identity());
                std::vector<avx::GT> finals;
                for (std::size_t k = 0; k < n; ++k) {
                    if (precomputed_b.b_is_zero[n + k] != 0) {
                        continue;
                    }
                    avx::GT::miller_loop_lines(
                        millers[k],
                        precomputed_b.prepared_lines[n + k].lines,
                        affine_a0);
                }
                middle_product_finalize_avx_millers(
                    out.first(n),
                    std::span<const avx::GT>(millers),
                    finals);
#endif
            }
        }
    }

    inline MiddleProductPrecomputedBScalar precompute_middle_product_b_scalar(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr)
    {
        BLST_ASSERT(n > 0);
        BLST_ASSERT(b.size() == 2 * n);
        return detail::precompute_middle_product_b_public<MiddleProductScalarContext>(
            b,
            n,
            timings);
    }

    #if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    inline MiddleProductPrecomputedBAvx precompute_middle_product_b_avx(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr)
    {
        BLST_ASSERT(n > 0);
        BLST_ASSERT(b.size() == 2 * n);
        return detail::precompute_middle_product_b_public<MiddleProductAvxContext>(
            b,
            n,
            timings);
    }
    #endif

    inline MiddleProductPrecomputedB precompute_middle_product_b(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr)
    {
#if BLST_DEFAULT_AVX
        return precompute_middle_product_b_avx(b, n, timings);
#else
        return precompute_middle_product_b_scalar(b, n, timings);
#endif
    }

    inline MiddleProductPrecomputedBMultipairScalar precompute_middle_product_b_multipair_scalar(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr)
    {
        BLST_ASSERT(n > 0);
        BLST_ASSERT(b.size() == 2 * n);
        return precompute_middle_product_b_multipair_impl<MiddleProductScalarContext>(b, n, timings);
    }

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    inline MiddleProductPrecomputedBMultipairAvx precompute_middle_product_b_multipair_avx(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr)
    {
        BLST_ASSERT(n > 0);
        BLST_ASSERT(b.size() == 2 * n);
        return precompute_middle_product_b_multipair_impl<MiddleProductAvxContext>(b, n, timings);
    }
#endif

    inline MiddleProductPrecomputedBMultipair precompute_middle_product_b_multipair(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr)
    {
#if BLST_DEFAULT_AVX
        return precompute_middle_product_b_multipair_avx(b, n, timings);
#else
        return precompute_middle_product_b_multipair_scalar(b, n, timings);
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

    inline void middle_product_multipair_scalar(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedBMultipairScalar& precomputed_b,
        std::span<GT> out,
        MiddleProductWorkspaceMultipairScalar& workspace,
        TimingStats* timings = nullptr)
    {
        BLST_ASSERT(precomputed_b.n == a.size());
        middle_product_multipair_public<MiddleProductScalarContext>(
            a,
            precomputed_b,
            out,
            workspace,
            timings);
    }

    #if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
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

    inline void middle_product_multipair_avx(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedBMultipairAvx& precomputed_b,
        std::span<GT> out,
        MiddleProductWorkspaceMultipairAvx& workspace,
        TimingStats* timings = nullptr)
    {
        BLST_ASSERT(precomputed_b.n == a.size());
        middle_product_multipair_public<MiddleProductAvxContext>(
            a,
            precomputed_b,
            out,
            workspace,
            timings);
    }
    #endif

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

    inline void middle_product_multipair(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedBMultipair& precomputed_b,
        std::span<GT> out,
        MiddleProductWorkspaceMultipair& workspace,
        TimingStats* timings = nullptr)
    {
#if BLST_DEFAULT_AVX
        middle_product_multipair_avx(a, precomputed_b, out, workspace, timings);
#else
        middle_product_multipair_scalar(a, precomputed_b, out, workspace, timings);
#endif
    }

    inline constexpr std::size_t middle_product_auto_multipair_max_n = 24;

    inline void middle_product_scalar(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        detail::middle_product_public<MiddleProductScalarContext>(a, b, out, timings);
    }

    #if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    inline void middle_product_avx(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        detail::middle_product_public<MiddleProductAvxContext>(a, b, out, timings);
    }
    #endif

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

    inline void middle_product_multipair_scalar(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        const auto precomputed = precompute_middle_product_b_multipair_scalar(b, a.size(), timings);
        MiddleProductWorkspaceMultipairScalar workspace;
        middle_product_multipair_scalar(a, precomputed, out, workspace, timings);
    }

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    inline void middle_product_multipair_avx(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        const auto precomputed = precompute_middle_product_b_multipair_avx(b, a.size(), timings);
        MiddleProductWorkspaceMultipairAvx workspace;
        middle_product_multipair_avx(a, precomputed, out, workspace, timings);
    }
#endif

    inline void middle_product_multipair(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
#if BLST_DEFAULT_AVX
        middle_product_multipair_avx(a, b, out, timings);
#else
        middle_product_multipair_scalar(a, b, out, timings);
#endif
    }

    inline void middle_product_scalar(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        MiddleProductAlgorithm algorithm,
        TimingStats* timings)
    {
        switch (algorithm) {
        case MiddleProductAlgorithm::Fft:
            middle_product_scalar(a, b, out, timings);
            break;
        case MiddleProductAlgorithm::Multipair:
            middle_product_multipair_scalar(a, b, out, timings);
            break;
        case MiddleProductAlgorithm::Auto:
            if (a.size() <= middle_product_auto_multipair_max_n) {
                middle_product_multipair_scalar(a, b, out, timings);
            } else {
                middle_product_scalar(a, b, out, timings);
            }
            break;
        default:
            BLST_ASSERT(false && "Unsupported middle product algorithm");
        }
    }

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    inline void middle_product_avx(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        MiddleProductAlgorithm algorithm,
        TimingStats* timings)
    {
        switch (algorithm) {
        case MiddleProductAlgorithm::Fft:
            middle_product_avx(a, b, out, timings);
            break;
        case MiddleProductAlgorithm::Multipair:
            middle_product_multipair_avx(a, b, out, timings);
            break;
        case MiddleProductAlgorithm::Auto:
            if (a.size() <= middle_product_auto_multipair_max_n) {
                middle_product_multipair_avx(a, b, out, timings);
            } else {
                middle_product_avx(a, b, out, timings);
            }
            break;
        default:
            BLST_ASSERT(false && "Unsupported middle product algorithm");
        }
    }
#endif

    inline void middle_product(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        MiddleProductAlgorithm algorithm,
        TimingStats* timings)
    {
        switch (algorithm) {
        case MiddleProductAlgorithm::Fft:
            middle_product(a, b, out, timings);
            break;
        case MiddleProductAlgorithm::Multipair:
            middle_product_multipair(a, b, out, timings);
            break;
        case MiddleProductAlgorithm::Auto:
            if (a.size() <= middle_product_auto_multipair_max_n) {
                middle_product_multipair(a, b, out, timings);
            } else {
                middle_product(a, b, out, timings);
            }
            break;
        default:
            BLST_ASSERT(false && "Unsupported middle product algorithm");
        }
    }

    inline void middle_product_scalar(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        MiddleProductAlgorithm algorithm)
    {
        middle_product_scalar(a, b, out, algorithm, nullptr);
    }

    #if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    inline void middle_product_avx(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        MiddleProductAlgorithm algorithm)
    {
        middle_product_avx(a, b, out, algorithm, nullptr);
    }
    #endif

    inline void middle_product(
        std::span<const G1Point> a,
        std::span<const G2Point> b,
        std::span<GT> out,
        MiddleProductAlgorithm algorithm)
    {
        middle_product(a, b, out, algorithm, nullptr);
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

    inline void middle_product_multipair_scalar(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedBMultipairScalar& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        MiddleProductWorkspaceMultipairScalar workspace;
        middle_product_multipair_scalar(a, precomputed_b, out, workspace, timings);
    }

    #if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    inline void middle_product_avx(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedBAvx& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        MiddleProductWorkspaceAvx workspace;
        middle_product_avx(a, precomputed_b, out, workspace, timings);
    }

    inline void middle_product_multipair_avx(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedBMultipairAvx& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        MiddleProductWorkspaceMultipairAvx workspace;
        middle_product_multipair_avx(a, precomputed_b, out, workspace, timings);
    }
    #endif

    inline void middle_product(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedB& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        MiddleProductWorkspace workspace;
        middle_product(a, precomputed_b, out, workspace, timings);
    }

    inline void middle_product_multipair(
        std::span<const G1Point> a,
        const MiddleProductPrecomputedBMultipair& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        MiddleProductWorkspaceMultipair workspace;
        middle_product_multipair(a, precomputed_b, out, workspace, timings);
    }

    inline void update_scalar(
        const G1Point& a0,
        const MiddleProductPrecomputedBScalar& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        detail::update_public<MiddleProductScalarContext>(a0, precomputed_b, out, timings);
    }

    inline void update_multipair_scalar(
        const G1Point& a0,
        const MiddleProductPrecomputedBMultipairScalar& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        update_multipair_impl<MiddleProductScalarContext>(
            MiddleProductScalarContext::to_backend_g1(a0),
            precomputed_b,
            out,
            timings);
    }

    #if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    inline void update_avx(
        const G1Point& a0,
        const MiddleProductPrecomputedBAvx& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        detail::update_public<MiddleProductAvxContext>(a0, precomputed_b, out, timings);
    }

    inline void update_multipair_avx(
        const G1Point& a0,
        const MiddleProductPrecomputedBMultipairAvx& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
        update_multipair_impl<MiddleProductAvxContext>(
            MiddleProductAvxContext::to_backend_g1(a0),
            precomputed_b,
            out,
            timings);
    }
    #endif

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

    inline void update_multipair(
        const G1Point& a0,
        const MiddleProductPrecomputedBMultipair& precomputed_b,
        std::span<GT> out,
        TimingStats* timings = nullptr)
    {
#if BLST_DEFAULT_AVX
        update_multipair_avx(a0, precomputed_b, out, timings);
#else
        update_multipair_scalar(a0, precomputed_b, out, timings);
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

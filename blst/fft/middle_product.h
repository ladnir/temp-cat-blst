// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// CatBLST middle-product implementation, locally composed from the BLST
// pairing and FFT primitives with scalar and MatrixNoK AVX backends.

#pragma once

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "blst/common/defines.h"
#include "blst/common/timer.h"
#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/domain.h"
#include "blst/fft/primitive_root.h"
#include "blst/fields/Fr.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/point.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/curve/GT.h"
#include "blst/avx/curve/map_to_G.h"
#include "blst/avx/fft/gt_context.h"
#endif

namespace blst {

    /**
     * @brief Pairing middle-product API and shared scalar/AVX implementation.
     *
     * For `a[0..n)` in G1 and `b[0..2*n)` in G2, evaluation returns GT values
     *
     * `out[k] = product_{i=0}^{n-1} e(a[i], b[n+k-i])`, for `0 <= k < n`.
     *
     * G1 and G2 mean the prime-order subgroups: callers must provide reduced,
     * on-curve subgroup points. Infinity is accepted as the identity. These
     * operations do not perform Release membership validation. Outputs are GT
     * elements; `combine` likewise requires GT inputs rather than arbitrary
     * Fp12 values.
     *
     * The algorithms are variable-time with respect to their public points and
     * sizes. Precomputation and convenience overloads allocate and can throw
     * `std::bad_alloc`; malformed shapes throw `std::invalid_argument`, and an
     * unsupported length or insufficient workspace capacity throws
     * `std::length_error`. FFT workspaces allocate fixed, 64-byte-aligned
     * storage at construction and never resize; selected-native evaluation
     * performs no allocation. Backend operations may still allocate where
     * documented by their concrete context.
     * The AVX multipair workspace also retains all six batch-final-exponentiation
     * buffers, so that stage performs no allocation after its capacity grows.
     * Selected-backend entry points consume their input spans directly; only
     * an explicit call to the non-selected scalar/AVX backend allocates a point
     * conversion vector. Validation and allocation failures occur before output
     * publication begins.
     *
     * A precomputed object is immutable after construction and may be shared by
     * concurrent evaluations when each evaluation owns a distinct workspace and
     * output. One workspace must not be used concurrently. Unless a function
     * explicitly says otherwise, output must not overlap its inputs. `combine`
     * supports exact output aliasing with either input but forbids partial
     * overlap. Precompute fields remain visible for expert inspection, but
     * callers must not modify them after construction.
     */

    // The shared FFT implementation uses signed-int indices. The largest
    // supported transform is therefore the largest power of two below
    // INT_MAX. A middle product of length n uses a 2*n transform.
    inline constexpr std::size_t middle_product_max_fft_size =
        std::size_t{1} << (std::numeric_limits<int>::digits - 1);
    inline constexpr std::size_t middle_product_max_n =
        middle_product_max_fft_size / 2;

    [[nodiscard]] inline constexpr bool middle_product_n_is_supported(
        std::size_t n) noexcept
    {
        return n != 0 && n <= middle_product_max_n;
    }

    [[nodiscard]] inline std::size_t middle_product_b_size(std::size_t n)
    {
        if (n == 0) {
            throw std::invalid_argument("middle product length must be nonzero");
        }
        if (n > middle_product_max_n) {
            throw std::length_error("middle product length exceeds the supported FFT index range");
        }
        return 2 * n;
    }

    inline void validate_middle_product_b_shape(
        std::size_t b_size,
        std::size_t n)
    {
        if (b_size != middle_product_b_size(n)) {
            throw std::invalid_argument("middle product requires exactly 2*n G2 inputs");
        }
    }

    inline void validate_middle_product_output_shape(
        std::size_t n,
        std::size_t out_size)
    {
        if (!middle_product_n_is_supported(n)) {
            (void)middle_product_b_size(n);
        }
        if (out_size < n) {
            throw std::invalid_argument("middle product output span is too small");
        }
    }

    /// Explicit implementation selection. `Auto` uses a benchmark-selected
    /// size threshold whose scalar and AVX recalibration is tracked separately.
    /// Values outside this enumeration violate the native contract and trigger
    /// a Debug assertion; Release performs no work for such a value.
    enum class MiddleProductAlgorithm {
        Fft,
        Multipair,
        Auto,
    };

    /// Stateless phases of one planned middle-product evaluation. Input and
    /// output use typed range functions; every other phase uses the common
    /// phase/stage range executor. Stage barriers and phase ordering belong to
    /// the caller and are deliberately not retained by CatBLST.
    enum class MiddleProductRangePhase : std::uint8_t {
        InputLoad,
        ForwardBitReverse,
        ForwardJ0,
        ForwardNonzero,
        AffineConversion,
        PointwisePairing,
        InverseBitReverse,
        InverseButterflies,
        Normalization,
        Extraction,
    };

    struct alignas(64) MiddleProductPreparedScalarLineTable {
        Fp6 lines[68];
    };

    inline std::size_t middle_product_basic_fft_size(std::size_t n)
    {
        if (n == 0) {
            throw std::invalid_argument("middle product length must be nonzero");
        }

        // The basic oracle transforms a polynomial of length 3*n-1.
        constexpr std::size_t max_n = (middle_product_max_fft_size + 1) / 3;
        if (n > max_n) {
            throw std::length_error("basic middle product length exceeds the supported FFT index range");
        }

        std::size_t fft_size = 1;
        const std::size_t target = (3 * n) - 1;
        while (fft_size < target) {
            fft_size <<= 1;
        }
        return fft_size;
    }

    inline std::size_t middle_product_fft_size(std::size_t n)
    {
        const std::size_t target = middle_product_b_size(n);

        std::size_t fft_size = 1;
        while (fft_size < target) {
            fft_size <<= 1;
        }
        return fft_size;
    }

    inline std::size_t fft_admissible_size(std::size_t x, std::size_t fft_size)
    {
        if (fft_size == 0 || (fft_size & (fft_size - 1)) != 0 ||
            fft_size > middle_product_max_fft_size) {
            throw std::invalid_argument("invalid middle-product FFT size");
        }

        constexpr std::size_t kRoundPow2 = 16;
        if (x == 0) {
            x = 1;
        }
        if (x >= fft_size) {
            return fft_size;
        }

        x = (((x - 1) / kRoundPow2) + 1) * kRoundPow2;

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

    /**
     * Build inverse powers from a complete ordered table
     * `[1, omega, ..., omega^(n-1)]`, where `omega` has order `n`.
     */
    inline std::vector<Fr> make_ifft_roots(std::span<const Fr> roots)
    {
        std::vector<Fr> inv_roots(roots.size());
        if (roots.empty()) {
            return inv_roots;
        }

        BLST_ASSERT((roots.size() & (roots.size() - 1)) == 0);
        BLST_ASSERT(Fr::equal(roots[0], Fr::one()));
        inv_roots[0] = roots[0];
        for (std::size_t i = 1; i < roots.size(); ++i) {
            inv_roots[i] = roots[roots.size() - i];
        }
        return inv_roots;
    }

    namespace detail {

        inline constexpr std::size_t middle_product_workspace_alignment = 64;

        template<class T, std::size_t Alignment>
        struct MiddleProductAlignedAllocator {
            static_assert(Alignment >= alignof(T));
            static_assert((Alignment & (Alignment - 1)) == 0);

            using value_type = T;
            using is_always_equal = std::true_type;

            constexpr MiddleProductAlignedAllocator() noexcept = default;

            template<class U>
            constexpr MiddleProductAlignedAllocator(
                const MiddleProductAlignedAllocator<U, Alignment>&) noexcept
            {
            }

            [[nodiscard]] T* allocate(const std::size_t count)
            {
                if (count > std::numeric_limits<std::size_t>::max() / sizeof(T)) {
                    throw std::bad_array_new_length();
                }
                return static_cast<T*>(::operator new(
                    count * sizeof(T),
                    std::align_val_t{Alignment}));
            }

            void deallocate(T* const allocation, std::size_t) noexcept
            {
                ::operator delete(allocation, std::align_val_t{Alignment});
            }

            template<class U>
            struct rebind {
                using other = MiddleProductAlignedAllocator<U, Alignment>;
            };
        };

        template<class T, class U, std::size_t Alignment>
        constexpr bool operator==(
            const MiddleProductAlignedAllocator<T, Alignment>&,
            const MiddleProductAlignedAllocator<U, Alignment>&) noexcept
        {
            return true;
        }

        template<class T>
        using MiddleProductAlignedVector = std::vector<
            T,
            MiddleProductAlignedAllocator<T, middle_product_workspace_alignment>>;

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
                    G::mult_vartime(
                        r,
                        a,
                        reinterpret_cast<const byte*>(&b_non_mont),
                        256);
                    return r;
                }
            }

            static auto mult(const G1Point& a, const G2Point& b)
            {
                G1AffinePoint p_aff = {};
                G1Point::to_affine(p_aff, a);
                G2AffinePoint q_aff = {};
                G2Point::to_affine(q_aff, b);

                MillerResult miller = {};
                GT::miller_loop(miller, &q_aff, &p_aff);
                GT pair = {};
                GT::final_exp(pair, miller);
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
                Group::mult_vartime(
                    r,
                    a,
                    reinterpret_cast<const byte*>(&b_non_mont),
                    256);
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
            static constexpr bool batch_affine_pointwise = false;

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
                return avx::AvxGtElement::to_default(avx::AvxGtElement::from_scalar(value));
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

            static Fr inverse_fft_size(
                const PreparedInverseRoots&,
                const std::size_t size)
            {
                Fr inverse = {};
                Fr::from_uint64(inverse, static_cast<std::uint64_t>(size));
                Fr::reciprocal(inverse, inverse);
                return inverse;
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

                scalar::MillerResult miller = {};
                scalar::GT::miller_loop(miller, &q_aff, &p_aff);
                GtFftElement pair = {};
                scalar::GT::final_exp(pair, miller);
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
                std::size_t length,
                const FftInverseTruncatedRangePlan& plan)
            {
                fft_inverse_truncated_range<GtFftElement, Fr, GtCoeffCtx>(
                    values,
                    std::span<const Fr>(inv_roots),
                    start,
                    length,
                    plan);
            }

            static void inverse_fft_stage_range(
                std::span<GtFftElement> values,
                const PreparedInverseRoots& inv_roots,
                const FftInverseTruncatedRangePlan& plan,
                std::size_t stage,
                std::size_t begin,
                std::size_t end)
            {
                fft_inverse_stage_range<GtFftElement, Fr, GtCoeffCtx>(
                    values,
                    std::span<const Fr>(inv_roots),
                    plan,
                    stage,
                    begin,
                    end);
            }

            static void normalize_inverse_range(
                std::span<GtFftElement> values,
                const PreparedInverseRoots&,
                const Fr& inverse_size,
                std::size_t begin,
                std::size_t end)
            {
                for (std::size_t idx = begin; idx < end; ++idx) {
                    values[idx] = GtCoeffCtx::mult(values[idx], inverse_size);
                }
            }

            template<class OutputGt>
            static void extract_middle_window(
                std::span<OutputGt> out,
                std::span<const GtFftElement> values,
                std::size_t start,
                std::size_t length)
            {
                static_assert(std::is_same_v<OutputGt, GT>);
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
            using GtFftElement = avx::AvxGtElement;
            using PreparedInverseRoots = avx::AvxGtPreparedRoots;
            using PreparedLineTable = avx::PreparedLineTable;
            // Production AVX policy: normalize each transformed batch once and
            // use prepared G2 lines for every pointwise pairing. This is a
            // compile-time backend property, not a runtime tuning switch.
            static constexpr bool batch_affine_pointwise = true;

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
                return avx::AvxGtElement::to_default(value);
            }

            static PreparedInverseRoots make_inverse_roots(std::span<const Fr> roots)
            {
                const std::vector<Fr> inv_roots = make_ifft_roots(roots);
                PreparedInverseRoots prepared{std::span<const Fr>(inv_roots)};
                prepared.prepare_inverse_scale();
                return prepared;
            }

            static std::size_t inverse_roots_size(const PreparedInverseRoots& roots)
            {
                return roots.size();
            }

            static Fr inverse_fft_size(
                const PreparedInverseRoots& roots,
                const std::size_t size)
            {
                (void)size;
                BLST_ASSERT(roots.has_scale);
                BLST_ASSERT(roots.size() == size);
                return roots.scale_scalar;
            }

            static void prepare_lines(PreparedLineTable& out, const G2Point& point)
            {
                G2AffinePoint point_affine = {};
                G2Point::to_affine(point_affine, point);
                avx::GT::precompute_lines(out.lines, point_affine);
            }

            static void prepare_lines_affine(
                PreparedLineTable& out,
                const G2AffinePoint& point)
            {
                avx::GT::precompute_lines(out.lines, point);
            }

            static GtFftElement pointwise_pair_direct(const G1Point& left, const G2Point& right)
            {
                G1AffinePoint p_aff = {};
                G1Point::to_affine(p_aff, left);
                G2AffinePoint q_aff = {};
                G2Point::to_affine(q_aff, right);

                avx::MillerResult miller = {};
                avx::GT::miller_loop(miller, &q_aff, &p_aff);
                avx::GT pair = {};
                avx::GT::final_exp(pair, miller);
                return GtFftElement(pair.data);
            }

            static GtFftElement pointwise_pair_prepared(
                const G1AffinePoint& left,
                const PreparedLineTable& right)
            {
                avx::MillerResult miller = {};
                avx::GT::miller_loop_lines(miller, right.lines, left);
                avx::GT pair = {};
                avx::GT::final_exp(pair, miller);
                return GtFftElement(pair.data);
            }

            static void batch_to_affine_g1(
                std::span<G1AffinePoint> out,
                std::span<const G1Point> in,
                std::span<typename G1Point::Field_type> prefix_products)
            {
                G1Point::to_affine_batch(out, in, prefix_products);
            }

            static void batch_to_affine_g2(
                std::span<G2AffinePoint> out,
                std::span<const G2Point> in,
                std::span<typename G2Point::Field_type> prefix_products)
            {
                G2Point::to_affine_batch(out, in, prefix_products);
            }

            static GtFftElement identity()
            {
                return GtFftElement::identity();
            }

            static void inverse_fft_middle_window(
                std::span<GtFftElement> values,
                const PreparedInverseRoots& inv_roots,
                std::size_t start,
                std::size_t length,
                const FftInverseTruncatedRangePlan& plan)
            {
                avx::AvxGtCoeffCtx::ScopedRoots scoped_roots(inv_roots);
                fft_inverse_truncated_range<GtFftElement, Fr, avx::AvxGtCoeffCtx>(
                    values,
                    inv_roots.scalar_span(),
                    start,
                    length,
                    plan);
            }

            static void inverse_fft_stage_range(
                std::span<GtFftElement> values,
                const PreparedInverseRoots& inv_roots,
                const FftInverseTruncatedRangePlan& plan,
                std::size_t stage,
                std::size_t begin,
                std::size_t end)
            {
                avx::AvxGtCoeffCtx::ScopedRoots scoped_roots(inv_roots);
                fft_inverse_stage_range<
                    GtFftElement,
                    Fr,
                    avx::AvxGtCoeffCtx>(
                        values,
                        inv_roots.scalar_span(),
                        plan,
                        stage,
                        begin,
                        end);
            }

            static void normalize_inverse_range(
                std::span<GtFftElement> values,
                const PreparedInverseRoots& inv_roots,
                const Fr& inverse_size,
                std::size_t begin,
                std::size_t end)
            {
                avx::AvxGtCoeffCtx::ScopedRoots scoped_roots(inv_roots);
                for (std::size_t idx = begin; idx < end; ++idx) {
                    values[idx] = avx::AvxGtCoeffCtx::mult(
                        values[idx], inverse_size);
                }
            }

            template<class OutputGt>
            static void extract_middle_window(
                std::span<OutputGt> out,
                std::span<const GtFftElement> values,
                std::size_t start,
                std::size_t length)
            {
                for (std::size_t i = 0; i < length; ++i) {
                    if constexpr (std::is_same_v<OutputGt, avx::GT>) {
                        out[i].data = values[start + i].data;
                    } else {
                        static_assert(std::is_same_v<OutputGt, GT>);
                        out[i] = to_public_gt(values[start + i]);
                    }
                }
            }
        };
    #endif

#if BLST_DEFAULT_AVX
        using DefaultMiddleProductContext = AvxMiddleProductContext;
#else
        using DefaultMiddleProductContext = ScalarMiddleProductContext;
#endif

        /// Immutable FFT precomputation for one length-`2*n` G2 operand.
        /// Fields are expert-visible for inspection and testing; modifying any
        /// field after construction invalidates the object's native contract.
        template<class Context>
        struct MiddleProductPrecomputedBImpl {
            std::size_t n = 0;
            std::size_t fft_size = 0;
            std::vector<typename Context::G2Point> b_tail;
            std::vector<typename Context::G2Point> b_fft;
            std::vector<typename Context::PreparedLineTable> prepared_lines;
            std::vector<uint8_t> b_fft_is_zero;
            std::vector<Fr> roots;
            typename Context::PreparedInverseRoots inv_roots;
            Fr inverse_fft_size = {};
            FftTruncatedForwardPlan forward_plan;
            FftInverseTruncatedRangePlan inverse_plan;
        };

        /// Mutable, non-thread-safe FFT evaluation storage reusable for every
        /// precompute no larger than its construction capacity. Buffer sizes
        /// and addresses remain fixed for the lifetime of the workspace.
        template<class Context>
        struct MiddleProductWorkspaceImpl {
        private:
            std::size_t max_fft_capacity_;

        public:
            using G1Point = typename Context::G1Point;
            using G1AffinePoint = typename Context::G1AffinePoint;
            using AffinePrefix = typename Context::G1Point::Field_type;
            using GtFftElement = typename Context::GtFftElement;

            explicit MiddleProductWorkspaceImpl(
                const std::size_t max_fft_capacity)
                : max_fft_capacity_(validate_capacity(max_fft_capacity)),
                  a_fft(max_fft_capacity_),
                  a_fft_affine(
                      Context::batch_affine_pointwise ? max_fft_capacity_ : 0),
                  affine_prefix_products(
                      Context::batch_affine_pointwise ? max_fft_capacity_ : 0),
                  c_fft(max_fft_capacity_)
            {
            }

            [[nodiscard]] std::size_t capacity() const noexcept
            {
                return max_fft_capacity_;
            }

            [[nodiscard]] bool has_capacity(
                const std::size_t required_fft_size) const noexcept
            {
                return required_fft_size <= max_fft_capacity_;
            }

            void require_capacity(const std::size_t required_fft_size) const
            {
                if (!has_capacity(required_fft_size)) {
                    throw std::length_error(
                        "middle-product workspace capacity is too small");
                }
            }

            MiddleProductAlignedVector<G1Point> a_fft;
            MiddleProductAlignedVector<G1AffinePoint> a_fft_affine;
            MiddleProductAlignedVector<AffinePrefix> affine_prefix_products;
            MiddleProductAlignedVector<GtFftElement> c_fft;

        private:
            [[nodiscard]] static std::size_t validate_capacity(
                const std::size_t max_fft_capacity)
            {
                if (!fft_detail::valid_transform_size(max_fft_capacity) ||
                    max_fft_capacity > middle_product_max_fft_size) {
                    throw std::invalid_argument(
                        "invalid middle-product workspace capacity");
                }
                return max_fft_capacity;
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
            validate_middle_product_b_shape(b.size(), n);

            MiddleProductPrecomputedBImpl<Context> precomputed;
            precomputed.n = n;
            precomputed.fft_size = middle_product_fft_size(n);
            precomputed.b_tail.assign(n, Context::g2_infinity());
            std::vector<typename Context::G2Point> b_fft(
                precomputed.fft_size, Context::g2_infinity());

            {
                ScopedTimer timer(timings, "middle_product_precompute.roots_ns");
                precomputed.roots = make_fft_roots(precomputed.fft_size);
            }
            precomputed.inv_roots = Context::make_inverse_roots(precomputed.roots);
            precomputed.inverse_fft_size = Context::inverse_fft_size(
                precomputed.inv_roots, precomputed.fft_size);

            std::vector<typename Context::G2Point> b1_fft(precomputed.fft_size, Context::g2_infinity());
            for (std::size_t i = 0; i < n; ++i) {
                b_fft[i] = b[i];
                precomputed.b_tail[i] = b[n + i];
                b1_fft[i] = precomputed.b_tail[i];
            }

            const std::size_t xn_b = fft_admissible_size(n, precomputed.fft_size);
            {
                ScopedTimer timer(timings, "middle_product_precompute.plans_ns");
                precomputed.forward_plan =
                    make_fft_truncated_forward_plan(precomputed.fft_size, xn_b);
                precomputed.inverse_plan =
                    make_fft_inverse_truncated_range_plan(precomputed.fft_size, n, n);
            }
            {
                ScopedTimer timer(timings, "middle_product_precompute.fft_b0_ns");
                fft_forward_truncated_prefix<typename Context::G2Point, Fr, typename Context::G2CoeffCtx>(
                    b_fft,
                    std::span<const Fr>(precomputed.roots),
                    xn_b,
                    precomputed.forward_plan,
                    timings,
                    "middle_product_precompute.fft_b0");
            }
            {
                ScopedTimer timer(timings, "middle_product_precompute.fft_b1_ns");
                fft_forward_truncated_prefix<typename Context::G2Point, Fr, typename Context::G2CoeffCtx>(
                    b1_fft,
                    std::span<const Fr>(precomputed.roots),
                    xn_b,
                    precomputed.forward_plan,
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
                    b_fft[j] = Context::G2CoeffCtx::plus(b_fft[j], zeta_b1);
                }
            }

            if constexpr (Context::batch_affine_pointwise) {
                ScopedTimer timer(
                    timings,
                    "middle_product_precompute.pointwise_prepare_ns");
                std::vector<typename Context::G2AffinePoint> b_fft_affine(
                    precomputed.fft_size);
                std::vector<typename Context::G2Point::Field_type> prefix_products(
                    precomputed.fft_size);
                Context::batch_to_affine_g2(
                    std::span<typename Context::G2AffinePoint>(b_fft_affine),
                    std::span<const typename Context::G2Point>(b_fft),
                    std::span<typename Context::G2Point::Field_type>(prefix_products));
                precomputed.prepared_lines.resize(precomputed.fft_size);
                precomputed.b_fft_is_zero.resize(precomputed.fft_size);
                for (std::size_t i = 0; i < precomputed.fft_size; ++i) {
                    const bool is_zero =
                        (Context::G2AffinePoint::is_inf(b_fft_affine[i]) != 0);
                    precomputed.b_fft_is_zero[i] =
                        static_cast<uint8_t>(is_zero ? 1 : 0);
                    if (!is_zero) {
                        Context::prepare_lines_affine(
                            precomputed.prepared_lines[i],
                            b_fft_affine[i]);
                    }
                }
            } else {
                precomputed.b_fft = std::move(b_fft);
            }

            return precomputed;
        }

        template<class Context>
        inline void validate_middle_product_range_precompute(
            const MiddleProductPrecomputedBImpl<Context>& precomputed)
        {
            const bool invalid_pointwise_state = [&]() {
                if constexpr (Context::batch_affine_pointwise) {
                    return
                        precomputed.prepared_lines.size() != precomputed.fft_size ||
                        precomputed.b_fft_is_zero.size() != precomputed.fft_size;
                } else {
                    return precomputed.b_fft.size() != precomputed.fft_size;
                }
            }();
            const bool supported_n = middle_product_n_is_supported(precomputed.n);
            const std::size_t convolution_size = supported_n ? 2 * precomputed.n : 0;
            const bool invalid_fft_size =
                !fft_detail::valid_transform_size(precomputed.fft_size) ||
                precomputed.fft_size > middle_product_max_fft_size ||
                precomputed.fft_size < convolution_size ||
                (precomputed.fft_size >> 1) >= convolution_size;
            const std::size_t expected_stages = invalid_fft_size
                ? 0
                : static_cast<std::size_t>(
                    std::countr_zero(precomputed.fft_size));
            if (!supported_n ||
                invalid_fft_size ||
                invalid_pointwise_state ||
                precomputed.roots.size() != precomputed.fft_size ||
                Context::inverse_roots_size(precomputed.inv_roots) !=
                    precomputed.fft_size ||
                precomputed.forward_plan.n != precomputed.fft_size ||
                precomputed.forward_plan.xn !=
                    fft_admissible_size(precomputed.n, precomputed.fft_size) ||
                precomputed.forward_plan.log2_n != expected_stages ||
                precomputed.forward_plan.stage_work.size() !=
                    precomputed.forward_plan.log2_n ||
                precomputed.inverse_plan.n != precomputed.fft_size ||
                precomputed.inverse_plan.start != precomputed.n ||
                precomputed.inverse_plan.length != precomputed.n ||
                precomputed.inverse_plan.stages != expected_stages ||
                precomputed.inverse_plan.stage_left.size() !=
                    precomputed.inverse_plan.stages) {
                throw std::invalid_argument(
                    "invalid ranged middle-product precompute state");
            }
        }

        [[nodiscard]] inline bool middle_product_range_is_single_stage_phase(
            const MiddleProductRangePhase phase) noexcept
        {
            switch (phase) {
            case MiddleProductRangePhase::ForwardJ0:
            case MiddleProductRangePhase::ForwardNonzero:
            case MiddleProductRangePhase::InverseButterflies:
                return false;
            default:
                return true;
            }
        }

        template<class Context>
        [[nodiscard]] inline std::size_t middle_product_range_work_count_impl(
            const MiddleProductPrecomputedBImpl<Context>& precomputed,
            const MiddleProductRangePhase phase,
            const std::size_t stage)
        {
            validate_middle_product_range_precompute(precomputed);
            if (middle_product_range_is_single_stage_phase(phase) && stage != 0) {
                throw std::invalid_argument(
                    "middle-product phase does not accept a stage");
            }

            switch (phase) {
            case MiddleProductRangePhase::InputLoad:
                return precomputed.fft_size;
            case MiddleProductRangePhase::ForwardBitReverse:
                return precomputed.forward_plan.bit_reverse_swaps.size();
            case MiddleProductRangePhase::ForwardJ0:
                if (stage >= precomputed.forward_plan.stage_work.size()) {
                    throw std::invalid_argument(
                        "middle-product forward stage is out of range");
                }
                return precomputed.forward_plan.stage_work[stage].left_j0.count;
            case MiddleProductRangePhase::ForwardNonzero:
                if (stage >= precomputed.forward_plan.stage_work.size()) {
                    throw std::invalid_argument(
                        "middle-product forward stage is out of range");
                }
                return precomputed.forward_plan.stage_work[stage]
                    .left_nonzero.count;
            case MiddleProductRangePhase::AffineConversion:
                if constexpr (Context::batch_affine_pointwise) {
                    return precomputed.fft_size;
                } else {
                    return 0;
                }
            case MiddleProductRangePhase::PointwisePairing:
                return precomputed.fft_size;
            case MiddleProductRangePhase::InverseBitReverse:
                return precomputed.inverse_plan.bit_reverse_swaps.size();
            case MiddleProductRangePhase::InverseButterflies:
                if (stage >= precomputed.inverse_plan.stage_left.size()) {
                    throw std::invalid_argument(
                        "middle-product inverse stage is out of range");
                }
                return precomputed.inverse_plan.stage_left[stage].count;
            case MiddleProductRangePhase::Normalization:
            case MiddleProductRangePhase::Extraction:
                return precomputed.n;
            }
            throw std::invalid_argument("invalid middle-product range phase");
        }

        inline void validate_middle_product_work_range(
            const std::size_t work_count,
            const std::size_t begin,
            const std::size_t end)
        {
            if (begin > end || end > work_count) {
                throw std::invalid_argument(
                    "middle-product work range is out of bounds");
            }
        }

        template<class Context>
        inline void middle_product_input_range_impl(
            std::span<const typename Context::G1Point> input,
            const MiddleProductPrecomputedBImpl<Context>& precomputed,
            MiddleProductWorkspaceImpl<Context>& workspace,
            const std::size_t begin,
            const std::size_t end)
        {
            if (input.size() != precomputed.n) {
                throw std::invalid_argument(
                    "middle-product input length does not match its precompute");
            }
            workspace.require_capacity(precomputed.fft_size);
            const std::size_t work_count = middle_product_range_work_count_impl(
                precomputed, MiddleProductRangePhase::InputLoad, 0);
            validate_middle_product_work_range(work_count, begin, end);

            typename Context::G1Point* const values = workspace.a_fft.data();
            const std::size_t copy_end = std::min(end, precomputed.n);
            for (std::size_t idx = begin; idx < copy_end; ++idx) {
                values[idx] = input[idx];
            }
            const std::size_t fill_begin = std::max(begin, precomputed.n);
            for (std::size_t idx = fill_begin; idx < end; ++idx) {
                values[idx] = Context::g1_infinity();
            }
        }

        template<class Context>
        inline void middle_product_execute_range_impl(
            const MiddleProductPrecomputedBImpl<Context>& precomputed,
            MiddleProductWorkspaceImpl<Context>& workspace,
            const MiddleProductRangePhase phase,
            const std::size_t stage,
            const std::size_t begin,
            const std::size_t end)
        {
            workspace.require_capacity(precomputed.fft_size);
            const std::size_t work_count =
                middle_product_range_work_count_impl(precomputed, phase, stage);
            validate_middle_product_work_range(work_count, begin, end);
            if (phase == MiddleProductRangePhase::InputLoad ||
                phase == MiddleProductRangePhase::Extraction) {
                throw std::invalid_argument(
                    "middle-product phase requires its typed range operation");
            }
            if (begin == end) {
                return;
            }

            std::span<typename Context::G1Point> a_fft(
                workspace.a_fft.data(), precomputed.fft_size);
            std::span<typename Context::GtFftElement> c_fft(
                workspace.c_fft.data(), precomputed.fft_size);

            switch (phase) {
            case MiddleProductRangePhase::ForwardBitReverse:
                fft_apply_bit_reverse_swaps_range(
                    a_fft,
                    std::span<const FftBitReverseSwap>(
                        precomputed.forward_plan.bit_reverse_swaps),
                    begin,
                    end);
                return;
            case MiddleProductRangePhase::ForwardJ0:
                fft_forward_stage_j0_range<
                    typename Context::G1Point,
                    typename Context::G1CoeffCtx>(
                        a_fft,
                        precomputed.forward_plan,
                        stage,
                        begin,
                        end);
                return;
            case MiddleProductRangePhase::ForwardNonzero:
                fft_forward_stage_nonzero_range<
                    typename Context::G1Point,
                    Fr,
                    typename Context::G1CoeffCtx>(
                        a_fft,
                        std::span<const Fr>(precomputed.roots),
                        precomputed.forward_plan,
                        stage,
                        begin,
                        end);
                return;
            case MiddleProductRangePhase::AffineConversion:
                if constexpr (Context::batch_affine_pointwise) {
                    const std::size_t length = end - begin;
                    Context::batch_to_affine_g1(
                        std::span<typename Context::G1AffinePoint>(
                            workspace.a_fft_affine.data() + begin,
                            length),
                        std::span<const typename Context::G1Point>(
                            workspace.a_fft.data() + begin,
                            length),
                        std::span<typename Context::G1Point::Field_type>(
                            workspace.affine_prefix_products.data() + begin,
                            length));
                }
                return;
            case MiddleProductRangePhase::PointwisePairing:
                if constexpr (Context::batch_affine_pointwise) {
                    for (std::size_t idx = begin; idx < end; ++idx) {
                        if (precomputed.b_fft_is_zero[idx] != 0 ||
                            Context::G1AffinePoint::is_inf(
                                workspace.a_fft_affine[idx]) != 0) {
                            c_fft[idx] = Context::identity();
                        } else {
                            c_fft[idx] = Context::pointwise_pair_prepared(
                                workspace.a_fft_affine[idx],
                                precomputed.prepared_lines[idx]);
                        }
                    }
                } else {
                    for (std::size_t idx = begin; idx < end; ++idx) {
                        c_fft[idx] = Context::pointwise_pair_direct(
                            a_fft[idx], precomputed.b_fft[idx]);
                    }
                }
                return;
            case MiddleProductRangePhase::InverseBitReverse:
                fft_apply_bit_reverse_swaps_range(
                    c_fft,
                    std::span<const FftBitReverseSwap>(
                        precomputed.inverse_plan.bit_reverse_swaps),
                    begin,
                    end);
                return;
            case MiddleProductRangePhase::InverseButterflies:
                Context::inverse_fft_stage_range(
                    c_fft,
                    precomputed.inv_roots,
                    precomputed.inverse_plan,
                    stage,
                    begin,
                    end);
                return;
            case MiddleProductRangePhase::Normalization:
                Context::normalize_inverse_range(
                    c_fft,
                    precomputed.inv_roots,
                    precomputed.inverse_fft_size,
                    precomputed.n + begin,
                    precomputed.n + end);
                return;
            case MiddleProductRangePhase::InputLoad:
            case MiddleProductRangePhase::Extraction:
                break;
            }
            throw std::invalid_argument("invalid middle-product range phase");
        }

        template<class Context, class OutputGt>
        inline void middle_product_output_range_impl(
            const MiddleProductPrecomputedBImpl<Context>& precomputed,
            MiddleProductWorkspaceImpl<Context>& workspace,
            std::span<OutputGt> output,
            const std::size_t begin,
            const std::size_t end)
        {
            validate_middle_product_output_shape(precomputed.n, output.size());
            workspace.require_capacity(precomputed.fft_size);
            const std::size_t work_count = middle_product_range_work_count_impl(
                precomputed, MiddleProductRangePhase::Extraction, 0);
            validate_middle_product_work_range(work_count, begin, end);
            Context::extract_middle_window(
                output.subspan(begin, end - begin),
                std::span<const typename Context::GtFftElement>(
                    workspace.c_fft.data(), precomputed.fft_size),
                precomputed.n + begin,
                end - begin);
        }

        // Sequential oracle for the stateless range boundary. This deliberately
        // remains separate from `middle_product_impl`, so the established
        // monolithic transform and the range path can be differentially tested
        // and benchmarked without sharing an executor.
        template<class Context, class OutputGt>
        inline void middle_product_ranged_serial_impl(
            std::span<const typename Context::G1Point> input,
            const MiddleProductPrecomputedBImpl<Context>& precomputed,
            std::span<OutputGt> output,
            MiddleProductWorkspaceImpl<Context>& workspace)
        {
            middle_product_input_range_impl<Context>(
                input, precomputed, workspace, 0, precomputed.fft_size);

            const auto execute_full = [&](const MiddleProductRangePhase phase,
                                          const std::size_t stage = 0) {
                const std::size_t count = middle_product_range_work_count_impl(
                    precomputed, phase, stage);
                middle_product_execute_range_impl<Context>(
                    precomputed, workspace, phase, stage, 0, count);
            };

            execute_full(MiddleProductRangePhase::ForwardBitReverse);
            for (std::size_t stage = 0;
                 stage < precomputed.forward_plan.stage_work.size();
                 ++stage) {
                execute_full(MiddleProductRangePhase::ForwardJ0, stage);
                execute_full(MiddleProductRangePhase::ForwardNonzero, stage);
            }
            execute_full(MiddleProductRangePhase::AffineConversion);
            execute_full(MiddleProductRangePhase::PointwisePairing);
            execute_full(MiddleProductRangePhase::InverseBitReverse);
            for (std::size_t stage = 0;
                 stage < precomputed.inverse_plan.stage_left.size();
                 ++stage) {
                execute_full(MiddleProductRangePhase::InverseButterflies, stage);
            }
            execute_full(MiddleProductRangePhase::Normalization);
            middle_product_output_range_impl<Context>(
                precomputed, workspace, output, 0, precomputed.n);
        }

        template<class Context, class OutputGt>
        inline void middle_product_impl(
            std::span<const typename Context::G1Point> a,
            const MiddleProductPrecomputedBImpl<Context>& precomputed_b,
            std::span<OutputGt> out,
            MiddleProductWorkspaceImpl<Context>& workspace,
            TimingStats* timings)
        {
            ScopedTimer total_timer(timings, "middle_product_pre.total_ns");
            const std::size_t n = a.size();
            validate_middle_product_output_shape(n, out.size());
            if (precomputed_b.n != n) {
                throw std::invalid_argument("middle product length does not match its precompute");
            }
            const bool invalid_pointwise_state = [&]() {
                if constexpr (Context::batch_affine_pointwise) {
                    return precomputed_b.prepared_lines.size() != precomputed_b.fft_size ||
                        precomputed_b.b_fft_is_zero.size() != precomputed_b.fft_size;
                } else {
                    return precomputed_b.b_fft.size() != precomputed_b.fft_size;
                }
            }();
            if (precomputed_b.fft_size != middle_product_fft_size(n) ||
                invalid_pointwise_state ||
                precomputed_b.roots.size() != precomputed_b.fft_size ||
                Context::inverse_roots_size(precomputed_b.inv_roots) != precomputed_b.fft_size ||
                precomputed_b.forward_plan.n != precomputed_b.fft_size ||
                precomputed_b.forward_plan.xn != fft_admissible_size(n, precomputed_b.fft_size) ||
                precomputed_b.inverse_plan.n != precomputed_b.fft_size ||
                precomputed_b.inverse_plan.start != n ||
                precomputed_b.inverse_plan.length != n) {
                throw std::invalid_argument("invalid middle-product precompute state");
            }

            workspace.require_capacity(precomputed_b.fft_size);
            std::span<typename Context::G1Point> a_fft(
                workspace.a_fft.data(), precomputed_b.fft_size);
            std::span<typename Context::GtFftElement> c_fft(
                workspace.c_fft.data(), precomputed_b.fft_size);

            {
                ScopedTimer timer(timings, "middle_product_pre.copy_a_ns");
                std::copy(a.begin(), a.end(), a_fft.begin());
            }
            {
                ScopedTimer timer(timings, "middle_product_pre.fill_a_fft_ns");
                std::fill(
                    a_fft.begin() + static_cast<std::ptrdiff_t>(n),
                    a_fft.end(),
                    Context::g1_infinity());
            }
            const std::size_t xn_a = fft_admissible_size(n, precomputed_b.fft_size);
            {
                ScopedTimer timer(timings, "middle_product_pre.fft_a_ns");
                fft_forward_truncated_prefix<typename Context::G1Point, Fr, typename Context::G1CoeffCtx>(
                    a_fft,
                    std::span<const Fr>(precomputed_b.roots),
                    xn_a,
                    precomputed_b.forward_plan,
                    timings,
                    "middle_product_pre.fft_a");
            }
            {
                ScopedTimer timer(timings, "middle_product_pre.pointwise_ns");
                if constexpr (Context::batch_affine_pointwise) {
                    Context::batch_to_affine_g1(
                        std::span<typename Context::G1AffinePoint>(
                            workspace.a_fft_affine.data(),
                            precomputed_b.fft_size),
                        std::span<const typename Context::G1Point>(a_fft),
                        std::span<typename Context::G1Point::Field_type>(
                            workspace.affine_prefix_products.data(),
                            precomputed_b.fft_size));
                    for (std::size_t i = 0; i < precomputed_b.fft_size; ++i) {
                        if (precomputed_b.b_fft_is_zero[i] != 0 ||
                            Context::G1AffinePoint::is_inf(workspace.a_fft_affine[i]) != 0) {
                            c_fft[i] = Context::identity();
                            continue;
                        }
                        c_fft[i] = Context::pointwise_pair_prepared(
                            workspace.a_fft_affine[i],
                            precomputed_b.prepared_lines[i]);
                    }
                } else {
                    for (std::size_t i = 0; i < precomputed_b.fft_size; ++i) {
                        c_fft[i] = Context::pointwise_pair_direct(
                            a_fft[i], precomputed_b.b_fft[i]);
                    }
                }
            }
            {
                ScopedTimer timer(timings, "middle_product_pre.ifft_ns");
                Context::inverse_fft_middle_window(
                    c_fft,
                    precomputed_b.inv_roots,
                    n,
                    n,
                    precomputed_b.inverse_plan);
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
            validate_middle_product_output_shape(n, out.size());
            if (precomputed_b.b_tail.size() != n) {
                throw std::invalid_argument("invalid middle-product precompute state");
            }

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
            validate_middle_product_b_shape(b.size(), n);
            if constexpr (std::is_same_v<typename Context::G2Point, blst::G2Point>) {
                return precompute_middle_product_b_impl<Context>(b, n, timings);
            } else {
                const std::vector<typename Context::G2Point> converted =
                    convert_g2_inputs<Context>(b);
                return precompute_middle_product_b_impl<Context>(
                    std::span<const typename Context::G2Point>(converted),
                    n,
                    timings);
            }
        }

        template<class Context>
        inline void middle_product_public(
            std::span<const blst::G1Point> a,
            const MiddleProductPrecomputedBImpl<Context>& precomputed_b,
            std::span<GT> out,
            MiddleProductWorkspaceImpl<Context>& workspace,
            TimingStats* timings)
        {
            validate_middle_product_output_shape(a.size(), out.size());
            if (precomputed_b.n != a.size()) {
                throw std::invalid_argument("middle product length does not match its precompute");
            }
            if constexpr (std::is_same_v<typename Context::G1Point, blst::G1Point>) {
                middle_product_impl<Context>(a, precomputed_b, out, workspace, timings);
            } else {
                const std::vector<typename Context::G1Point> converted_a =
                    convert_g1_inputs<Context>(a);
                middle_product_impl<Context>(
                    std::span<const typename Context::G1Point>(converted_a),
                    precomputed_b,
                    out,
                    workspace,
                    timings);
            }
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
            validate_middle_product_b_shape(b.size(), n);
            validate_middle_product_output_shape(n, out.size());

            const auto precomputed =
                precompute_middle_product_b_public<Context>(b, n, timings);
            MiddleProductWorkspaceImpl<Context> workspace(precomputed.fft_size);
            middle_product_public<Context>(a, precomputed, out, workspace, timings);
        }

    } // namespace detail

    using MiddleProductScalarContext = detail::ScalarMiddleProductContext;
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    using MiddleProductAvxContext = detail::AvxMiddleProductContext;
    static_assert(
        MiddleProductAvxContext::batch_affine_pointwise,
        "the production AVX middle product requires prepared batch-affine pointwise pairing");
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

    /// Immutable multipair line precomputation. Fields are expert-visible but
    /// must not be modified after successful construction.
    template<class Context>
    struct MiddleProductPrecomputedBMultipairImpl {
        std::size_t n = 0;
        std::vector<typename Context::PreparedLineTable> prepared_lines;
        std::vector<uint8_t> b_is_zero;
    };

    /// Mutable, non-thread-safe multipair evaluation storage. Reuse one object
    /// sequentially to retain temporary capacities between calls.
    template<class Context>
    struct MiddleProductWorkspaceMultipairImpl {
        std::vector<typename Context::G1AffinePoint> affine_a_nonzero;
        std::vector<std::size_t> a_indices_nonzero;
    };

    template<>
    struct MiddleProductWorkspaceMultipairImpl<MiddleProductScalarContext> {
        std::vector<MiddleProductScalarContext::G1AffinePoint> affine_a_nonzero;
        std::vector<std::size_t> a_indices_nonzero;
        std::vector<scalar::PreparedG1LineScale> a_line_scales;
        std::vector<const Fp6*> line_batch_scalar;
        std::vector<scalar::PreparedG1LineScale> line_scale_batch_scalar;
    };

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    template<>
    struct MiddleProductWorkspaceMultipairImpl<MiddleProductAvxContext> {
        std::vector<MiddleProductAvxContext::G1AffinePoint> affine_a_nonzero;
        std::vector<std::size_t> a_indices_nonzero;
        std::vector<const MiddleProductAvxContext::G1AffinePoint*> p_batch;
        std::vector<const avx::PreparedLine*> line_batch_avx;
        std::vector<avx::MillerResult> miller_batch;
        std::vector<avx::GT> final_batch;
        avx::GT::FinalExpScratch final_exp_scratch;
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
        std::vector<scalar::PreparedG1LineScale>& out_scales,
        std::span<const scalar::G1AffinePoint> input)
    {
        out_scales.clear();
        out_scales.reserve(input.size());
        for (const scalar::G1AffinePoint& point : input) {
            scalar::PreparedG1LineScale scale = {};
            scalar::Fp::add(scale.neg_x2, point.x, point.x);
            scalar::Fp::neg(scale.neg_x2, scale.neg_x2);
            scalar::Fp::add(scale.y2, point.y, point.y);
            out_scales.push_back(scale);
        }
    }

    template<class Context>
    inline MiddleProductPrecomputedBMultipairImpl<Context>
    precompute_middle_product_b_multipair_backend_impl(
        std::span<const typename Context::G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr)
    {
        ScopedTimer total_timer(timings, "middle_product_multipair_precompute.total_ns");
        validate_middle_product_b_shape(b.size(), n);
        const std::size_t b_size = middle_product_b_size(n);

        MiddleProductPrecomputedBMultipairImpl<Context> precomputed;
        precomputed.n = n;
        precomputed.prepared_lines.resize(b_size);
        precomputed.b_is_zero.resize(b_size);

        for (std::size_t j = 0; j < b.size(); ++j) {
            const bool is_zero = (Context::G2Point::is_inf(b[j]) != 0);
            precomputed.b_is_zero[j] = static_cast<uint8_t>(is_zero ? 1 : 0);
            if (is_zero) {
                continue;
            }
            Context::prepare_lines(precomputed.prepared_lines[j], b[j]);
        }

        return precomputed;
    }

    template<class Context>
    inline MiddleProductPrecomputedBMultipairImpl<Context> precompute_middle_product_b_multipair_impl(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr)
    {
        validate_middle_product_b_shape(b.size(), n);
        if constexpr (std::is_same_v<typename Context::G2Point, G2Point>) {
            return precompute_middle_product_b_multipair_backend_impl<Context>(
                b,
                n,
                timings);
        } else {
            const std::vector<typename Context::G2Point> converted =
                middle_product_convert_multipair_g2_inputs<Context>(b);
            return precompute_middle_product_b_multipair_backend_impl<Context>(
                std::span<const typename Context::G2Point>(converted),
                n,
                timings);
        }
    }

    inline GT middle_product_multipair_scalar_window_dense(
        std::span<const scalar::PreparedG1LineScale> a_line_scales,
        const MiddleProductPrecomputedBMultipairScalar& precomputed_b,
        std::size_t output_index,
        std::vector<const Fp6*>& line_batch)
    {
        const std::size_t n = precomputed_b.n;
        line_batch.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            line_batch[i] = precomputed_b.prepared_lines[n + output_index - i].lines;
        }

        scalar::MillerResult miller = scalar::MillerResult::identity();
        scalar::GT::miller_loop_lines_n_preprocessed(
            miller,
            line_batch.data(),
            a_line_scales.data(),
            line_batch.size());

        scalar::GT final_pair = scalar::GT::identity();
        scalar::GT::final_exp(final_pair, miller);
        return MiddleProductScalarContext::to_public_gt(final_pair);
    }

    inline GT middle_product_multipair_scalar_window_sparse(
        std::span<const scalar::PreparedG1LineScale> a_line_scales,
        std::span<const std::size_t> a_indices_nonzero,
        const MiddleProductPrecomputedBMultipairScalar& precomputed_b,
        std::size_t output_index,
        std::vector<const Fp6*>& line_batch,
        std::vector<scalar::PreparedG1LineScale>& scale_batch)
    {
        line_batch.clear();
        scale_batch.clear();

        for (std::size_t nz = 0; nz < a_line_scales.size(); ++nz) {
            const std::size_t i = a_indices_nonzero[nz];
            const std::size_t j = precomputed_b.n + output_index - i;
            if (precomputed_b.b_is_zero[j] != 0) {
                continue;
            }
            scale_batch.push_back(a_line_scales[nz]);
            line_batch.push_back(precomputed_b.prepared_lines[j].lines);
        }

        if (line_batch.empty()) {
            return GT::identity();
        }

        scalar::MillerResult miller = scalar::MillerResult::identity();
        scalar::GT::miller_loop_lines_n_preprocessed(
            miller,
            line_batch.data(),
            scale_batch.data(),
            line_batch.size());

        scalar::GT final_pair = scalar::GT::identity();
        scalar::GT::final_exp(final_pair, miller);
        return MiddleProductScalarContext::to_public_gt(final_pair);
    }

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    inline avx::MillerResult middle_product_multipair_avx_window_dense(
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

        avx::MillerResult miller = avx::MillerResult::identity();
        avx::GT::miller_loop_lines_n_preprocessed(
            miller,
            line_batch.data(),
            a_affine_nonzero.data(),
            line_batch.size());
        return miller;
    }

    inline avx::MillerResult middle_product_multipair_avx_window_sparse(
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
            return avx::MillerResult::identity();
        }

        avx::MillerResult miller = avx::MillerResult::identity();
        avx::GT::miller_loop_lines_n(
            miller,
            line_batch.data(),
            p_batch.data(),
            line_batch.size());
        return miller;
    }

    inline void middle_product_finalize_avx_millers(
        std::span<GT> out,
        std::span<const avx::MillerResult> millers,
        std::vector<avx::GT>& final_pairs,
        avx::GT::FinalExpScratch& final_exp_scratch)
    {
        BLST_ASSERT(out.size() >= millers.size());
        final_pairs.resize(millers.size());
        avx::GT::final_exp_n(
            final_pairs.data(),
            millers.data(),
            millers.size(),
            final_exp_scratch);
        for (std::size_t i = 0; i < millers.size(); ++i) {
            out[i] = MiddleProductAvxContext::to_public_gt(avx::AvxGtElement(final_pairs[i].data));
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
        validate_middle_product_output_shape(n, out.size());
        if (precomputed_b.n != n) {
            throw std::invalid_argument("middle product length does not match its precompute");
        }
        const std::size_t b_size = middle_product_b_size(n);
        if (precomputed_b.prepared_lines.size() != b_size ||
            precomputed_b.b_is_zero.size() != b_size) {
            throw std::invalid_argument("invalid multipair middle-product precompute state");
        }

        middle_product_collect_multipair_g1_affine<Context>(
            workspace.affine_a_nonzero,
            workspace.a_indices_nonzero,
            a);

        {
            ScopedTimer timer(timings, "middle_product_multipair.pointwise_ns");
            if constexpr (std::is_same_v<Context, MiddleProductScalarContext>) {
                middle_product_prepare_scalar_line_eval_points(
                    workspace.a_line_scales,
                    std::span<const scalar::G1AffinePoint>(workspace.affine_a_nonzero));

                const bool dense_a = (workspace.a_indices_nonzero.size() == n);
                const bool dense_b = middle_product_is_dense_mask(std::span<const uint8_t>(precomputed_b.b_is_zero));

                if (dense_a && dense_b) {
                    for (std::size_t k = 0; k < n; ++k) {
                        out[k] = middle_product_multipair_scalar_window_dense(
                            std::span<const scalar::PreparedG1LineScale>(workspace.a_line_scales),
                            precomputed_b,
                            k,
                            workspace.line_batch_scalar);
                    }
                } else {
                    workspace.line_batch_scalar.reserve(
                        workspace.affine_a_nonzero.size());
                    workspace.line_scale_batch_scalar.reserve(
                        workspace.affine_a_nonzero.size());
                    for (std::size_t k = 0; k < n; ++k) {
                        out[k] = middle_product_multipair_scalar_window_sparse(
                            std::span<const scalar::PreparedG1LineScale>(workspace.a_line_scales),
                            std::span<const std::size_t>(workspace.a_indices_nonzero),
                            precomputed_b,
                            k,
                            workspace.line_batch_scalar,
                            workspace.line_scale_batch_scalar);
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
                    std::span<const avx::MillerResult>(workspace.miller_batch),
                    workspace.final_batch,
                    workspace.final_exp_scratch);
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
        validate_middle_product_output_shape(a.size(), out.size());
        if (precomputed_b.n != a.size()) {
            throw std::invalid_argument("middle product length does not match its precompute");
        }
        if constexpr (std::is_same_v<typename Context::G1Point, G1Point>) {
            middle_product_multipair_impl<Context>(
                a,
                precomputed_b,
                out,
                workspace,
                timings);
        } else {
            const std::vector<typename Context::G1Point> converted_a =
                middle_product_convert_multipair_g1_inputs<Context>(a);
            middle_product_multipair_impl<Context>(
                std::span<const typename Context::G1Point>(converted_a),
                precomputed_b,
                out,
                workspace,
                timings);
        }
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
        validate_middle_product_output_shape(n, out.size());
        const std::size_t b_size = middle_product_b_size(n);
        if (precomputed_b.prepared_lines.size() != b_size ||
            precomputed_b.b_is_zero.size() != b_size) {
            throw std::invalid_argument("invalid multipair middle-product precompute state");
        }

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
                    scalar::MillerResult miller = scalar::MillerResult::identity();
                    scalar::GT::miller_loop_lines_finite(
                        miller,
                        precomputed_b.prepared_lines[n + k].lines,
                        affine_a0);
                    scalar::GT final_pair = scalar::GT::identity();
                    scalar::GT::final_exp(final_pair, miller);
                    out[k] = MiddleProductScalarContext::to_public_gt(final_pair);
                }
            } else {
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
                std::vector<avx::MillerResult> millers(
                    n, avx::MillerResult::identity());
                std::vector<avx::GT> finals;
                avx::GT::FinalExpScratch final_exp_scratch;
                for (std::size_t k = 0; k < n; ++k) {
                    if (precomputed_b.b_is_zero[n + k] != 0) {
                        continue;
                    }
                    avx::GT::miller_loop_lines_finite(
                        millers[k],
                        precomputed_b.prepared_lines[n + k].lines,
                        affine_a0);
                }
                middle_product_finalize_avx_millers(
                    out.first(n),
                    std::span<const avx::MillerResult>(millers),
                    finals,
                    final_exp_scratch);
#endif
            }
        }
    }

    // Keep the workspace entry point separate from the convenience
    // implementation. This preserves the measured layout of the outlined AVX
    // pairing and final-exponentiation kernels while allowing hot callers to
    // retain all temporary capacity.
    template<class Context>
    inline void update_multipair_workspace_impl(
        const typename Context::G1Point& a0,
        const MiddleProductPrecomputedBMultipairImpl<Context>& precomputed_b,
        std::span<GT> out,
        MiddleProductWorkspaceMultipairImpl<Context>& workspace,
        TimingStats* timings = nullptr)
    {
        if constexpr (std::is_same_v<Context, MiddleProductScalarContext>) {
            static_cast<void>(workspace);
            update_multipair_impl<Context>(a0, precomputed_b, out, timings);
        } else {
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
            ScopedTimer total_timer(timings, "middle_product_multipair_update.total_ns");
            const std::size_t n = precomputed_b.n;
            validate_middle_product_output_shape(n, out.size());
            const std::size_t b_size = middle_product_b_size(n);
            if (precomputed_b.prepared_lines.size() != b_size ||
                precomputed_b.b_is_zero.size() != b_size) {
                throw std::invalid_argument("invalid multipair middle-product precompute state");
            }

            if (Context::G1Point::is_inf(a0) != 0) {
                std::fill(out.begin(), out.begin() + n, GT::identity());
                return;
            }

            typename Context::G1AffinePoint affine_a0 = {};
            Context::G1Point::to_affine(affine_a0, a0);

            ScopedTimer timer(timings, "middle_product_multipair_update.pointwise_ns");
            workspace.miller_batch.resize(n);
            for (std::size_t k = 0; k < n; ++k) {
                if (precomputed_b.b_is_zero[n + k] != 0) {
                    workspace.miller_batch[k] = avx::MillerResult::identity();
                    continue;
                }
                avx::GT::miller_loop_lines_finite(
                    workspace.miller_batch[k],
                    precomputed_b.prepared_lines[n + k].lines,
                    affine_a0);
            }
            middle_product_finalize_avx_millers(
                out.first(n),
                std::span<const avx::MillerResult>(workspace.miller_batch),
                workspace.final_batch,
                workspace.final_exp_scratch);
#endif
        }
    }

    /**
     * @brief Precompute a length-`2*n` G2 operand for repeated evaluation.
     * @pre `1 <= n <= middle_product_max_n` and `b.size() == 2*n`.
     * @pre Every element of @p b is a reduced G2 subgroup point.
     * @post The returned object is immutable and records logical length @p n.
     * @throws std::invalid_argument for an invalid shape, std::length_error for
     * an unsupported length, and std::bad_alloc on allocation failure.
     */
    inline MiddleProductPrecomputedBScalar precompute_middle_product_b_scalar(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr)
    {
        validate_middle_product_b_shape(b.size(), n);
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
        validate_middle_product_b_shape(b.size(), n);
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
        validate_middle_product_b_shape(b.size(), n);
        return precompute_middle_product_b_multipair_impl<MiddleProductScalarContext>(b, n, timings);
    }

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    inline MiddleProductPrecomputedBMultipairAvx precompute_middle_product_b_multipair_avx(
        std::span<const G2Point> b,
        std::size_t n,
        TimingStats* timings = nullptr)
    {
        validate_middle_product_b_shape(b.size(), n);
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
        validate_middle_product_b_shape(b.size(), n);
        validate_middle_product_output_shape(n, out.size());

        const std::size_t fft_size = middle_product_basic_fft_size(n);
        const std::size_t b_size = middle_product_b_size(n);
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
        for (std::size_t i = 0; i < b_size; ++i) {
            b_fft[i] = b[i];
        }

        const std::size_t xn_a = fft_admissible_size(n, fft_size);
        const std::size_t xn_b = fft_admissible_size(b_size, fft_size);
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

    /**
     * @brief Evaluate a precomputed pairing middle product using reusable
     * caller-owned workspace.
     * @pre `a.size() == precomputed_b.n` and `out.size() >= a.size()`.
     * @pre `workspace.capacity() >= precomputed_b.fft_size`.
     * @pre Every input is a reduced G1 subgroup point and output is disjoint
     * from all inputs. No Release membership validation is performed.
     * @post The first `a.size()` outputs contain GT elements.
     */
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

    /**
     * Benchmark-selected one-shot dispatch for the scalar backend.
     *
     * The FFT cost changes discontinuously when `2*n` crosses a power of two,
     * so a single monotone cutoff is not representative. These ranges are the
     * measured multipair winners on the supported Peach deployment CPU.
     */
    [[nodiscard]] inline constexpr bool middle_product_auto_uses_multipair_scalar(
        std::size_t n) noexcept
    {
        if (n <= 32) {
            return n <= 31;
        }
        if (n <= 64) {
            return n <= 49;
        }
        if (n <= 128) {
            return n <= 76;
        }
        return false;
    }

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    /** Benchmark-selected one-shot dispatch for the native AVX backend. */
    [[nodiscard]] inline constexpr bool middle_product_auto_uses_multipair_avx(
        std::size_t n) noexcept
    {
        if (n <= 32) {
            return n <= 31;
        }
        if (n <= 64) {
            return n <= 50;
        }
        if (n <= 128) {
            return n <= 76;
        }
        return false;
    }
#endif

    [[nodiscard]] inline constexpr bool middle_product_auto_uses_multipair(
        std::size_t n) noexcept
    {
#if BLST_DEFAULT_AVX
        return middle_product_auto_uses_multipair_avx(n);
#else
        return middle_product_auto_uses_multipair_scalar(n);
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
            if (middle_product_auto_uses_multipair_scalar(a.size())) {
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
            if (middle_product_auto_uses_multipair_avx(a.size())) {
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
            if (middle_product_auto_uses_multipair(a.size())) {
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
        MiddleProductWorkspaceScalar workspace(precomputed_b.fft_size);
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
        MiddleProductWorkspaceAvx workspace(precomputed_b.fft_size);
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
        MiddleProductWorkspace workspace(precomputed_b.fft_size);
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

    /**
     * @brief Evaluate the `a[0]` contribution against a precomputed B tail.
     * @pre @p a0 is a reduced G1 subgroup point and @p out has at least
     * `precomputed_b.n` disjoint elements.
     * @post The first `n` outputs are `e(a0, b[n+k])` in GT.
     */
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
        MiddleProductWorkspaceMultipairScalar& workspace,
        TimingStats* timings = nullptr)
    {
        update_multipair_workspace_impl<MiddleProductScalarContext>(
            MiddleProductScalarContext::to_backend_g1(a0),
            precomputed_b,
            out,
            workspace,
            timings);
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
        MiddleProductWorkspaceMultipairAvx& workspace,
        TimingStats* timings = nullptr)
    {
        update_multipair_workspace_impl<MiddleProductAvxContext>(
            MiddleProductAvxContext::to_backend_g1(a0),
            precomputed_b,
            out,
            workspace,
            timings);
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
        MiddleProductWorkspaceMultipair& workspace,
        TimingStats* timings = nullptr)
    {
#if BLST_DEFAULT_AVX
        update_multipair_avx(a0, precomputed_b, out, workspace, timings);
#else
        update_multipair_scalar(a0, precomputed_b, out, workspace, timings);
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

    /**
     * @brief Combine two GT middle-product vectors elementwise.
     * @pre Inputs contain reduced GT elements, sizes match, and output has
     * sufficient length. Exact output aliasing is supported; partial overlap
     * is forbidden. No Release membership validation is performed.
     */
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

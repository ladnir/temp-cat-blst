/*
 * Copyright Supranational LLC
 * Copyright 2026 Simon Langowski
 * Copyright 2026 Category Labs
 * SPDX-License-Identifier: Apache-2.0 AND MIT
 *
 * CatBLST AVX target-group and pairing owner. The public group contracts and
 * cubed pairing convention follow BLST. The RNS field, prepared-line, and
 * Miller arithmetic is adapted from VROOM commits d6ecd534 and 26bd9aee.
 */
#pragma once

#include "blst/avx/fields/conversion_inversion.h"
#include "blst/avx/fields/final_exponentiation.h"
#include "blst/avx/fields/matrixnok_fp12.h"
#include "blst/avx/curve/point.h"
#include "blst/avx/gt/glv4.h"
#include "blst/avx/pairing/matrixnok_pairing.h"
#include "blst/common/bytes.h"
#include "blst/scalar/curve/GT.h"

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace blst::avx {

namespace avx_gt_detail {

inline constexpr size_t kMillerLoopNMax = 16;

using Ring = matrixnok::FpRing;
using FP12 = matrixnok::Fp12Ring;
using PairingFp2 = matrixnok::G2Fp2Ring;
using PairingFp2Element = PairingFp2::StandardElement;
using PairingFp2Prepared = FP2N<typename Ring::StandardElement, typename Ring::StandardElement, typename Ring::StandardElement>;
using PairingFp2LeftPrepared = decltype(std::declval<PairingFp2Element>().prep_left(std::declval<const Ring&>()));

inline const Ring& ring()
{
    return avx_curve_detail::kCurveRing;
}

inline const PairingFp2& fp2ring()
{
    return avx_curve_detail::kCurveFp2Ring;
}

inline constexpr FP12 kFp12Ring{};

inline const FP12& fp12ring()
{
    return kFp12Ring;
}

inline std::array<uint64_t, 4> scalar_to_u64_le(const byte* scalar, size_t nbits)
{
    BLST_ASSERT(nbits <= 256);
    std::array<uint64_t, 4> out = {};
    const size_t nbytes = (nbits + 7) / 8;
    for (size_t i = 0; i < nbytes; ++i) {
        out[i / 8] |= static_cast<uint64_t>(scalar[i]) << ((i % 8) * 8);
    }
    const size_t top_limb_bits = nbits & 63U;
    if (top_limb_bits != 0) {
        out[nbits / 64] &= (uint64_t{1} << top_limb_bits) - 1U;
    }
    return out;
}

inline blst::Fp12 to_scalar(const typename FP12::StandardElement& value)
{
    return avx_blst_direct_bridge::ring_fp12_to_scalar(ring(), value);
}

inline typename FP12::StandardElement from_scalar(const blst::Fp12& value)
{
    return avx_blst_direct_bridge::scalar_fp12_to_ring(ring(), value);
}

// The established CatBLST GT wire order interleaves the two Fp6 halves for
// each Fp2 coefficient. Keep this order stable while canonicalizing each
// native MatrixNoK coefficient independently at the serialization boundary.
inline constexpr std::array<size_t, 12> kGtWireToNativeCoefficient = {
    0, 1, 6, 7, 2, 3, 8, 9, 4, 5, 10, 11,
};

inline void fp_to_bendian(
    unsigned char out[48], const Ring::StandardElement& value)
{
    blst::Fp montgomery =
        avx_blst_direct_bridge::ring_to_scalar_fp(ring(), value);
    blst::Fp canonical = {};
    blst::Fp::from_mont(canonical, montgomery);
    blst::be_bytes_from_limbs(out, canonical.l, sizeof(canonical));
}

[[nodiscard]] inline bool fp_from_bendian(
    Ring::StandardElement& out, const unsigned char in[48])
{
    blst::Fp canonical = {};
    blst::limbs_from_be_bytes(canonical.l, in, sizeof(canonical));
    if (!blst::Fp::is_redc(canonical))
        return false;
    out = Ring::from_canonical_limbs(canonical);
    return true;
}

inline bool_t fp12_equal(
    const FP12::StandardElement& a,
    const FP12::StandardElement& b)
{
    bool_t equal = 1;
    BLST_UNROLL
    for (size_t i = 0; i < a.size(); ++i) {
        const blst::Fp a_canonical =
            avx_blst_direct_bridge::ring_to_scalar_fp(ring(), a[i]);
        const blst::Fp b_canonical =
            avx_blst_direct_bridge::ring_to_scalar_fp(ring(), b[i]);
        equal &= blst::Fp::equal(a_canonical, b_canonical);
    }
    return equal;
}

inline bool_t fp12_is_zero(const FP12::StandardElement& value)
{
    return fp12_equal(value, FP12::StandardElement{});
}

constexpr AVXVector<8> literal_vector(const std::array<uint64_t, 8>& lanes)
{
    return AVXVector<8>(_vv512_const_u64(
        lanes[0], lanes[1], lanes[2], lanes[3],
        lanes[4], lanes[5], lanes[6], lanes[7]));
}

constexpr Ring::StandardElement literal_standard_element(
    const std::array<std::array<uint64_t, 8>, 2>& value)
{
    return Ring::StandardElement(
        Ring::StandardBoundedElement(literal_vector(value[0])),
        Ring::StandardBoundedElement(literal_vector(value[1])));
}

constexpr FP12::StandardElement gt_generator_value()
{
    const auto& value = matrixnok::generated::gt_generator;
    return FP12::StandardElement{{
        literal_standard_element(value[0]),
        literal_standard_element(value[1]),
        literal_standard_element(value[2]),
        literal_standard_element(value[3]),
        literal_standard_element(value[4]),
        literal_standard_element(value[5]),
        literal_standard_element(value[6]),
        literal_standard_element(value[7]),
        literal_standard_element(value[8]),
        literal_standard_element(value[9]),
        literal_standard_element(value[10]),
        literal_standard_element(value[11]),
    }};
}

constexpr FP12::StandardElement fp12_identity_value()
{
    return FP12::StandardElement{{
        literal_standard_element({{
            matrixnok::generated::one_m1,
            matrixnok::generated::one_m2,
        }}),
        {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {},
    }};
}

template<int64_t N>
BLST_FORCEINLINE typename FP12::StandardElement narrow_to_standard(const typename FP12::template CommonElement<N>& value)
{
    return fp12ring().narrow_to_standard(value, ring());
}

inline PairingFp2Prepared prep_left_helper(const PairingFp2Element& point)
{
    return PairingFp2Prepared(point.x(), point.y(), ring().standard_negate(point.y()));
}

template<class T>
inline PairingFp2Element normalize_fp2(const T& value)
{
    return fp2ring().batch_reduce_expand(value.prep_left(ring()).mul_fp(ring().one()))[0];
}

template<class T>
inline T ctime_select_value(bool_t selection, const T& when_true, const T& when_false)
{
    return when_true.ctime_select(selection, when_false);
}

template<class T, size_t N>
inline std::array<T, N> ctime_select_value(
    bool_t selection,
    const std::array<T, N>& when_true,
    const std::array<T, N>& when_false)
{
    std::array<T, N> out = {};
    for (size_t i = 0; i < N; ++i) {
        out[i] = ctime_select_value(selection, when_true[i], when_false[i]);
    }
    return out;
}

struct PreparedLine
{
    // Default initialization is overwrite-before-read. Value initialization
    // (`PreparedLine{}`) creates initialized zero carriers without reductions.
    // Raw zero storage is not an encoded identity line: a producer must
    // overwrite all three members before evaluation.
    PairingFp2Element c;
    PairingFp2LeftPrepared x;
    PairingFp2LeftPrepared y;
};

struct alignas(64) PreparedLineTable
{
    // A default table is only overwrite-before-read storage. `precompute_lines`
    // fills every entry, including the encoded pairing-identity table.
    PreparedLine lines[68];
};

inline typename FP12::StandardElement make_line_value(
    const PairingFp2Element& c,
    const PairingFp2Element& x,
    const PairingFp2Element& y)
{
    return {
        c.x(),
        c.y(),
        x.x(),
        x.y(),
        ring().zero(),
        ring().zero(),
        ring().zero(),
        ring().zero(),
        y.x(),
        y.y(),
        ring().zero(),
        ring().zero(),
    };
}

inline auto scale_line_xy(const PreparedLine& in, const G1AffinePoint& p)
{
    return fp2ring().batch_reduce_expand(
        in.x.mul_fp(p.X),
        in.y.mul_fp(p.Y));
}

inline PreparedLine identity_prepared_line()
{
    const PairingFp2Element zero = fp2ring().zero();
    const PairingFp2LeftPrepared prepared_zero = zero.prep_left(ring());
    return PreparedLine{fp2ring().one(), prepared_zero, prepared_zero};
}

inline void apply_prepared_line(
    typename FP12::StandardElement& data,
    const PreparedLine& line,
    const G1AffinePoint& p)
{
    const auto scaled = scale_line_xy(line, p);
    fp12ring().flat_mulxy00z0(data, data, line.c, scaled[0], scaled[1], ring());
}

inline void line_add(PreparedLine& line, G2Point& T, const G2Point& R, const G2AffinePoint& Q)
{
    static constexpr auto kThree = std::integral_constant<int, 3>{};

    const PairingFp2Prepared qx = prep_left_helper(Q.X);
    const PairingFp2Prepared qy = prep_left_helper(Q.Y);

    const auto x2x1 = qx * R.X;
    const auto x2y1 = qx * R.Y;
    const auto x2z1 = qx * R.Z;
    const auto y2x1 = qy * R.X;
    const auto y2y1 = qy * R.Y;
    const auto y2z1 = qy * R.Z;

    const auto [x2x1_reduced, x2y1_reduced, x2z1_reduced, y2x1_reduced, y2y1_reduced, y2z1_reduced] =
        fp2ring().batch_reduce_expand(x2x1, x2y1, x2z1, y2x1, y2y1, y2z1);

    line.c = normalize_fp2(x2y1_reduced - y2x1_reduced);
    line.x = normalize_fp2(y2z1_reduced - R.Y).prep_left(ring());
    line.y = normalize_fp2(R.X - x2z1_reduced).prep_left(ring());

    auto x2y1_plus_y2x1 = (x2y1_reduced + y2x1_reduced).prep_left(ring());
    const auto y2z1_plus_y1 = (y2z1_reduced + R.Y).prep(ring());
    const auto x2z1_plus_x1 = x2z1_reduced + R.X;
    const auto three_x2x1 = (x2x1_reduced * kThree).prep(ring());
    const auto b3z1 = fp2ring().mul_3b(R.Z).prep(ring());
    auto y2y1_plus_b3z1 = (y2y1_reduced + b3z1).prep_left(ring());
    const auto y2y1_minus_b3z1 = (y2y1_reduced - b3z1).prep(ring());
    auto b3_x2z1_plus_x1 = fp2ring().mul_3b(x2z1_plus_x1).prep_left(ring());

    const auto t2_19 = x2y1_plus_y2x1 * y2y1_minus_b3z1;
    const auto x3_20 = t2_19.mul_sub_accumulate(b3_x2z1_plus_x1, y2z1_plus_y1, ring());
    const auto y3_21 = b3_x2z1_plus_x1 * three_x2x1;
    const auto y3_23 = y3_21.mul_accumulate(y2y1_plus_b3z1, y2y1_minus_b3z1);
    const auto z3_25 = y2y1_plus_b3z1 * y2z1_plus_y1;
    const auto z3_26 = z3_25.mul_accumulate(x2y1_plus_y2x1, three_x2x1);

    const auto [x3, y3, z3] = fp2ring().batch_reduce_expand(x3_20, y3_23, z3_26);
    T.X = x3;
    T.Y = y3;
    T.Z = z3;
}

inline void line_dbl(PreparedLine& line, G2Point& T, const G2Point& Q)
{
    static constexpr auto kTwo = std::integral_constant<int, 2>{};
    static constexpr auto kThree = std::integral_constant<int, 3>{};
    static constexpr auto kEight = std::integral_constant<int, 8>{};

    const auto px = Q.X.prep_left(ring());
    const auto py = Q.Y.prep_left(ring());
    const auto pz = Q.Z.prep_left(ring());

    const auto yy = py.square();
    const auto zz = pz.square();
    const auto xx = px.square();
    const auto xy = px * py.to_fp2();
    const auto yz = py * pz.to_fp2();

    const auto [xx_reduced, yy_reduced, zz_reduced, xy_reduced, yz_reduced] =
        fp2ring().batch_reduce_expand(xx, yy, zz, xy, yz);

    const auto eight_yy = (yy_reduced * kEight).prep_left(ring());
    const auto b3_zz = fp2ring().mul_3b(zz_reduced).prep(ring());
    const auto yy_plus_b3_zz = (yy_reduced + b3_zz).prep(ring());
    const auto b9_zz = b3_zz * kThree;
    const auto yy_minus_b9_zz = (yy_reduced - b9_zz).prep_left(ring());
    const auto two_xy = (xy_reduced * kTwo).prep(ring());

    line.c = normalize_fp2(b3_zz - yy_reduced);
    line.x = normalize_fp2(xx_reduced * kThree).prep_left(ring());
    line.y = fp2ring().standard_negate(normalize_fp2(yz_reduced * kTwo)).prep_left(ring());

    const auto y3_14 = yy_minus_b9_zz * yy_plus_b3_zz;
    const auto x3_18 = yy_minus_b9_zz * two_xy;
    const auto z3_10 = eight_yy * yz_reduced;
    const auto y3_15 = y3_14.mul_accumulate(eight_yy, b3_zz);

    const auto [x3, y3, z3] = fp2ring().batch_reduce_expand(x3_18, y3_15, z3_10);
    T.X = x3;
    T.Y = y3;
    T.Z = z3;
}

inline void pre_add_n_dbl(PreparedLine lines[], G2Point& t, const G2AffinePoint& q, size_t n)
{
    line_add(lines[0], t, t, q);
    ++lines;
    while (n-- != 0) {
        line_dbl(lines[0], t, t);
        ++lines;
    }
}

inline void apply_add_n_dbl(
    typename FP12::StandardElement& data,
    const PreparedLine lines[],
    const G1AffinePoint& p,
    size_t n)
{
    apply_prepared_line(data, lines[0], p);
    ++lines;

    while (n-- != 0) {
        data = fp12ring().fp12_flat_sqr(data, ring());
        apply_prepared_line(data, lines[0], p);
        ++lines;
    }
}

inline void start_lines_n(
    typename FP12::StandardElement& data,
    const PreparedLine* const Qlines[],
    const G1AffinePoint P[],
    size_t n)
{
    BLST_ASSERT(n != 0);

    const auto first = scale_line_xy(Qlines[0][0], P[0]);
    data = make_line_value(Qlines[0][0].c, first[0], first[1]);

    for (size_t i = 1; i < n; ++i) {
        apply_prepared_line(data, Qlines[i][0], P[i]);
    }
}

inline void start_lines_n(
    typename FP12::StandardElement& data,
    const PreparedLine* const Qlines[],
    const G1AffinePoint* const P[],
    size_t n)
{
    BLST_ASSERT(n != 0);

    const auto first = scale_line_xy(Qlines[0][0], *P[0]);
    data = make_line_value(Qlines[0][0].c, first[0], first[1]);

    for (size_t i = 1; i < n; ++i) {
        apply_prepared_line(data, Qlines[i][0], *P[i]);
    }
}

inline void post_add_n_dbl_n(
    typename FP12::StandardElement& data,
    const PreparedLine* const Qlines[],
    const G1AffinePoint P[],
    size_t n,
    size_t start,
    size_t k)
{
    size_t line_idx = start;

    for (size_t i = 0; i < n; ++i) {
        apply_prepared_line(data, Qlines[i][line_idx], P[i]);
    }

    while (k-- != 0) {
        ++line_idx;
        data = fp12ring().fp12_flat_sqr(data, ring());
        for (size_t i = 0; i < n; ++i) {
            apply_prepared_line(data, Qlines[i][line_idx], P[i]);
        }
    }
}

inline void post_add_n_dbl_n(
    typename FP12::StandardElement& data,
    const PreparedLine* const Qlines[],
    const G1AffinePoint* const P[],
    size_t n,
    size_t start,
    size_t k)
{
    size_t line_idx = start;

    for (size_t i = 0; i < n; ++i) {
        apply_prepared_line(data, Qlines[i][line_idx], *P[i]);
    }

    while (k-- != 0) {
        ++line_idx;
        data = fp12ring().fp12_flat_sqr(data, ring());
        for (size_t i = 0; i < n; ++i) {
            apply_prepared_line(data, Qlines[i][line_idx], *P[i]);
        }
    }
}

inline void start_dbl_n(
    typename FP12::StandardElement& data,
    G2Point T[],
    const G1AffinePoint P[],
    size_t n)
{
    BLST_ASSERT(n != 0);

    PreparedLine line = {};
    line_dbl(line, T[0], T[0]);
    const auto first = scale_line_xy(line, P[0]);
    data = make_line_value(line.c, first[0], first[1]);

    for (size_t i = 1; i < n; ++i) {
        line_dbl(line, T[i], T[i]);
        apply_prepared_line(data, line, P[i]);
    }
}

inline void add_n_dbl_n(
    typename FP12::StandardElement& data,
    G2Point T[],
    const G2AffinePoint Q[],
    const G1AffinePoint P[],
    size_t n,
    size_t k)
{
    PreparedLine line = {};

    for (size_t i = 0; i < n; ++i) {
        line_add(line, T[i], T[i], Q[i]);
        apply_prepared_line(data, line, P[i]);
    }

    while (k-- != 0) {
        data = fp12ring().fp12_flat_sqr(data, ring());
        for (size_t i = 0; i < n; ++i) {
            line_dbl(line, T[i], T[i]);
            apply_prepared_line(data, line, P[i]);
        }
    }
}

}  // namespace avx_gt_detail

using PreparedLine = avx_gt_detail::PreparedLine;
using PreparedLineTable = avx_gt_detail::PreparedLineTable;

/**
 * Native AVX pairing accumulator before final exponentiation.
 *
 * This is not a GT element. Only multiplication, conjugation, and the
 * multiplicative identity are exposed because those are the operations used
 * to combine pairing equations before `GT::final_exp`.
 */
struct MillerResult
{
    using Field_type = avx_gt_detail::FP12::StandardElement;

    Field_type data;

    constexpr MillerResult() = default;
    constexpr explicit MillerResult(const Field_type& value) : data(value) {}

    /** Return the multiplicative identity in the pre-final Miller domain. */
    static constexpr MillerResult identity()
    {
        return MillerResult(avx_gt_detail::fp12_identity_value());
    }

    /** Multiply two Miller-domain values. Exact output aliasing is supported. */
    static void mul(MillerResult& out, const MillerResult& a, const MillerResult& b)
    {
        out.data = avx_gt_detail::fp12ring().fp12_flat_mul(
            a.data, b.data, avx_gt_detail::ring());
    }

    /** Apply Fp12 conjugation in place. */
    static void conjugate(MillerResult& value)
    {
        value.data = avx_gt_detail::fp12ring().fp12_flat_conjugate(
            value.data, avx_gt_detail::ring());
    }
};

/**
 * @brief Native AVX target-group representation.
 *
 * Storage does not enforce membership. Group operations assume their operands
 * already lie in the prime-order target group; Release validation occurs only
 * in explicitly named predicates. Ordinary operations support exact output
 * aliasing and forbid partial overlap.
 */
struct GT
{
    using Field_type = avx_gt_detail::FP12::StandardElement;

    /**
     * Caller-owned storage for batch final exponentiation. Its six structure-
     * of-arrays buffers retain capacity across calls and are not thread-safe.
     */
    struct FinalExpScratch
    {
        /** Establish capacity for `n` elements; may throw `std::bad_alloc`. */
        void reserve(size_t n)
        {
            f6_terms_.reserve(n);
            f2_c0_c1_c2_terms_.reserve(n);
            f2_terms_.reserve(n);
            f1_terms_.reserve(n);
            prefix_products_.reserve(n);
            f1_inverses_.reserve(n);
        }

        /** Return true when a batch of `n` needs no scratch growth. */
        [[nodiscard]] bool has_capacity(size_t n) const noexcept
        {
            return f6_terms_.capacity() >= n &&
                f2_c0_c1_c2_terms_.capacity() >= n &&
                f2_terms_.capacity() >= n &&
                f1_terms_.capacity() >= n &&
                prefix_products_.capacity() >= n &&
                f1_inverses_.capacity() >= n;
        }

    private:
        friend struct GT;

        using FieldElement = avx_gt_detail::Ring::StandardElement;
        using Fp6Element = std::array<FieldElement, 6>;
        using Fp2Element = std::array<FieldElement, 2>;

        std::vector<Fp6Element> f6_terms_;
        std::vector<Fp6Element> f2_c0_c1_c2_terms_;
        std::vector<Fp2Element> f2_terms_;
        std::vector<FieldElement> f1_terms_;
        std::vector<FieldElement> prefix_products_;
        std::vector<FieldElement> f1_inverses_;

        template<class T>
        static void grow_to(std::vector<T>& values, size_t n)
        {
            if (values.size() < n) {
                values.resize(n);
            }
        }

        void ensure_size(size_t n)
        {
            grow_to(f6_terms_, n);
            grow_to(f2_c0_c1_c2_terms_, n);
            grow_to(f2_terms_, n);
            grow_to(f1_terms_, n);
            grow_to(prefix_products_, n);
            grow_to(f1_inverses_, n);
        }
    };

    Field_type data = {};

    constexpr GT() = default;
    explicit constexpr GT(const Field_type& value) : data(value) {}

    auto operator[](size_t idx) -> decltype(data[idx]) { return data[idx]; }
    auto operator[](size_t idx) const -> decltype(data[idx]) { return data[idx]; }

    /** GT group addition, represented as Fp12 multiplication. */
    static void add(GT& out, const GT& a, const GT& b)
    {
        out.data = avx_gt_detail::fp12ring().fp12_flat_mul(a.data, b.data, avx_gt_detail::ring());
    }

    /** Compatibility spelling for `add`; no mixed-addition precondition. */
    static void dadd(GT& out, const GT& a, const GT& b)
    {
        add(out, a, b);
    }

    /** GT group subtraction, represented as multiplication by the conjugate. */
    static void sub(GT& out, const GT& a, const GT& b)
    {
        const auto b_inv = avx_gt_detail::fp12ring().fp12_flat_conjugate(b.data, avx_gt_detail::ring());
        out.data = avx_gt_detail::fp12ring().fp12_flat_mul(a.data, b_inv, avx_gt_detail::ring());
    }

    /** GT group doubling. Exact output aliasing is supported. */
    static BLST_FORCEINLINE void doubling(GT& out, const GT& a)
    {
        // The public GT representation is standard. Internal callers that can
        // keep CommonElement<6> should prefer fp12_flat_cyclotomic_sqr and
        // delay normalization, as the GLV exponentiation path does.
        add(out, a, a);
    }

    /** Replace `out` with `[2^n](out + a)`. `n` is public. */
    static void add_n_dbl(GT& out, const GT& a, size_t n)
    {
        add(out, out, a);
        for (size_t i = 0; i < n; ++i) {
            doubling(out, out);
        }
    }

    /// Constant-time GT scalar multiplication for canonical scalars at fixed
    /// public nbits, where `nbits <= 256`. `scalar` points to
    /// `ceil(nbits/8)` readable little-endian
    /// bytes and may be null only when nbits is zero. Noncanonical inputs are
    /// accepted modulo r, but timing may reveal that reduction was required.
    /// The input must already lie in GT; Release does not validate it. Exact
    /// out/a aliasing is supported.
    static void mult(GT& out, const GT& a, const byte* scalar, size_t nbits)
    {
        BLST_ASSERT(nbits <= 256);
        BLST_ASSERT(nbits == 0 || scalar != nullptr);
        BLST_ASSERT(in_group(a));
        auto words = avx_gt_detail::scalar_to_u64_le(scalar, nbits);
        out.data = avx_cyclotomic_exp_ct(
            a.data, words.data(), nbits, avx_gt_detail::fp12ring(), avx_gt_detail::ring());
        blst::vec_zero(words.data(), sizeof(words));
    }

    /// Optimized variable-time GT multiplication for public scalars only.
    /// `nbits` must not exceed 256. `scalar` points to `ceil(nbits/8)` readable
    /// little-endian bytes and may be null only when nbits is zero. Inputs are
    /// accepted modulo r. The group
    /// operand must already lie in GT; Release does not validate it. Exact
    /// out/a aliasing is supported.
    static void mult_vartime(GT& out, const GT& a, const byte* scalar, size_t nbits)
    {
        BLST_ASSERT(nbits <= 256);
        BLST_ASSERT(nbits == 0 || scalar != nullptr);
        BLST_ASSERT(in_group(a));
        auto words = avx_gt_detail::scalar_to_u64_le(scalar, nbits);
        if (nbits > 128 && nbits <= 256) {
            out.data = avx_cyclotomic_exp_glv4_unsigned_vartime(
                a.data, words.data(), avx_gt_detail::fp12ring(), avx_gt_detail::ring());
            return;
        }
        out.data = avx_cyclotomic_exp_basic(a.data, words.data(), avx_gt_detail::fp12ring(), avx_gt_detail::ring());
    }

    /** Constant-time conditional negation. cbit must be zero or one. */
    static void cneg(GT& p, limb_t cbit)
    {
        BLST_ASSERT(cbit <= 1);
        for (size_t i = 6; i < p.data.size(); ++i) {
            const auto neg = avx_gt_detail::ring().standard_negate(p.data[i]);
            p.data[i] = avx_gt_detail::ctime_select_value(cbit, neg, p.data[i]);
        }
    }

    /** GT negation. `out` may exactly alias `p`. */
    static void neg(GT& out, const GT& p)
    {
        out.data = avx_gt_detail::fp12ring().fp12_flat_conjugate(p.data, avx_gt_detail::ring());
    }

    /** Constant-time conditional swap. `cbit` must be zero or one. */
    static void cswap(GT& a, GT& b, bool_t cbit)
    {
        BLST_ASSERT(cbit <= 1);
        const GT aa = a;
        const GT bb = b;
        a.data = avx_gt_detail::ctime_select_value(cbit, bb.data, aa.data);
        b.data = avx_gt_detail::ctime_select_value(cbit, aa.data, bb.data);
    }

    /** Constant-time conditional copy. `cbit` must be zero or one. */
    static void ccopy(GT& a, const GT& b, bool_t cbit)
    {
        BLST_ASSERT(cbit <= 1);
        a.data = avx_gt_detail::ctime_select_value(cbit, b.data, a.data);
    }

    /** Return one exactly when `a` and `b` represent the same GT element. */
    static limb_t equal(const GT& a, const GT& b)
    {
        return avx_gt_detail::fp12_equal(a.data, b.data);
    }

    /** Return the canonical generator of the prime-order target group. */
    static const GT& generator();

    /** Return the multiplicative identity of GT. */
    static constexpr GT identity()
    {
        return GT(avx_gt_detail::fp12_identity_value());
    }

    /** Return one exactly when `p` is the canonical GT identity. */
    static int is_identity(const GT& p)
    {
        return equal(p, identity());
    }

    /** Perform a complete prime-order target-group membership check. */
    static bool_t in_group(const GT& p)
    {
        if (avx_gt_detail::fp12_is_zero(p.data) || !is_cyclotomic(p))
            return 0;

        const auto frobenius =
            avx_gt_detail::fp12ring().fp12_flat_frobenius_map1(
                p.data, avx_gt_detail::ring());
        const auto raised = avx_gt_detail::narrow_to_standard(
            raise_to_z(
                p.data, avx_gt_detail::fp12ring(), avx_gt_detail::ring()));
        return avx_gt_detail::fp12_equal(frobenius, raised);
    }

    /** Check the wider cyclotomic relation; this does not establish GT. */
    static bool_t is_cyclotomic(const GT& p)
    {
        const auto frobenius2 =
            avx_gt_detail::fp12ring().fp12_flat_frobenius_map2(
                p.data, avx_gt_detail::ring());
        const auto frobenius4 =
            avx_gt_detail::fp12ring().fp12_flat_frobenius_map2(
                frobenius2, avx_gt_detail::ring());
        const auto product =
            avx_gt_detail::fp12ring().fp12_flat_mul(
                frobenius4, p.data, avx_gt_detail::ring());
        return avx_gt_detail::fp12_equal(
            avx_gt_detail::narrow_to_standard(frobenius2), product);
    }

    /**
     * Serialize a GT element as twelve canonical big-endian 48-byte field
     * coefficients. Wire order is Fp2 index, then Fp6 half, then Fp
     * component, matching scalar `Fp12::to_bendian`; it is independent of the
     * native MatrixNoK layout. `serialized` points to 576 writable bytes.
     */
    static void serialize(unsigned char serialized[12 * 48], const GT& p)
    {
        BLST_UNROLL
        for (size_t wire = 0;
             wire < avx_gt_detail::kGtWireToNativeCoefficient.size();
             ++wire) {
            avx_gt_detail::fp_to_bendian(
                serialized + wire * 48,
                p.data[avx_gt_detail::kGtWireToNativeCoefficient[wire]]);
        }
    }

    /**
     * Decode twelve canonical big-endian field coefficients in the scalar
     * CatBLST order and require GT membership. The identity is accepted.
     * `serialized` points to 576 readable bytes. Failure leaves `out`
     * unchanged; input bytes may overlap `out`. Membership validation is
     * variable-time with respect to this public input.
     */
    [[nodiscard]] static bool deserialize(
        GT& out, const unsigned char serialized[12 * 48])
    {
        GT decoded = {};
        if (!deserialize_unchecked(decoded, serialized) || !in_group(decoded))
            return false;
        out = decoded;
        return true;
    }

    /**
     * Decode twelve canonical big-endian field coefficients without a GT
     * membership check. `serialized` points to 576 readable bytes. On success
     * `out` is only a general Fp12 value; combine candidates through the Fp12
     * API and establish membership before any ordinary GT operation. Failure
     * leaves `out` unchanged; input bytes may overlap `out`.
     */
    [[nodiscard]] static bool deserialize_unchecked(
        GT& out, const unsigned char serialized[12 * 48])
    {
        Field_type decoded = {};
        BLST_UNROLL
        for (size_t wire = 0;
             wire < avx_gt_detail::kGtWireToNativeCoefficient.size();
             ++wire) {
            if (!avx_gt_detail::fp_from_bendian(
                    decoded[avx_gt_detail::kGtWireToNativeCoefficient[wire]],
                    serialized + wire * 48))
                return false;
        }
        out.data = decoded;
        return true;
    }

    /**
     * Map a valid nonzero Miller result into GT using variable-time inversion.
     * `out` and `f` must not overlap. No allocation is performed.
     */
    static void final_exp(GT& out, const MillerResult& f)
    {
        matrixnok::PairingInverter inverter;
        out.data = blst::avx::final_exp(
            f.data, inverter, avx_gt_detail::fp12ring(), avx_gt_detail::ring());
    }

    /**
     * Batch-map `n` valid nonzero Miller results into GT. The arrays contain at
     * least `n` elements and must not overlap. Zero is a successful no-op.
     * `scratch` must be exclusively owned by this call. Its six buffers grow
     * to `n` and retain capacity for reuse; growth may throw `std::bad_alloc`.
     * No output is published before scratch growth and inversion preparation
     * complete. After sufficient capacity is established, this operation
     * performs no allocation.
     */
    static void final_exp_n(
        GT out[],
        const MillerResult in[],
        size_t n,
        FinalExpScratch& scratch)
    {
        if (n == 0 || out == nullptr || in == nullptr) {
            return;
        }

        const auto& fp12 = avx_gt_detail::fp12ring();
        const auto& ring = avx_gt_detail::ring();
        ScalarFpVartimeInverter<avx_gt_detail::Ring> inverter;
        static_assert(std::is_same_v<
            decltype(fp12.fp12_flat_inverse_to_fp6(std::declval<Field_type>(), ring)),
            FinalExpScratch::Fp6Element>);
        static_assert(std::is_same_v<
            decltype(fp12.fp6_flat_inverse_to_fp2_c0_c1_c2(
                std::declval<FinalExpScratch::Fp6Element>(), ring)),
            FinalExpScratch::Fp6Element>);
        static_assert(std::is_same_v<
            decltype(fp12.fp6_flat_inverse_to_fp2_ret(
                std::declval<FinalExpScratch::Fp6Element>(),
                std::declval<FinalExpScratch::Fp6Element>(),
                ring)),
            FinalExpScratch::Fp2Element>);

        scratch.ensure_size(n);
        auto& f6_terms = scratch.f6_terms_;
        auto& f2_c0_c1_c2_terms = scratch.f2_c0_c1_c2_terms_;
        auto& f2_terms = scratch.f2_terms_;
        auto& f1_terms = scratch.f1_terms_;
        auto& prefix_products = scratch.prefix_products_;
        auto& f1_inverses = scratch.f1_inverses_;
        for (size_t i = 0; i < n; ++i) {
            f6_terms[i] = fp12.fp12_flat_inverse_to_fp6(in[i].data, ring);
            f2_c0_c1_c2_terms[i] = fp12.fp6_flat_inverse_to_fp2_c0_c1_c2(f6_terms[i], ring);
            f2_terms[i] = fp12.fp6_flat_inverse_to_fp2_ret(f2_c0_c1_c2_terms[i], f6_terms[i], ring);
            f1_terms[i] = fp12.fp2_flat_inverse_to_fp(f2_terms[i], ring)[0];
        }

        // Batch the base-field inversion that dominates fp12 inversion in final exponentiation.
        prefix_products[0] = f1_terms[0];
        for (size_t i = 1; i < n; ++i) {
            prefix_products[i] = ring.modmul(prefix_products[i - 1], f1_terms[i]);
        }

        FinalExpScratch::FieldElement running_inverse =
            inverter.invert(prefix_products[n - 1], ring);
        for (size_t i = n; i-- > 1;) {
            f1_inverses[i] = ring.modmul(prefix_products[i - 1], running_inverse);
            running_inverse = ring.modmul(running_inverse, f1_terms[i]);
        }
        f1_inverses[0] = running_inverse;

        for (size_t i = 0; i < n; ++i) {
            std::array<FinalExpScratch::FieldElement, 1> f1_inv = {f1_inverses[i]};
            const auto f2_inv = fp12.fp_flat_inverse_to_fp2(f2_terms[i], f1_inv, ring);
            const auto f6_inv = fp12.fp2_flat_inverse_to_fp6(f2_c0_c1_c2_terms[i], f2_inv, ring);
            const auto f_inverse = fp12.fp6_flat_inverse_to_fp12(in[i].data, f6_inv, ring);
            out[i].data = final_exp_after_inverse(in[i].data, f_inverse, fp12, ring);
        }
    }

    /**
     * Pairing construction APIs. Non-identity inputs must already be validated
     * prime-order subgroup points; infinity is accepted as the subgroup
     * identity. No Release validation is performed. Outputs are MillerResult,
     * never GT, until passed through `final_exp`. The pointer-array overload
     * uses the compact continuation convention documented by the C API. These
     * operations are allocation-free and variable-time with respect to affine
     * infinity observations and public batch shape.
     */
    /** A null Q or P selects that group's canonical generator. */
    static void miller_loop(MillerResult& ret, const G2AffinePoint* Q, const G1AffinePoint* P)
    {
        const G2AffinePoint& q = Q ? *Q : G2AffinePoint::generator();
        const G1AffinePoint& p = P ? *P : G1AffinePoint::generator();

        if (G2AffinePoint::is_inf(q) || G1AffinePoint::is_inf(p)) {
            ret = MillerResult::identity();
            return;
        }

        const avx_gt_detail::PairingFp2 fp2ring(avx_gt_detail::ring());
        const matrixnok::MillerEngine miller(
            {p.X, p.Y},
            {q.X, q.Y},
            avx_gt_detail::ring());
        ret.data = miller.miller_loop(
            avx_gt_detail::ring(),
            fp2ring,
            avx_gt_detail::fp12ring());
    }

    /**
     * Multiply `n` pairings from contiguous arrays. For nonzero `n`, both
     * arrays contain at least `n` elements. Infinity pairs contribute identity.
     */
    static void miller_loop_n(
        MillerResult& out,
        const G2AffinePoint Q[],
        const G1AffinePoint P[],
        size_t n)
    {
        out = MillerResult::identity();
        if (n == 0) {
            return;
        }

        size_t offset = 0;
        while (offset < n) {
            G2Point T[avx_gt_detail::kMillerLoopNMax] = {};
            G2AffinePoint q_batch[avx_gt_detail::kMillerLoopNMax] = {};
            G1AffinePoint p_batch[avx_gt_detail::kMillerLoopNMax] = {};
            size_t batch = 0;

            while (offset < n && batch < avx_gt_detail::kMillerLoopNMax) {
                if (!G2AffinePoint::is_inf(Q[offset]) && !G1AffinePoint::is_inf(P[offset])) {
                    q_batch[batch] = Q[offset];
                    p_batch[batch] = P[offset];
                    T[batch].X = Q[offset].X;
                    T[batch].Y = Q[offset].Y;
                    T[batch].Z = avx_gt_detail::fp2ring().one();
                    ++batch;
                }
                ++offset;
            }

            if (batch == 0) {
                continue;
            }

            typename GT::Field_type chunk = {};
            avx_gt_detail::start_dbl_n(chunk, T, p_batch, batch);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 2);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 3);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 9);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 32);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 16);
            chunk = avx_gt_detail::fp12ring().fp12_flat_conjugate(chunk, avx_gt_detail::ring());
            out.data = avx_gt_detail::fp12ring().fp12_flat_mul(
                out.data,
                chunk,
                avx_gt_detail::ring());
        }
    }

    /**
     * Multiply `n` pairings from compact pointer arrays. The first pointer in
     * each array is nonnull; a later null pointer advances one element from the
     * previous resolved pointer. Infinity pairs contribute identity.
     */
    static void miller_loop_n(
        MillerResult& out,
        const G2AffinePoint* const Qs[],
        const G1AffinePoint* const Ps[],
        size_t n)
    {
        out = MillerResult::identity();
        if (n == 0) {
            return;
        }
        if (Qs == nullptr || Ps == nullptr || Qs[0] == nullptr || Ps[0] == nullptr) {
            out = {};
            return;
        }

        const G2AffinePoint* qptr = nullptr;
        const G1AffinePoint* pptr = nullptr;
        G2Point T[avx_gt_detail::kMillerLoopNMax] = {};
        G2AffinePoint q_batch[avx_gt_detail::kMillerLoopNMax] = {};
        G1AffinePoint p_batch[avx_gt_detail::kMillerLoopNMax] = {};
        size_t batch = 0;

        auto flush_batch = [&]() {
            if (batch == 0) {
                return;
            }
            typename GT::Field_type chunk = {};
            avx_gt_detail::start_dbl_n(chunk, T, p_batch, batch);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 2);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 3);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 9);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 32);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 16);
            chunk = avx_gt_detail::fp12ring().fp12_flat_conjugate(chunk, avx_gt_detail::ring());
            out.data = avx_gt_detail::fp12ring().fp12_flat_mul(
                out.data,
                chunk,
                avx_gt_detail::ring());
            batch = 0;
        };

        for (size_t i = 0; i < n; ++i) {
            qptr = Qs[i] ? Qs[i] : qptr + 1;
            pptr = Ps[i] ? Ps[i] : pptr + 1;
            BLST_ASSERT(qptr != nullptr);
            BLST_ASSERT(pptr != nullptr);

            if (G2AffinePoint::is_inf(*qptr) || G1AffinePoint::is_inf(*pptr)) {
                continue;
            }

            q_batch[batch] = *qptr;
            p_batch[batch] = *pptr;
            T[batch].X = qptr->X;
            T[batch].Y = qptr->Y;
            T[batch].Z = avx_gt_detail::fp2ring().one();
            ++batch;

            if (batch == avx_gt_detail::kMillerLoopNMax) {
                flush_batch();
            }
        }

        flush_batch();
    }

    /**
     * Fill all 68 native prepared lines for a subgroup G2 point. G2 infinity
     * produces the encoded identity table. Every output entry is overwritten.
     */
    static void precompute_lines(PreparedLine Qlines[68], const G2AffinePoint& Q)
    {
        if (G2AffinePoint::is_inf(Q)) {
            const PreparedLine identity_line = avx_gt_detail::identity_prepared_line();
            for (size_t i = 0; i < 68; ++i) {
                Qlines[i] = identity_line;
            }
            return;
        }

        G2Point t = {};
        t.X = Q.X;
        t.Y = Q.Y;
        t.Z = avx_gt_detail::fp2ring().one();

        avx_gt_detail::line_dbl(Qlines[0], t, t);
        avx_gt_detail::pre_add_n_dbl(&Qlines[1], t, Q, 2);
        avx_gt_detail::pre_add_n_dbl(&Qlines[4], t, Q, 3);
        avx_gt_detail::pre_add_n_dbl(&Qlines[8], t, Q, 9);
        avx_gt_detail::pre_add_n_dbl(&Qlines[18], t, Q, 32);
        avx_gt_detail::pre_add_n_dbl(&Qlines[51], t, Q, 16);
    }

    /**
     * Fill the scalar compatibility table for a subgroup G2 point. G2
     * infinity produces the encoded identity table.
     */
    static void precompute_lines(blst::scalar::Fp6 Qlines[68], const G2AffinePoint& Q)
    {
        blst::scalar::G2AffinePoint scalar_q = {};
        scalar_q.x = avx_curve_detail::ring_fp2_to_scalar(Q.X);
        scalar_q.y = avx_curve_detail::ring_fp2_to_scalar(Q.Y);
        blst::scalar::GT::precompute_lines(Qlines, scalar_q);
    }

    /**
     * Evaluate a complete native line table at a subgroup G1 point. G1
     * infinity returns the Miller identity. `Qlines` and `ret` do not overlap.
     */
    static void miller_loop_lines(
        MillerResult& ret,
        const PreparedLine Qlines[68],
        const G1AffinePoint& P)
    {
        if (G1AffinePoint::is_inf(P)) {
            ret = MillerResult::identity();
            return;
        }

        miller_loop_lines_finite(ret, Qlines, P);
    }

    /** Internal finite-point evaluator. The caller must already have excluded
     * G1 infinity; Release performs no observation or validation. */
    static void miller_loop_lines_finite(
        MillerResult& ret,
        const PreparedLine Qlines[68],
        const G1AffinePoint& P)
    {
        BLST_ASSERT(!G1AffinePoint::is_inf(P));
        const auto first_line = avx_gt_detail::scale_line_xy(Qlines[0], P);
        ret.data = avx_gt_detail::make_line_value(Qlines[0].c, first_line[0], first_line[1]);
        avx_gt_detail::apply_add_n_dbl(ret.data, &Qlines[1], P, 2);
        avx_gt_detail::apply_add_n_dbl(ret.data, &Qlines[4], P, 3);
        avx_gt_detail::apply_add_n_dbl(ret.data, &Qlines[8], P, 9);
        avx_gt_detail::apply_add_n_dbl(ret.data, &Qlines[18], P, 32);
        avx_gt_detail::apply_add_n_dbl(ret.data, &Qlines[51], P, 16);
        ret.data = avx_gt_detail::fp12ring().fp12_flat_conjugate(ret.data, avx_gt_detail::ring());
    }

    /**
     * Multiply prepared pairings from compact line-table and G1 pointer arrays.
     * The first pointer in each array is nonnull; a later null line pointer
     * advances by 68 lines and a later null point pointer advances by one point.
     * G1 infinity contributes identity. This operation is allocation-free.
     */
    static void miller_loop_lines_n(
        MillerResult& out,
        const PreparedLine* const Qlines[],
        const G1AffinePoint* const Ps[],
        size_t n)
    {
        out = MillerResult::identity();
        if (n == 0 || Qlines == nullptr || Ps == nullptr || Qlines[0] == nullptr || Ps[0] == nullptr) {
            return;
        }

        const PreparedLine* line_batch[avx_gt_detail::kMillerLoopNMax] = {};
        G1AffinePoint p_batch[avx_gt_detail::kMillerLoopNMax] = {};
        const PreparedLine* lines_ptr = nullptr;
        const G1AffinePoint* pptr = nullptr;
        size_t batch = 0;

        const auto flush_batch = [&]() {
            if (batch == 0) {
                return;
            }
            MillerResult chunk = {};
            miller_loop_lines_n_preprocessed(chunk, line_batch, p_batch, batch);
            MillerResult::mul(out, out, chunk);
            batch = 0;
        };

        for (size_t i = 0; i < n; ++i) {
            lines_ptr = Qlines[i] ? Qlines[i] : lines_ptr + 68;
            pptr = Ps[i] ? Ps[i] : pptr + 1;
            BLST_ASSERT(lines_ptr != nullptr);
            BLST_ASSERT(pptr != nullptr);

            if (G1AffinePoint::is_inf(*pptr)) {
                continue;
            }

            line_batch[batch] = lines_ptr;
            p_batch[batch] = *pptr;
            ++batch;
            if (batch == avx_gt_detail::kMillerLoopNMax) {
                flush_batch();
            }
        }

        flush_batch();
    }

    /**
     * Internal branch-free batch evaluator. Every line pointer is nonnull and
     * references 68 initialized entries; `P` contains `n` finite subgroup G1
     * points. Callers own filtering and batching. No allocation is performed.
     */
    static void miller_loop_lines_n_preprocessed(
        MillerResult& out,
        const PreparedLine* const Qlines[],
        const G1AffinePoint P[],
        size_t n)
    {
        out = MillerResult::identity();
        if (n == 0 || Qlines == nullptr || P == nullptr || Qlines[0] == nullptr) {
            return;
        }

        for (size_t i = 0; i < n; ++i) {
            BLST_ASSERT(Qlines[i] != nullptr);
            BLST_ASSERT(!G1AffinePoint::is_inf(P[i]));
        }

        avx_gt_detail::start_lines_n(out.data, Qlines, P, n);
        avx_gt_detail::post_add_n_dbl_n(out.data, Qlines, P, n, 1, 2);
        avx_gt_detail::post_add_n_dbl_n(out.data, Qlines, P, n, 4, 3);
        avx_gt_detail::post_add_n_dbl_n(out.data, Qlines, P, n, 8, 9);
        avx_gt_detail::post_add_n_dbl_n(out.data, Qlines, P, n, 18, 32);
        avx_gt_detail::post_add_n_dbl_n(out.data, Qlines, P, n, 51, 16);
        out.data = avx_gt_detail::fp12ring().fp12_flat_conjugate(out.data, avx_gt_detail::ring());
    }

    /** Evaluate a scalar compatibility line table at a subgroup G1 point. */
    static void miller_loop_lines(
        MillerResult& ret,
        const blst::scalar::Fp6 Qlines[68],
        const G1AffinePoint& P)
    {
        blst::scalar::G1AffinePoint scalar_p = {};
        scalar_p.x = avx_curve_detail::ring_fp_to_scalar(P.X);
        scalar_p.y = avx_curve_detail::ring_fp_to_scalar(P.Y);

        blst::scalar::MillerResult scalar_ret = {};
        blst::scalar::GT::miller_loop_lines(scalar_ret, Qlines, scalar_p);
        ret.data = avx_gt_detail::from_scalar(scalar_ret.data);
    }

    /** Multiply a GT element by the fixed signed BLS parameter z. */
    static void mult_by_z(GT& out, const GT& a)
    {
        out.data = avx_gt_detail::narrow_to_standard(
            raise_to_z(a.data, avx_gt_detail::fp12ring(), avx_gt_detail::ring()));
    }

    /** Multiply a GT element by the exact fixed exponent z/2. */
    static void mult_by_z_div_by_2(GT& out, const GT& a)
    {
        out.data = avx_gt_detail::narrow_to_standard(
            raise_to_z_div_by_2(a.data, avx_gt_detail::fp12ring(), avx_gt_detail::ring()));
    }

    /**
     * Apply Frobenius power `n`, where `n` is one, two, or three. The input
     * must lie in GT; Release does not validate the domain or exponent. Exact
     * output aliasing is supported.
     */
    static void frobenius_map(GT& out, const GT& a, size_t n)
    {
        BLST_ASSERT(n >= 1 && n <= 3);
        switch (n) {
        case 1:
            out.data = avx_gt_detail::fp12ring().fp12_flat_frobenius_map1(a.data, avx_gt_detail::ring());
            break;
        case 2:
            out.data = avx_gt_detail::fp12ring().fp12_flat_frobenius_map2(a.data, avx_gt_detail::ring());
            break;
        case 3:
            out.data = avx_gt_detail::fp12ring().fp12_flat_frobenius_map3(a.data, avx_gt_detail::ring());
            break;
        default:
            BLST_ASSERT(false && "Unsupported Frobenius power");
        }
    }

    /** Negate a GT element in place through Fp12 conjugation. */
    static void conjugate(GT& a)
    {
        a.data = avx_gt_detail::fp12ring().fp12_flat_conjugate(a.data, avx_gt_detail::ring());
    }
};

inline constexpr GT kGtGenerator(avx_gt_detail::gt_generator_value());
inline constexpr GT kGtIdentity(avx_gt_detail::fp12_identity_value());

inline const GT& GT::generator()
{
    return kGtGenerator;
}

}  // namespace blst::avx

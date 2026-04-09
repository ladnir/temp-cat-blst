#pragma once

#include <array>
#include <cstddef>
#include <cstring>
#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>

#include "blst/avx/curve/point.h"
#include "blst/avx/curve/scalar_mult.h"
#include "blst/common/defines.h"
#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/primitveRoot.h"
#include "blst/fields/Fr.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst::avx {

template<class Group>
struct DenseTransformCoeffCtx {
    static auto plus(const Group& a, const Group& b)
    {
        Group r = {};
        Group::dadd(r, a, b);
        return r;
    }

    static auto minus(const Group& a, const Group& b)
    {
        Group neg_b = {};
        Group::neg(neg_b, b);
        Group r = {};
        Group::dadd(r, a, neg_b);
        return r;
    }

    static auto mult(const Group& a, const Fr& b)
    {
        Group r = {};
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

template<std::size_t M>
struct G1DenseTransformPlan {
    static_assert(M > 0);
    static_assert((M & (M - 1)) == 0);
    static_assert(M <= 32);

    std::array<Fr, M> roots = {};
    std::array<std::array<uint8_t, 32>, M> root_scalar_bytes = {};
    std::array<uint8_t, M * M> twiddle_indexes = {};
};

namespace detail {

template<std::size_t M>
inline std::array<uint8_t, 32> fr_to_scalar_bytes(const Fr& value)
{
    Fr value_non_mont = {};
    Fr::from_mont(value_non_mont, value);

    std::array<uint8_t, 32> out = {};
    static_assert(sizeof(value_non_mont) == out.size());
    std::memcpy(out.data(), &value_non_mont, out.size());
    return out;
}

template<std::size_t M>
inline uint8_t find_root_index(
    const std::array<Fr, M>& roots,
    const Fr& value)
{
    for (std::size_t i = 0; i < M; ++i) {
        if (Fr::equal(roots[i], value)) {
            return static_cast<uint8_t>(i);
        }
    }

    BLST_ASSERT(false && "dense transform coefficient was not a root of unity");
    return 0;
}

template<std::size_t M>
inline auto make_precomputed_tables(std::span<const G1Point> in)
{
    using Curve = avx_curve_detail::CurveG1;
    using Ring = avx_curve_detail::CurveRing;
    using Table = PrecomputedTable<Curve, Ring, 4>;
    using ProjectivePoint = typename Curve::ProjPoint;

    const Curve& curve = avx_curve_detail::curve_g1();
    const Ring& ring = avx_curve_detail::curve_ring();

    std::vector<Table> out;
    out.reserve(M);
    for (std::size_t i = 0; i < M; ++i) {
        out.emplace_back(ProjectivePoint(in[i].X, in[i].Y, in[i].Z), curve, ring);
    }
    return out;
}

template<std::size_t M, std::size_t... I>
inline auto make_row_scalars_impl(
    const G1DenseTransformPlan<M>& plan,
    std::size_t row,
    std::index_sequence<I...>)
{
    return std::array<ScalarForExponent<256>, M>{
        ScalarForExponent<256>(
            plan.root_scalar_bytes[plan.twiddle_indexes[(row * M) + I]])...
    };
}

template<std::size_t M>
inline auto make_row_scalars(
    const G1DenseTransformPlan<M>& plan,
    std::size_t row)
{
    return make_row_scalars_impl<M>(plan, row, std::make_index_sequence<M>{});
}

template<std::size_t M>
inline auto dense_table_msm(
    const std::array<ScalarForExponent<256>, M>& scalars,
    std::span<const PrecomputedTable<avx_curve_detail::CurveG1, avx_curve_detail::CurveRing, 4>> tables)
{
    using Curve = avx_curve_detail::CurveG1;
    using Ring = avx_curve_detail::CurveRing;
    using ProjectivePoint = typename Curve::ProjPoint;

    constexpr int kTableLogSize = 4;
    constexpr int kScalarBits = 256;
    constexpr int kNumWindows = (kScalarBits + kTableLogSize) / (kTableLogSize + 1);
    constexpr bool kCarryPossible = (kScalarBits % (kTableLogSize + 1)) == 0;

    const Curve& curve = avx_curve_detail::curve_g1();
    const Ring& ring = avx_curve_detail::curve_ring();

    std::array<std::array<std::pair<uint8_t, bool>, kNumWindows>, M> booth_indexes = {};
    ProjectivePoint result = curve.zero(ring);

    for (std::size_t i = 0; i < M; ++i) {
        bool carry = false;
        for (int window = 0; window < kNumWindows; ++window) {
            uint64_t window_val =
                scalars[i].template get_bits<kTableLogSize + 1>(window * (kTableLogSize + 1));
            window_val += carry ? 1u : 0u;
            booth_indexes[i][window] = tables[i].booth_encode(window_val);
            carry = booth_indexes[i][window].second;
        }

        if constexpr (kCarryPossible) {
            ProjectivePoint carry_point = curve.zero(ring);
            carry_point = carry_point.ctime_replace_if(
                carry,
                tables[i].data[(1 << kTableLogSize) - 1]);
            result = curve.add_point(result, carry_point, ring);
        }
    }

    if constexpr (kCarryPossible) {
        result = curve.double_point(result, ring);
    }

    for (int window = kNumWindows - 1; window >= 0; --window) {
        for (std::size_t i = 0; i < M; ++i) {
            result = curve.add_point(
                result,
                tables[i].gather_booth(
                    booth_indexes[i][window].first,
                    booth_indexes[i][window].second,
                    curve,
                    ring),
                ring);
        }

        if (window != 0) {
            for (int bit = 0; bit < (kTableLogSize + 1); ++bit) {
                result = curve.double_point(result, ring);
            }
        }
    }

    return result;
}

} // namespace detail

template<std::size_t M>
inline G1DenseTransformPlan<M> make_g1_dense_transform_plan()
{
    G1DenseTransformPlan<M> plan = {};

    const Fr primitive = primRootOfUnity<Fr>(M);
    Fr current_root = Fr::one();
    for (std::size_t i = 0; i < M; ++i) {
        plan.roots[i] = current_root;
        plan.root_scalar_bytes[i] = detail::fr_to_scalar_bytes<M>(current_root);
        Fr::mul(current_root, current_root, primitive);
    }

    for (std::size_t column = 0; column < M; ++column) {
        std::vector<Fr> basis(M);
        basis[column] = Fr::one();
        fft<Fr, Fr, FrCoeffCtx>(
            std::span<Fr>(basis.data(), basis.size()),
            std::span<const Fr>(plan.roots.data(), plan.roots.size()),
            false);

        for (std::size_t row = 0; row < M; ++row) {
            plan.twiddle_indexes[(row * M) + column] =
                detail::find_root_index(plan.roots, basis[row]);
        }
    }

    return plan;
}

template<std::size_t M>
inline void g1_dense_transform(
    std::span<const G1Point> in,
    const G1DenseTransformPlan<M>& plan,
    std::span<G1Point> out)
{
    BLST_ASSERT(in.size() == M);
    BLST_ASSERT(out.size() >= M);

    const auto tables = detail::make_precomputed_tables<M>(in);
    for (std::size_t row = 0; row < M; ++row) {
        const auto scalars = detail::make_row_scalars(plan, row);
        const auto raw = detail::dense_table_msm<M>(
            scalars,
            std::span<const PrecomputedTable<avx_curve_detail::CurveG1, avx_curve_detail::CurveRing, 4>>(
                tables.data(),
                tables.size()));
        out[row] = G1Point(raw.x, raw.y, raw.z);
    }
}

} // namespace blst::avx

#endif

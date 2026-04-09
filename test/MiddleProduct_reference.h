#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>

#include "blst/common/timer.h"
#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fft/middle_product.h"
#include "blst/scalar/curve/GT.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/curve/GT.h"
#include "blst/avx/fft/gt_context.h"
#endif

namespace blst::test_support {

template<class Context>
using DirectPairingGt = std::conditional_t<
    std::is_same_v<typename Context::G1Point, blst::scalar::G1Point>,
    blst::scalar::GT,
    blst::avx::GT>;

template<class Context>
inline GT direct_pairing_gt_to_public(const DirectPairingGt<Context>& value)
{
    if constexpr (std::is_same_v<DirectPairingGt<Context>, blst::scalar::GT>) {
        return Context::to_public_gt(value);
    } else {
        return blst::avx::AvxGt52::to_default(blst::avx::AvxGt52(value.data));
    }
}

template<class Context>
struct DirectMiddleProductPrecomputedB {
    std::size_t n = 0;
    MiddleProductPairingMode pairing_mode = MiddleProductPairingMode::Direct;
    std::vector<typename Context::G2Point> b;
    std::vector<typename Context::G2AffinePoint> b_affine;
    std::vector<typename Context::PreparedLineTable> prepared_lines;
    std::vector<uint8_t> b_is_zero;
};

template<class Context>
inline std::vector<typename Context::G1Point> direct_convert_g1_inputs(
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
inline std::vector<typename Context::G2Point> direct_convert_g2_inputs(
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
inline std::vector<uint8_t> direct_make_g2_affine(
    std::vector<typename Context::G2AffinePoint>& out,
    std::span<const typename Context::G2Point> input)
{
    std::vector<uint8_t> is_zero(input.size(), 0);
    out.resize(input.size());
    for (std::size_t j = 0; j < input.size(); ++j) {
        const bool zero = (Context::G2Point::is_inf(input[j]) != 0);
        is_zero[j] = static_cast<uint8_t>(zero ? 1 : 0);
        if (!zero) {
            Context::G2Point::to_affine(out[j], input[j]);
        }
    }
    return is_zero;
}

template<class Context>
inline void direct_make_g1_affine(
    std::vector<typename Context::G1AffinePoint>& out_affine,
    std::vector<uint8_t>& out_is_zero,
    std::span<const typename Context::G1Point> input)
{
    out_affine.resize(input.size());
    out_is_zero.resize(input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
        const bool zero = (Context::G1Point::is_inf(input[i]) != 0);
        out_is_zero[i] = static_cast<uint8_t>(zero ? 1 : 0);
        if (!zero) {
            Context::G1Point::to_affine(out_affine[i], input[i]);
        }
    }
}

template<class Context>
inline DirectMiddleProductPrecomputedB<Context> precompute_middle_product_b_direct(
    std::span<const G2Point> b,
    std::size_t n,
    TimingStats* timings = nullptr,
    MiddleProductPairingMode pairing_mode = MiddleProductPairingMode::Direct)
{
    ScopedTimer total_timer(timings, "middle_product_direct_precompute.total_ns");
    BLST_ASSERT(n > 0);
    BLST_ASSERT(b.size() == 2 * n);

    DirectMiddleProductPrecomputedB<Context> precomputed;
    precomputed.n = n;
    precomputed.pairing_mode = pairing_mode;
    precomputed.b = direct_convert_g2_inputs<Context>(b);

    {
        ScopedTimer timer(timings, "middle_product_direct_precompute.zero_mask_ns");
        precomputed.b_is_zero =
            direct_make_g2_affine<Context>(precomputed.b_affine, precomputed.b);
    }

    if (pairing_mode == MiddleProductPairingMode::PrecomputedLines) {
        ScopedTimer timer(timings, "middle_product_direct_precompute.lines_ns");
        precomputed.prepared_lines.resize(precomputed.b.size());
        for (std::size_t j = 0; j < precomputed.b.size(); ++j) {
            if (precomputed.b_is_zero[j] != 0) {
                continue;
            }
            Context::prepare_lines(precomputed.prepared_lines[j], precomputed.b[j]);
        }
    }

    return precomputed;
}

template<class Context>
inline GT direct_pair_term(
    const typename Context::G1Point& left,
    const DirectMiddleProductPrecomputedB<Context>& precomputed_b,
    std::size_t b_index)
{
    if (precomputed_b.b_is_zero[b_index] != 0) {
        return GT::identity();
    }

    if (precomputed_b.pairing_mode == MiddleProductPairingMode::PrecomputedLines &&
        precomputed_b.prepared_lines.size() == precomputed_b.b.size()) {
        return Context::to_public_gt(
            Context::pointwise_pair_lines(precomputed_b.prepared_lines[b_index], left));
    }

    return Context::to_public_gt(
        Context::pointwise_pair_direct(left, precomputed_b.b[b_index]));
}

template<class Context>
inline GT direct_multipair_window(
    std::span<const typename Context::G1AffinePoint> a_affine_nonzero,
    std::span<const std::size_t> a_indices_nonzero,
    const DirectMiddleProductPrecomputedB<Context>& precomputed_b,
    std::size_t output_index,
    std::vector<typename Context::G1AffinePoint>& p_batch,
    std::vector<typename Context::G2AffinePoint>& q_batch)
{
    p_batch.clear();
    q_batch.clear();

    for (std::size_t nz = 0; nz < a_affine_nonzero.size(); ++nz) {
        const std::size_t i = a_indices_nonzero[nz];
        const std::size_t j = precomputed_b.n + output_index - i;
        if (precomputed_b.b_is_zero[j] != 0) {
            continue;
        }
        p_batch.push_back(a_affine_nonzero[nz]);
        q_batch.push_back(precomputed_b.b_affine[j]);
    }

    if (p_batch.empty()) {
        return GT::identity();
    }

    using PairingGt = DirectPairingGt<Context>;

    PairingGt miller = PairingGt::identity();
    PairingGt::miller_loop_n(
        miller,
        q_batch.data(),
        p_batch.data(),
        q_batch.size());

    PairingGt final_pair = PairingGt::identity();
    PairingGt::final_mult(final_pair, miller);
    return direct_pairing_gt_to_public<Context>(final_pair);
}

template<class Context>
inline GT direct_multipair_lines_window(
    std::span<const typename Context::G1AffinePoint> a_affine_nonzero,
    std::span<const std::size_t> a_indices_nonzero,
    const DirectMiddleProductPrecomputedB<Context>& precomputed_b,
    std::size_t output_index)
{
    using PairingGt = DirectPairingGt<Context>;

    PairingGt miller_acc = PairingGt::identity();
    bool has_term = false;

    for (std::size_t nz = 0; nz < a_affine_nonzero.size(); ++nz) {
        const std::size_t i = a_indices_nonzero[nz];
        const std::size_t j = precomputed_b.n + output_index - i;
        if (precomputed_b.b_is_zero[j] != 0) {
            continue;
        }

        PairingGt miller_term = PairingGt::identity();
        PairingGt::miller_loop_lines(
            miller_term,
            precomputed_b.prepared_lines[j].lines,
            a_affine_nonzero[nz]);
        PairingGt::add(miller_acc, miller_acc, miller_term);
        has_term = true;
    }

    if (!has_term) {
        return GT::identity();
    }

    PairingGt final_pair = PairingGt::identity();
    PairingGt::final_mult(final_pair, miller_acc);
    return direct_pairing_gt_to_public<Context>(final_pair);
}

template<class Context>
inline void middle_product_direct(
    std::span<const G1Point> a,
    const DirectMiddleProductPrecomputedB<Context>& precomputed_b,
    std::span<GT> out,
    TimingStats* timings = nullptr)
{
    ScopedTimer total_timer(timings, "middle_product_direct.total_ns");
    const std::size_t n = a.size();
    BLST_ASSERT(n > 0);
    BLST_ASSERT(precomputed_b.n == n);
    BLST_ASSERT(precomputed_b.b.size() == 2 * n);
    BLST_ASSERT(out.size() >= n);

    {
        ScopedTimer timer(timings, "middle_product_direct.fill_ns");
        std::fill(out.begin(), out.begin() + n, GT::identity());
    }

    const std::vector<typename Context::G1Point> converted_a = direct_convert_g1_inputs<Context>(a);
    std::vector<typename Context::G1AffinePoint> affine_a;
    std::vector<uint8_t> a_is_zero;
    direct_make_g1_affine<Context>(affine_a, a_is_zero, std::span<const typename Context::G1Point>(converted_a));
    std::vector<typename Context::G1AffinePoint> affine_a_nonzero;
    std::vector<std::size_t> a_indices_nonzero;
    affine_a_nonzero.reserve(n);
    a_indices_nonzero.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (a_is_zero[i] != 0) {
            continue;
        }
        affine_a_nonzero.push_back(affine_a[i]);
        a_indices_nonzero.push_back(i);
    }

    {
        ScopedTimer timer(timings, "middle_product_direct.pointwise_ns");
        if (precomputed_b.pairing_mode == MiddleProductPairingMode::Direct) {
            std::vector<typename Context::G1AffinePoint> p_batch;
            std::vector<typename Context::G2AffinePoint> q_batch;
            p_batch.reserve(affine_a_nonzero.size());
            q_batch.reserve(affine_a_nonzero.size());
            for (std::size_t k = 0; k < n; ++k) {
                out[k] = direct_multipair_window<Context>(
                    std::span<const typename Context::G1AffinePoint>(affine_a_nonzero),
                    std::span<const std::size_t>(a_indices_nonzero),
                    precomputed_b,
                    k,
                    p_batch,
                    q_batch);
            }
        } else {
            for (std::size_t k = 0; k < n; ++k) {
                out[k] = direct_multipair_lines_window<Context>(
                    std::span<const typename Context::G1AffinePoint>(affine_a_nonzero),
                    std::span<const std::size_t>(a_indices_nonzero),
                    precomputed_b,
                    k);
            }
        }
    }
}

template<class Context>
inline void update_direct(
    const G1Point& a0,
    const DirectMiddleProductPrecomputedB<Context>& precomputed_b,
    std::span<GT> out,
    TimingStats* timings = nullptr)
{
    ScopedTimer total_timer(timings, "middle_product_direct_update.total_ns");
    BLST_ASSERT(precomputed_b.n > 0);
    BLST_ASSERT(precomputed_b.b.size() == 2 * precomputed_b.n);
    BLST_ASSERT(out.size() >= precomputed_b.n);

    const typename Context::G1Point converted_a0 = Context::to_backend_g1(a0);
    if (Context::G1Point::is_inf(converted_a0) != 0) {
        std::fill(out.begin(), out.begin() + precomputed_b.n, GT::identity());
        return;
    }

    {
        ScopedTimer timer(timings, "middle_product_direct_update.pointwise_ns");
        for (std::size_t k = 0; k < precomputed_b.n; ++k) {
            out[k] = direct_pair_term<Context>(converted_a0, precomputed_b, precomputed_b.n + k);
        }
    }
}

template<class Context>
inline void combine_direct(
    std::span<const GT> mid,
    std::span<const GT> cnst,
    std::span<GT> out)
{
    combine(mid, cnst, out);
}

} // namespace blst::test_support

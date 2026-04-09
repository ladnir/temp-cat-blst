#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "blst/common/timer.h"
#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fft/middle_product.h"

namespace blst::test_support {

template<class Context>
struct DirectMiddleProductPrecomputedB {
    std::size_t n = 0;
    MiddleProductPairingMode pairing_mode = MiddleProductPairingMode::Direct;
    std::vector<typename Context::G2Point> b;
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
    precomputed.b_is_zero.resize(precomputed.b.size());

    {
        ScopedTimer timer(timings, "middle_product_direct_precompute.zero_mask_ns");
        for (std::size_t j = 0; j < precomputed.b.size(); ++j) {
            const bool is_zero = (Context::G2Point::is_inf(precomputed.b[j]) != 0);
            precomputed.b_is_zero[j] = static_cast<uint8_t>(is_zero ? 1 : 0);
        }
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

    {
        ScopedTimer timer(timings, "middle_product_direct.pointwise_ns");
        for (std::size_t i = 0; i < n; ++i) {
            if (Context::G1Point::is_inf(converted_a[i]) != 0) {
                continue;
            }

            for (std::size_t k = 0; k < n; ++k) {
                const std::size_t j = n + k - i;
                GT acc = out[k];
                const GT term = direct_pair_term<Context>(converted_a[i], precomputed_b, j);
                GT::add(acc, acc, term);
                out[k] = acc;
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

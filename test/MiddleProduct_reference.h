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
inline void direct_collect_nonzero_g1_affine(
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

inline void direct_prepare_scalar_line_eval_points(
    std::vector<blst::scalar::G1AffinePoint>& out_px2,
    std::span<const blst::scalar::G1AffinePoint> input)
{
    out_px2.clear();
    out_px2.reserve(input.size());
    for (const blst::scalar::G1AffinePoint& point : input) {
        blst::scalar::G1AffinePoint px2 = {};
        blst::scalar::Fp::add(px2.X, point.X, point.X);
        blst::scalar::Fp::neg(px2.X, px2.X);
        blst::scalar::Fp::add(px2.Y, point.Y, point.Y);
        out_px2.push_back(px2);
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
inline GT direct_pair_term_affine(
    const typename Context::G1AffinePoint& left_affine,
    const DirectMiddleProductPrecomputedB<Context>& precomputed_b,
    std::size_t b_index)
{
    if (precomputed_b.b_is_zero[b_index] != 0) {
        return GT::identity();
    }

    using PairingGt = DirectPairingGt<Context>;

    if (precomputed_b.pairing_mode == MiddleProductPairingMode::PrecomputedLines &&
        precomputed_b.prepared_lines.size() == precomputed_b.b.size()) {
        PairingGt miller = PairingGt::identity();
        PairingGt::miller_loop_lines(
            miller,
            precomputed_b.prepared_lines[b_index].lines,
            left_affine);
        PairingGt final_pair = PairingGt::identity();
        PairingGt::final_mult(final_pair, miller);
        return direct_pairing_gt_to_public<Context>(final_pair);
    }

    PairingGt miller = PairingGt::identity();
    PairingGt::miller_loop(
        miller,
        &precomputed_b.b_affine[b_index],
        &left_affine);
    PairingGt final_pair = PairingGt::identity();
    PairingGt::final_mult(final_pair, miller);
    return direct_pairing_gt_to_public<Context>(final_pair);
}

template<class Context>
inline GT direct_multipair_window(
    std::span<const typename Context::G1AffinePoint> a_affine_nonzero,
    std::span<const std::size_t> a_indices_nonzero,
    const DirectMiddleProductPrecomputedB<Context>& precomputed_b,
    std::size_t output_index,
    std::vector<const typename Context::G1AffinePoint*>& p_batch,
    std::vector<const typename Context::G2AffinePoint*>& q_batch)
{
    p_batch.clear();
    q_batch.clear();

    for (std::size_t nz = 0; nz < a_affine_nonzero.size(); ++nz) {
        const std::size_t i = a_indices_nonzero[nz];
        const std::size_t j = precomputed_b.n + output_index - i;
        if (precomputed_b.b_is_zero[j] != 0) {
            continue;
        }
        p_batch.push_back(&a_affine_nonzero[nz]);
        q_batch.push_back(&precomputed_b.b_affine[j]);
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
    std::size_t output_index,
    std::vector<const typename Context::G1AffinePoint*>& p_batch,
    std::vector<const blst::scalar::Fp6*>* line_batch_scalar = nullptr)
{
    using PairingGt = DirectPairingGt<Context>;

    p_batch.clear();

    if constexpr (std::is_same_v<PairingGt, blst::scalar::GT>) {
        BLST_ASSERT(line_batch_scalar != nullptr);
        auto& line_batch = *line_batch_scalar;
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

        if (p_batch.empty()) {
            return GT::identity();
        }

        PairingGt miller = PairingGt::identity();
        PairingGt::miller_loop_lines_n(
            miller,
            line_batch.data(),
            p_batch.data(),
            p_batch.size());

        PairingGt final_pair = PairingGt::identity();
        PairingGt::final_mult(final_pair, miller);
        return direct_pairing_gt_to_public<Context>(final_pair);
    } else {
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
}

inline GT direct_multipair_lines_window_scalar_preprocessed(
    std::span<const blst::scalar::G1AffinePoint> a_line_points,
    std::span<const std::size_t> a_indices_nonzero,
    const DirectMiddleProductPrecomputedB<blst::MiddleProductScalarContext>& precomputed_b,
    std::size_t output_index,
    std::vector<const blst::scalar::Fp6*>& line_batch,
    std::vector<blst::scalar::G1AffinePoint>& px2_batch)
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

    blst::scalar::GT miller = blst::scalar::GT::identity();
    blst::scalar::GT::miller_loop_lines_n_preprocessed(
        miller,
        line_batch.data(),
        px2_batch.data(),
        line_batch.size());

    blst::scalar::GT final_pair = blst::scalar::GT::identity();
    blst::scalar::GT::final_mult(final_pair, miller);
    return direct_pairing_gt_to_public<blst::MiddleProductScalarContext>(final_pair);
}

inline GT direct_multipair_lines_window_scalar_dense(
    std::span<const blst::scalar::G1AffinePoint> a_line_points,
    const DirectMiddleProductPrecomputedB<blst::MiddleProductScalarContext>& precomputed_b,
    std::size_t output_index,
    std::vector<const blst::scalar::Fp6*>& line_batch)
{
    const std::size_t n = precomputed_b.n;
    line_batch.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        line_batch[i] = precomputed_b.prepared_lines[n + output_index - i].lines;
    }

    blst::scalar::GT miller = blst::scalar::GT::identity();
    blst::scalar::GT::miller_loop_lines_n_preprocessed(
        miller,
        line_batch.data(),
        a_line_points.data(),
        line_batch.size());

    blst::scalar::GT final_pair = blst::scalar::GT::identity();
    blst::scalar::GT::final_mult(final_pair, miller);
    return direct_pairing_gt_to_public<blst::MiddleProductScalarContext>(final_pair);
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

    const std::vector<typename Context::G1Point> converted_a = direct_convert_g1_inputs<Context>(a);
    std::vector<typename Context::G1AffinePoint> affine_a_nonzero;
    std::vector<std::size_t> a_indices_nonzero;
    direct_collect_nonzero_g1_affine<Context>(
        affine_a_nonzero,
        a_indices_nonzero,
        std::span<const typename Context::G1Point>(converted_a));

    {
        ScopedTimer timer(timings, "middle_product_direct.pointwise_ns");
        if (precomputed_b.pairing_mode == MiddleProductPairingMode::Direct) {
            std::vector<const typename Context::G1AffinePoint*> p_batch;
            std::vector<const typename Context::G2AffinePoint*> q_batch;
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
            if constexpr (std::is_same_v<DirectPairingGt<Context>, blst::scalar::GT>) {
                std::vector<blst::scalar::G1AffinePoint> a_line_points;
                direct_prepare_scalar_line_eval_points(
                    a_line_points,
                    std::span<const blst::scalar::G1AffinePoint>(affine_a_nonzero));

                std::vector<const blst::scalar::Fp6*> line_batch_scalar;
                line_batch_scalar.reserve(affine_a_nonzero.size());
                const bool dense_a = (a_indices_nonzero.size() == n);
                const bool dense_b =
                    std::find(
                        precomputed_b.b_is_zero.begin(),
                        precomputed_b.b_is_zero.end(),
                        static_cast<uint8_t>(1)) == precomputed_b.b_is_zero.end();

                if (dense_a && dense_b) {
                    for (std::size_t k = 0; k < n; ++k) {
                        out[k] = direct_multipair_lines_window_scalar_dense(
                            std::span<const blst::scalar::G1AffinePoint>(a_line_points),
                            precomputed_b,
                            k,
                            line_batch_scalar);
                    }
                } else {
                    std::vector<blst::scalar::G1AffinePoint> px2_batch_scalar;
                    px2_batch_scalar.reserve(affine_a_nonzero.size());
                    for (std::size_t k = 0; k < n; ++k) {
                        out[k] = direct_multipair_lines_window_scalar_preprocessed(
                            std::span<const blst::scalar::G1AffinePoint>(a_line_points),
                            std::span<const std::size_t>(a_indices_nonzero),
                            precomputed_b,
                            k,
                            line_batch_scalar,
                            px2_batch_scalar);
                    }
                }
            } else {
                std::vector<const typename Context::G1AffinePoint*> p_batch;
                p_batch.reserve(affine_a_nonzero.size());
                std::vector<const blst::scalar::Fp6*> line_batch_scalar;
                for (std::size_t k = 0; k < n; ++k) {
                    out[k] = direct_multipair_lines_window<Context>(
                        std::span<const typename Context::G1AffinePoint>(affine_a_nonzero),
                        std::span<const std::size_t>(a_indices_nonzero),
                        precomputed_b,
                        k,
                        p_batch,
                        &line_batch_scalar);
                }
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

    typename Context::G1AffinePoint affine_a0 = {};
    Context::G1Point::to_affine(affine_a0, converted_a0);

    {
        ScopedTimer timer(timings, "middle_product_direct_update.pointwise_ns");
        for (std::size_t k = 0; k < precomputed_b.n; ++k) {
            out[k] = direct_pair_term_affine<Context>(affine_a0, precomputed_b, precomputed_b.n + k);
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

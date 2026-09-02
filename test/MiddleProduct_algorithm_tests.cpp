#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include "blst/fft/middle_product.h"

namespace {

using blst::GT;
using blst::G1Point;
using blst::G2Point;
using blst::MiddleProductAlgorithm;
using blst::MiddleProductScalarContext;

G1Point g1(uint64_t scalar)
{
    G1Point out = {};
    blst::Fr s = {};
    blst::Fr::from_uint64(s, scalar);
    blst::Fr s_non_mont = {};
    blst::Fr::from_mont(s_non_mont, s);
    G1Point::mult(
        out,
        G1Point::generator(),
        reinterpret_cast<const blst::byte*>(&s_non_mont),
        256);
    return out;
}

G2Point g2(uint64_t scalar)
{
    G2Point out = {};
    blst::Fr s = {};
    blst::Fr::from_uint64(s, scalar);
    blst::Fr s_non_mont = {};
    blst::Fr::from_mont(s_non_mont, s);
    G2Point::mult(
        out,
        G2Point::generator(),
        reinterpret_cast<const blst::byte*>(&s_non_mont),
        256);
    return out;
}

std::vector<G1Point> make_g1_data(std::size_t n)
{
    std::vector<G1Point> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        if ((i % 7) == 3) {
            out.push_back(G1Point::infinity());
            continue;
        }
        out.push_back(g1(static_cast<uint64_t>((37 * i + 11) % 257 + 1)));
    }
    return out;
}

std::vector<G2Point> make_g2_data(std::size_t n)
{
    std::vector<G2Point> out;
    out.reserve(2 * n);
    for (std::size_t i = 0; i < (2 * n); ++i) {
        if ((i % 11) == 5) {
            out.push_back(G2Point::infinity());
            continue;
        }
        out.push_back(g2(static_cast<uint64_t>((43 * i + 17) % 263 + 1)));
    }
    return out;
}

template<class Function>
void for_shuffled_ranges(
    const std::size_t count,
    std::uint64_t seed,
    Function&& function)
{
    std::vector<std::pair<std::size_t, std::size_t>> ranges;
    for (std::size_t begin = 0; begin < count;) {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        const std::size_t length = 1 + static_cast<std::size_t>(seed % 7);
        const std::size_t end = std::min(count, begin + length);
        ranges.emplace_back(begin, end);
        begin = end;
    }
    for (std::size_t i = ranges.size(); i > 1; --i) {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        std::swap(ranges[i - 1], ranges[seed % i]);
    }
    for (const auto [begin, end] : ranges) {
        function(begin, end);
    }
}

template<class Function>
bool for_parallel_ranges(const std::size_t count, Function&& function)
{
    if (count == 0) {
        return true;
    }
    constexpr std::size_t max_workers = 4;
    const std::size_t worker_count = std::min(count, max_workers);
    std::array<std::atomic<bool>, max_workers> succeeded;
    for (std::atomic<bool>& value : succeeded) {
        value.store(true, std::memory_order_relaxed);
    }
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (std::size_t worker = 0; worker < worker_count; ++worker) {
        const std::size_t begin = (count * worker) / worker_count;
        const std::size_t end = (count * (worker + 1)) / worker_count;
        workers.emplace_back([&, worker, begin, end]() {
            try {
                function(begin, end);
            } catch (...) {
                succeeded[worker].store(false, std::memory_order_relaxed);
            }
        });
    }
    for (std::thread& worker : workers) {
        worker.join();
    }
    for (std::size_t worker = 0; worker < worker_count; ++worker) {
        if (!succeeded[worker].load(std::memory_order_relaxed)) {
            return false;
        }
    }
    return true;
}

template<class Context>
void expect_partitioned_ranges_match_monolithic(const std::size_t n)
{
    const std::vector<G1Point> public_a = make_g1_data(n);
    const std::vector<G2Point> public_b = make_g2_data(n);
    const auto precomputed =
        blst::detail::precompute_middle_product_b_public<Context>(
            public_b, n, nullptr);
    const std::vector<typename Context::G1Point> native_a =
        blst::detail::convert_g1_inputs<Context>(public_a);
    blst::detail::MiddleProductWorkspaceImpl<Context> monolithic_workspace(
        precomputed.fft_size);
    blst::detail::MiddleProductWorkspaceImpl<Context> ranged_workspace(
        precomputed.fft_size);
    std::vector<GT> expected(n, GT::identity());
    std::vector<GT> actual(n, GT::identity());

    blst::detail::middle_product_impl<Context>(
        std::span<const typename Context::G1Point>(native_a),
        precomputed,
        std::span<GT>(expected),
        monolithic_workspace,
        nullptr);

    for_shuffled_ranges(
        precomputed.fft_size,
        0x123456789abcdef0ULL + n,
        [&](const std::size_t begin, const std::size_t end) {
            blst::detail::middle_product_input_range_impl<Context>(
                native_a, precomputed, ranged_workspace, begin, end);
        });

    const auto execute_partitioned = [&](
        const blst::MiddleProductRangePhase phase,
        const std::size_t stage,
        const std::uint64_t seed) {
        const std::size_t count =
            blst::detail::middle_product_range_work_count_impl<Context>(
                precomputed, phase, stage);
        for_shuffled_ranges(
            count,
            seed,
            [&](const std::size_t begin, const std::size_t end) {
                blst::detail::middle_product_execute_range_impl<Context>(
                    precomputed,
                    ranged_workspace,
                    phase,
                    stage,
                    begin,
                    end);
            });
    };

    execute_partitioned(
        blst::MiddleProductRangePhase::ForwardBitReverse,
        0,
        0x201 + n);
    for (std::size_t stage = 0;
         stage < precomputed.forward_plan.stage_work.size();
         ++stage) {
        const auto execute_j0 = [&]() {
            execute_partitioned(
                blst::MiddleProductRangePhase::ForwardJ0,
                stage,
                0x300 + stage + n);
        };
        const auto execute_nonzero = [&]() {
            execute_partitioned(
                blst::MiddleProductRangePhase::ForwardNonzero,
                stage,
                0x400 + stage + n);
        };
        if ((stage & 1) == 0) {
            execute_j0();
            execute_nonzero();
        } else {
            execute_nonzero();
            execute_j0();
        }
    }
    execute_partitioned(
        blst::MiddleProductRangePhase::AffineConversion,
        0,
        0x501 + n);
    execute_partitioned(
        blst::MiddleProductRangePhase::PointwisePairing,
        0,
        0x601 + n);
    execute_partitioned(
        blst::MiddleProductRangePhase::InverseBitReverse,
        0,
        0x701 + n);
    for (std::size_t stage = 0;
         stage < precomputed.inverse_plan.stage_left.size();
         ++stage) {
        execute_partitioned(
            blst::MiddleProductRangePhase::InverseButterflies,
            stage,
            0x800 + stage + n);
    }
    execute_partitioned(
        blst::MiddleProductRangePhase::Normalization,
        0,
        0x901 + n);

    for_shuffled_ranges(
        n,
        0xa01 + n,
        [&](const std::size_t begin, const std::size_t end) {
            blst::detail::middle_product_output_range_impl<Context>(
                precomputed,
                ranged_workspace,
                std::span<GT>(actual),
                begin,
                end);
        });

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(GT::equal(expected[i], actual[i]))
            << "n=" << n << ", i=" << i;
    }
}

template<class Context>
void expect_parallel_ranges_match_monolithic(const std::size_t n)
{
    const std::vector<G1Point> public_a = make_g1_data(n);
    const std::vector<G2Point> public_b = make_g2_data(n);
    const auto precomputed =
        blst::detail::precompute_middle_product_b_public<Context>(
            public_b, n, nullptr);
    const std::vector<typename Context::G1Point> native_a =
        blst::detail::convert_g1_inputs<Context>(public_a);
    blst::detail::MiddleProductWorkspaceImpl<Context> expected_workspace(
        precomputed.fft_size);
    blst::detail::MiddleProductWorkspaceImpl<Context> parallel_workspace(
        precomputed.fft_size);
    std::vector<GT> expected(n, GT::identity());
    std::vector<GT> actual(n, GT::identity());
    blst::detail::middle_product_impl<Context>(
        std::span<const typename Context::G1Point>(native_a),
        precomputed,
        std::span<GT>(expected),
        expected_workspace,
        nullptr);

    ASSERT_TRUE(for_parallel_ranges(
        precomputed.fft_size,
        [&](const std::size_t begin, const std::size_t end) {
            blst::detail::middle_product_input_range_impl<Context>(
                native_a, precomputed, parallel_workspace, begin, end);
        }));

    const auto execute_parallel = [&](
        const blst::MiddleProductRangePhase phase,
        const std::size_t stage = 0) {
        const std::size_t count =
            blst::detail::middle_product_range_work_count_impl<Context>(
                precomputed, phase, stage);
        return for_parallel_ranges(
            count,
            [&](const std::size_t begin, const std::size_t end) {
                blst::detail::middle_product_execute_range_impl<Context>(
                    precomputed,
                    parallel_workspace,
                    phase,
                    stage,
                    begin,
                    end);
            });
    };

    ASSERT_TRUE(execute_parallel(
        blst::MiddleProductRangePhase::ForwardBitReverse));
    for (std::size_t stage = 0;
         stage < precomputed.forward_plan.stage_work.size();
         ++stage) {
        ASSERT_TRUE(execute_parallel(
            blst::MiddleProductRangePhase::ForwardJ0, stage));
        ASSERT_TRUE(execute_parallel(
            blst::MiddleProductRangePhase::ForwardNonzero, stage));
    }
    ASSERT_TRUE(execute_parallel(
        blst::MiddleProductRangePhase::AffineConversion));
    ASSERT_TRUE(execute_parallel(
        blst::MiddleProductRangePhase::PointwisePairing));
    ASSERT_TRUE(execute_parallel(
        blst::MiddleProductRangePhase::InverseBitReverse));
    for (std::size_t stage = 0;
         stage < precomputed.inverse_plan.stage_left.size();
         ++stage) {
        ASSERT_TRUE(execute_parallel(
            blst::MiddleProductRangePhase::InverseButterflies, stage));
    }
    ASSERT_TRUE(execute_parallel(
        blst::MiddleProductRangePhase::Normalization));
    ASSERT_TRUE(for_parallel_ranges(
        n,
        [&](const std::size_t begin, const std::size_t end) {
            blst::detail::middle_product_output_range_impl<Context>(
                precomputed,
                parallel_workspace,
                std::span<GT>(actual),
                begin,
                end);
        }));

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(GT::equal(expected[i], actual[i])) << i;
    }
}

template<class Context>
auto fft_precompute_for_context(
    std::span<const G2Point> b,
    std::size_t n,
    blst::TimingStats* timings)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        return blst::precompute_middle_product_b_scalar(b, n, timings);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    } else {
        return blst::precompute_middle_product_b_avx(b, n, timings);
#endif
    }
}

template<class Context>
auto multipair_precompute_for_context(
    std::span<const G2Point> b,
    std::size_t n,
    blst::TimingStats* timings)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        return blst::precompute_middle_product_b_multipair_scalar(b, n, timings);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    } else {
        return blst::precompute_middle_product_b_multipair_avx(b, n, timings);
#endif
    }
}

template<class Context, class PrecomputedB>
void fft_middle_product_for_context(
    std::span<const G1Point> a,
    const PrecomputedB& precomputed_b,
    std::span<GT> out)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        blst::middle_product_scalar(a, precomputed_b, out);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    } else {
        blst::middle_product_avx(a, precomputed_b, out);
#endif
    }
}

template<class Context, class PrecomputedB>
void multipair_middle_product_for_context(
    std::span<const G1Point> a,
    const PrecomputedB& precomputed_b,
    std::span<GT> out)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        blst::middle_product_multipair_scalar(a, precomputed_b, out);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    } else {
        blst::middle_product_multipair_avx(a, precomputed_b, out);
#endif
    }
}

template<class Context>
void algorithm_middle_product_for_context(
    std::span<const G1Point> a,
    std::span<const G2Point> b,
    std::span<GT> out,
    MiddleProductAlgorithm algorithm,
    blst::TimingStats* timings = nullptr)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        blst::middle_product_scalar(a, b, out, algorithm, timings);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    } else if constexpr (std::is_same_v<Context, blst::MiddleProductAvxContext>) {
        blst::middle_product_avx(a, b, out, algorithm, timings);
#endif
    } else {
        static_assert(std::is_same_v<Context, blst::MiddleProductScalarContext>);
    }
}

template<class Context, class PrecomputedB>
void multipair_update_for_context(
    const G1Point& a0,
    const PrecomputedB& precomputed_b,
    std::span<GT> out)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        blst::MiddleProductWorkspaceMultipairScalar workspace;
        blst::update_multipair_scalar(a0, precomputed_b, out, workspace);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    } else {
        blst::MiddleProductWorkspaceMultipairAvx workspace;
        blst::update_multipair_avx(a0, precomputed_b, out, workspace);
#endif
    }
}

template<class Context>
void expect_algorithm_dispatch_matches_explicit(
    std::span<const G1Point> a,
    std::span<const G2Point> b,
    MiddleProductAlgorithm algorithm)
{
    const std::size_t n = a.size();
    std::vector<GT> dispatch_out(n, GT::identity());
    std::vector<GT> explicit_out(n, GT::identity());

    algorithm_middle_product_for_context<Context>(a, b, std::span<GT>(dispatch_out), algorithm);

    switch (algorithm) {
    case MiddleProductAlgorithm::Fft: {
        const auto fft_precomputed = fft_precompute_for_context<Context>(b, n, nullptr);
        fft_middle_product_for_context<Context>(a, fft_precomputed, std::span<GT>(explicit_out));
        break;
    }
    case MiddleProductAlgorithm::Multipair: {
        const auto multipair_precomputed = multipair_precompute_for_context<Context>(b, n, nullptr);
        multipair_middle_product_for_context<Context>(
            a,
            multipair_precomputed,
            std::span<GT>(explicit_out));
        break;
    }
    case MiddleProductAlgorithm::Auto:
    default:
        BLST_ASSERT(false && "Explicit dispatch helper expects Fft or Multipair");
    }

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(GT::equal(dispatch_out[i], explicit_out[i]))
            << "algorithm dispatch mismatch at n=" << n << ", i=" << i;
    }
}

template<class Context>
void expect_auto_dispatch_uses_expected_algorithm(
    std::span<const G1Point> a,
    std::span<const G2Point> b)
{
    const std::size_t n = a.size();
    const bool use_multipair = [&] {
        if constexpr (std::is_same_v<Context, MiddleProductScalarContext>) {
            return blst::middle_product_auto_uses_multipair_scalar(n);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
        } else {
            return blst::middle_product_auto_uses_multipair_avx(n);
#endif
        }
    }();
    const MiddleProductAlgorithm expected_algorithm = use_multipair
        ? MiddleProductAlgorithm::Multipair
        : MiddleProductAlgorithm::Fft;

    blst::TimingStats timings;
    std::vector<GT> out(n, GT::identity());

    algorithm_middle_product_for_context<Context>(
        a,
        b,
        std::span<GT>(out),
        MiddleProductAlgorithm::Auto,
        &timings);

    if (expected_algorithm == MiddleProductAlgorithm::Multipair) {
        EXPECT_GT(timings.get_ns("middle_product_multipair.total_ns"), 0);
        EXPECT_EQ(timings.get_ns("middle_product.total_ns"), 0);
    } else {
        EXPECT_GT(timings.get_ns("middle_product.total_ns"), 0);
        EXPECT_EQ(timings.get_ns("middle_product_multipair.total_ns"), 0);
    }
}

template<class Context>
void expect_multipair_matches_fft(
    std::span<const G1Point> a,
    std::span<const G2Point> b)
{
    const std::size_t n = a.size();
    const auto multipair_precomputed =
        multipair_precompute_for_context<Context>(b, n, nullptr);
    const auto fft_precomputed =
        fft_precompute_for_context<Context>(b, n, nullptr);

    std::vector<GT> multipair_out(n, GT::identity());
    std::vector<GT> fft_out(n, GT::identity());

    multipair_middle_product_for_context<Context>(
        a,
        multipair_precomputed,
        std::span<GT>(multipair_out));
    fft_middle_product_for_context<Context>(
        a,
        fft_precomputed,
        std::span<GT>(fft_out));

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(GT::equal(multipair_out[i], fft_out[i]))
            << "multipair mismatch at n=" << n << ", i=" << i;
    }
}

template<class Context>
void expect_multipair_update_and_combine_matches_full(
    std::span<const G1Point> a,
    std::span<const G2Point> b)
{
    const std::size_t n = a.size();
    const auto multipair_precomputed =
        multipair_precompute_for_context<Context>(b, n, nullptr);

    std::vector<GT> expected(n, GT::identity());
    std::vector<GT> mid(n, GT::identity());
    std::vector<GT> cnst(n, GT::identity());
    std::vector<GT> combined(n, GT::identity());

    std::vector<G1Point> a_without_constant(a.begin(), a.end());
    a_without_constant[0] = G1Point::infinity();

    multipair_middle_product_for_context<Context>(
        a,
        multipair_precomputed,
        std::span<GT>(expected));
    multipair_middle_product_for_context<Context>(
        std::span<const G1Point>(a_without_constant),
        multipair_precomputed,
        std::span<GT>(mid));
    multipair_update_for_context<Context>(a[0], multipair_precomputed, std::span<GT>(cnst));
    blst::combine(
        std::span<const GT>(mid),
        std::span<const GT>(cnst),
        std::span<GT>(combined));

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(GT::equal(expected[i], combined[i]))
            << "multipair update/combine mismatch at n=" << n << ", i=" << i;
    }
}

} // namespace

TEST(MiddleProductSizing, RejectsUnsupportedLengthsWithoutOverflow)
{
    EXPECT_THROW(blst::middle_product_fft_size(0), std::invalid_argument);
    EXPECT_EQ(
        blst::middle_product_fft_size(blst::middle_product_max_n),
        blst::middle_product_max_fft_size);
    EXPECT_THROW(
        blst::middle_product_fft_size(blst::middle_product_max_n + 1),
        std::length_error);
    EXPECT_THROW(
        blst::middle_product_fft_size(std::numeric_limits<std::size_t>::max()),
        std::length_error);
}

TEST(MiddleProductSizing, BasicOracleAndAdmissibleSizesAreChecked)
{
    constexpr std::size_t basic_max_n =
        (blst::middle_product_max_fft_size + 1) / 3;

    EXPECT_EQ(
        blst::middle_product_basic_fft_size(basic_max_n),
        blst::middle_product_max_fft_size);
    EXPECT_THROW(
        blst::middle_product_basic_fft_size(basic_max_n + 1),
        std::length_error);
    EXPECT_EQ(
        blst::fft_admissible_size(std::numeric_limits<std::size_t>::max(), 1024),
        1024);
    EXPECT_THROW(blst::fft_admissible_size(1, 0), std::invalid_argument);
    EXPECT_THROW(blst::fft_admissible_size(1, 3), std::invalid_argument);
}

TEST(MiddleProductPrecompute, InverseRootsReverseForwardPowers)
{
    EXPECT_TRUE(blst::make_ifft_roots(std::span<const blst::Fr>{}).empty());

    for (const std::size_t size : {std::size_t{1}, std::size_t{2}, std::size_t{8}}) {
        const std::vector<blst::Fr> roots = blst::make_fft_roots(size);
        const std::vector<blst::Fr> inv_roots = blst::make_ifft_roots(roots);

        ASSERT_EQ(inv_roots.size(), roots.size());
        EXPECT_TRUE(blst::Fr::equal(inv_roots[0], roots[0]));
        for (std::size_t i = 1; i < size; ++i) {
            EXPECT_TRUE(blst::Fr::equal(inv_roots[i], roots[size - i]));
            blst::Fr product = {};
            blst::Fr::mul(product, roots[i], inv_roots[i]);
            EXPECT_TRUE(blst::Fr::equal(product, blst::Fr::one()));
        }
    }
}

TEST(MiddleProductPrecompute, ScalarPlansRemainStableAcrossEvaluations)
{
    constexpr std::size_t n = 3;
    const std::vector<G1Point> a = make_g1_data(n);
    const std::vector<G2Point> b = make_g2_data(n);
    const auto precomputed = blst::precompute_middle_product_b_scalar(b, n);

    EXPECT_EQ(precomputed.forward_plan.n, precomputed.fft_size);
    EXPECT_EQ(
        precomputed.forward_plan.xn,
        blst::fft_admissible_size(n, precomputed.fft_size));
    EXPECT_EQ(precomputed.inverse_plan.n, precomputed.fft_size);
    EXPECT_EQ(precomputed.inverse_plan.start, n);
    EXPECT_EQ(precomputed.inverse_plan.length, n);

    const auto* const forward_storage = precomputed.forward_plan.left_nonzero.data();
    const auto* const inverse_storage = precomputed.inverse_plan.left_indices.data();
    std::vector<GT> first(n, GT::identity());
    std::vector<GT> second(n, GT::identity());
    blst::MiddleProductWorkspaceScalar first_workspace(precomputed.fft_size);
    blst::MiddleProductWorkspaceScalar second_workspace(precomputed.fft_size);

    blst::middle_product_scalar(a, precomputed, first, first_workspace);
    blst::middle_product_scalar(a, precomputed, second, second_workspace);

    EXPECT_EQ(precomputed.forward_plan.left_nonzero.data(), forward_storage);
    EXPECT_EQ(precomputed.inverse_plan.left_indices.data(), inverse_storage);
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(GT::equal(first[i], second[i]));
    }
}

TEST(MiddleProductWorkspace, ScalarCapacityIsAlignedStableAndReusableAcrossAnchors)
{
    constexpr std::size_t small_n = 5;
    constexpr std::size_t large_n = 17;
    const std::vector<G1Point> small_a = make_g1_data(small_n);
    const std::vector<G2Point> small_b = make_g2_data(small_n);
    const std::vector<G1Point> large_a = make_g1_data(large_n);
    const std::vector<G2Point> large_b = make_g2_data(large_n);
    const auto small_precomputed =
        blst::precompute_middle_product_b_scalar(small_b, small_n);
    const auto large_precomputed =
        blst::precompute_middle_product_b_scalar(large_b, large_n);
    ASSERT_LT(small_precomputed.fft_size, large_precomputed.fft_size);

    blst::MiddleProductWorkspaceScalar workspace(large_precomputed.fft_size);
    EXPECT_EQ(workspace.capacity(), large_precomputed.fft_size);
    ASSERT_NE(workspace.a_fft.data(), nullptr);
    ASSERT_NE(workspace.c_fft.data(), nullptr);
    EXPECT_EQ(
        reinterpret_cast<std::uintptr_t>(workspace.a_fft.data()) % 64,
        std::uintptr_t{0});
    EXPECT_EQ(
        reinterpret_cast<std::uintptr_t>(workspace.c_fft.data()) % 64,
        std::uintptr_t{0});
    const auto* const a_address = workspace.a_fft.data();
    const auto* const c_address = workspace.c_fft.data();

    std::vector<GT> small_expected(small_n, GT::identity());
    std::vector<GT> small_actual(small_n, GT::identity());
    std::vector<GT> large_out(large_n, GT::identity());
    blst::middle_product_scalar(small_a, small_precomputed, small_expected);
    blst::middle_product_scalar(
        large_a, large_precomputed, large_out, workspace);
    blst::middle_product_scalar(
        small_a, small_precomputed, small_actual, workspace);

    EXPECT_EQ(workspace.a_fft.data(), a_address);
    EXPECT_EQ(workspace.c_fft.data(), c_address);
    for (std::size_t i = 0; i < small_n; ++i) {
        EXPECT_TRUE(GT::equal(small_expected[i], small_actual[i])) << i;
    }
}

TEST(MiddleProductWorkspace, CapacityMismatchFailsBeforeOutputPublication)
{
    constexpr std::size_t small_n = 5;
    constexpr std::size_t large_n = 17;
    const std::vector<G1Point> large_a = make_g1_data(large_n);
    const std::vector<G2Point> large_b = make_g2_data(large_n);
    const auto large_precomputed =
        blst::precompute_middle_product_b_scalar(large_b, large_n);
    const std::size_t small_capacity = blst::middle_product_fft_size(small_n);
    ASSERT_LT(small_capacity, large_precomputed.fft_size);
    blst::MiddleProductWorkspaceScalar workspace(small_capacity);
    const std::vector<GT> canary(large_n, GT::generator());
    std::vector<GT> output = canary;

    EXPECT_THROW(
        blst::middle_product_scalar(
            large_a, large_precomputed, output, workspace),
        std::length_error);
    for (std::size_t i = 0; i < output.size(); ++i) {
        EXPECT_TRUE(GT::equal(output[i], canary[i])) << i;
    }

    EXPECT_THROW(blst::MiddleProductWorkspaceScalar(0), std::invalid_argument);
    EXPECT_THROW(blst::MiddleProductWorkspaceScalar(3), std::invalid_argument);
}

TEST(MiddleProductRanges, ScalarRandomPartitionsMatchMonolithic)
{
    for (const std::size_t n : {std::size_t{1}, std::size_t{5}, std::size_t{17}}) {
        expect_partitioned_ranges_match_monolithic<
            blst::MiddleProductScalarContext>(n);
    }
}

TEST(MiddleProductRanges, ScalarConcurrentDisjointRangesMatchMonolithic)
{
    expect_parallel_ranges_match_monolithic<
        blst::MiddleProductScalarContext>(17);
}

TEST(MiddleProductRanges, RejectsInvalidPhaseStageRangeAndCapacity)
{
    constexpr std::size_t n = 17;
    const std::vector<G1Point> public_a = make_g1_data(n);
    const std::vector<G2Point> public_b = make_g2_data(n);
    const auto precomputed = blst::detail::precompute_middle_product_b_public<
        blst::MiddleProductScalarContext>(public_b, n, nullptr);
    const std::vector<blst::MiddleProductScalarContext::G1Point> native_a =
        blst::detail::convert_g1_inputs<blst::MiddleProductScalarContext>(
            public_a);
    blst::MiddleProductWorkspaceScalar workspace(precomputed.fft_size);

    EXPECT_NO_THROW(
        blst::detail::middle_product_execute_range_impl<
            blst::MiddleProductScalarContext>(
                precomputed,
                workspace,
                blst::MiddleProductRangePhase::ForwardBitReverse,
                0,
                0,
                0));
    EXPECT_THROW(
        (void)blst::detail::middle_product_range_work_count_impl<
            blst::MiddleProductScalarContext>(
                precomputed,
                blst::MiddleProductRangePhase::ForwardBitReverse,
                1),
        std::invalid_argument);
    EXPECT_THROW(
        (void)blst::detail::middle_product_range_work_count_impl<
            blst::MiddleProductScalarContext>(
                precomputed,
                blst::MiddleProductRangePhase::ForwardJ0,
                precomputed.forward_plan.stage_work.size()),
        std::invalid_argument);
    const std::size_t swaps =
        blst::detail::middle_product_range_work_count_impl<
            blst::MiddleProductScalarContext>(
                precomputed,
                blst::MiddleProductRangePhase::ForwardBitReverse,
                0);
    EXPECT_THROW(
        blst::detail::middle_product_execute_range_impl<
            blst::MiddleProductScalarContext>(
                precomputed,
                workspace,
                blst::MiddleProductRangePhase::ForwardBitReverse,
                0,
                swaps,
                swaps + 1),
        std::invalid_argument);
    EXPECT_THROW(
        blst::detail::middle_product_execute_range_impl<
            blst::MiddleProductScalarContext>(
                precomputed,
                workspace,
                blst::MiddleProductRangePhase::InputLoad,
                0,
                0,
                0),
        std::invalid_argument);

    blst::MiddleProductWorkspaceScalar undersized(
        blst::middle_product_fft_size(5));
    EXPECT_THROW(
        blst::detail::middle_product_input_range_impl<
            blst::MiddleProductScalarContext>(
                native_a, precomputed, undersized, 0, 0),
        std::length_error);
}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
TEST(MiddleProductWorkspace, AvxCapacityAlignsEveryNativeBuffer)
{
    constexpr std::size_t small_n = 5;
    constexpr std::size_t large_n = 17;
    const std::vector<G1Point> small_a = make_g1_data(small_n);
    const std::vector<G2Point> small_b = make_g2_data(small_n);
    const std::vector<G1Point> large_a = make_g1_data(large_n);
    const std::vector<G2Point> large_b = make_g2_data(large_n);
    const auto small_precomputed =
        blst::precompute_middle_product_b_avx(small_b, small_n);
    const auto large_precomputed =
        blst::precompute_middle_product_b_avx(large_b, large_n);
    blst::MiddleProductWorkspaceAvx workspace(large_precomputed.fft_size);

    const auto expect_aligned = [](const auto* const address) {
        ASSERT_NE(address, nullptr);
        EXPECT_EQ(
            reinterpret_cast<std::uintptr_t>(address) % 64,
            std::uintptr_t{0});
    };
    expect_aligned(workspace.a_fft.data());
    expect_aligned(workspace.a_fft_affine.data());
    expect_aligned(workspace.affine_prefix_products.data());
    expect_aligned(workspace.c_fft.data());
    const auto* const a_address = workspace.a_fft.data();
    const auto* const affine_address = workspace.a_fft_affine.data();
    const auto* const prefix_address = workspace.affine_prefix_products.data();
    const auto* const c_address = workspace.c_fft.data();

    std::vector<GT> small_expected(small_n, GT::identity());
    std::vector<GT> small_actual(small_n, GT::identity());
    std::vector<GT> large_out(large_n, GT::identity());
    blst::middle_product_avx(small_a, small_precomputed, small_expected);
    blst::middle_product_avx(
        large_a, large_precomputed, large_out, workspace);
    blst::middle_product_avx(
        small_a, small_precomputed, small_actual, workspace);

    EXPECT_EQ(workspace.a_fft.data(), a_address);
    EXPECT_EQ(workspace.a_fft_affine.data(), affine_address);
    EXPECT_EQ(workspace.affine_prefix_products.data(), prefix_address);
    EXPECT_EQ(workspace.c_fft.data(), c_address);
    for (std::size_t i = 0; i < small_n; ++i) {
        EXPECT_TRUE(GT::equal(small_expected[i], small_actual[i])) << i;
    }
}


TEST(MiddleProductRanges, AvxRandomPartitionsMatchMonolithic)
{
    for (const std::size_t n : {std::size_t{1}, std::size_t{5}, std::size_t{17}}) {
        expect_partitioned_ranges_match_monolithic<
            blst::MiddleProductAvxContext>(n);
    }
}


TEST(MiddleProductRanges, AvxConcurrentDisjointRangesMatchMonolithic)
{
    expect_parallel_ranges_match_monolithic<
        blst::MiddleProductAvxContext>(17);
}
#endif

TEST(MiddleProductMultipair, ScalarMatchesFft)
{
    const std::vector<std::size_t> sizes = { 3, 5, 9, 16 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_multipair_matches_fft<MiddleProductScalarContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b));
    }
}

TEST(MiddleProductMultipair, ScalarSparseWorkspaceReservesBeforePublication)
{
    constexpr std::size_t n = 10;
    const std::vector<G1Point> a = make_g1_data(n);
    const std::vector<G2Point> b = make_g2_data(n);
    const auto precomputed =
        blst::precompute_middle_product_b_multipair_scalar(b, n);
    blst::MiddleProductWorkspaceMultipairScalar workspace;
    std::vector<GT> out(n, GT::identity());

    blst::middle_product_multipair_scalar(a, precomputed, out, workspace);

    std::size_t nonzero_a = 0;
    for (const G1Point& point : a) {
        nonzero_a += (G1Point::is_inf(point) == 0) ? 1 : 0;
    }
    ASSERT_LT(nonzero_a, n);
    EXPECT_GE(workspace.line_batch_scalar.capacity(), nonzero_a);
    EXPECT_GE(workspace.line_scale_batch_scalar.capacity(), nonzero_a);
}

TEST(MiddleProductMultipair, ScalarUpdateAndCombineReconstructsFullResult)
{
    const std::vector<std::size_t> sizes = { 4, 6, 10, 18 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_multipair_update_and_combine_matches_full<MiddleProductScalarContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b));
    }
}

TEST(MiddleProductAlgorithm, ScalarFftMatchesExplicitFft)
{
    const std::vector<std::size_t> sizes = { 3, 5, 9, 16, 33 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_algorithm_dispatch_matches_explicit<MiddleProductScalarContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b),
            MiddleProductAlgorithm::Fft);
    }
}

TEST(MiddleProductAlgorithm, ScalarMultipairMatchesExplicitMultipair)
{
    const std::vector<std::size_t> sizes = { 3, 5, 9, 16, 24 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_algorithm_dispatch_matches_explicit<MiddleProductScalarContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b),
            MiddleProductAlgorithm::Multipair);
    }
}

TEST(MiddleProductAlgorithm, ScalarAutoUsesExpectedAlgorithm)
{
    const std::vector<std::size_t> sizes = { 3, 32 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_auto_dispatch_uses_expected_algorithm<MiddleProductScalarContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b));
    }
}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
template<class Point>
Point avx_point(uint8_t scalar)
{
    std::array<blst::byte, 1> scalar_bytes = {scalar};
    Point out = {};
    Point::mult(out, Point::generator(), scalar_bytes.data(), 8);
    return out;
}

TEST(AvxBatchAffine, G1MatchesIndividualConversion)
{
    using Affine = blst::avx::G1AffinePoint;
    using Point = blst::avx::G1Point;

    const std::vector<Point> points = {
        Point::infinity(),
        avx_point<Point>(1),
        avx_point<Point>(2),
        Point::infinity(),
        avx_point<Point>(17),
        avx_point<Point>(255),
    };
    std::vector<Affine> actual(points.size());
    std::vector<Point::Field_type> prefix_products(points.size());
    Point::to_affine_batch(
        std::span<Affine>(actual),
        std::span<const Point>(points),
        std::span<Point::Field_type>(prefix_products));

    for (std::size_t i = 0; i < points.size(); ++i) {
        Affine expected = {};
        Point::to_affine(expected, points[i]);
        EXPECT_TRUE(Affine::equal(expected, actual[i])) << "G1 index " << i;
    }

    std::vector<Point> all_infinity(4, Point::infinity());
    actual.assign(all_infinity.size(), Affine::generator());
    prefix_products.resize(all_infinity.size());
    Point::to_affine_batch(
        std::span<Affine>(actual),
        std::span<const Point>(all_infinity),
        std::span<Point::Field_type>(prefix_products));
    for (const Affine& point : actual) {
        EXPECT_TRUE(Affine::is_inf(point));
    }

    Point::to_affine_batch({}, {}, {});
}

TEST(AvxBatchAffine, G2MatchesIndividualConversion)
{
    using Affine = blst::avx::G2AffinePoint;
    using Point = blst::avx::G2Point;

    const std::vector<Point> points = {
        Point::infinity(),
        avx_point<Point>(1),
        avx_point<Point>(2),
        Point::infinity(),
        avx_point<Point>(17),
        avx_point<Point>(255),
    };
    std::vector<Affine> actual(points.size());
    std::vector<Point::Field_type> prefix_products(points.size());
    Point::to_affine_batch(
        std::span<Affine>(actual),
        std::span<const Point>(points),
        std::span<Point::Field_type>(prefix_products));

    for (std::size_t i = 0; i < points.size(); ++i) {
        Affine expected = {};
        Point::to_affine(expected, points[i]);
        EXPECT_TRUE(Affine::equal(expected, actual[i])) << "G2 index " << i;
    }

    std::vector<Point> all_infinity(4, Point::infinity());
    actual.assign(all_infinity.size(), Affine::generator());
    prefix_products.resize(all_infinity.size());
    Point::to_affine_batch(
        std::span<Affine>(actual),
        std::span<const Point>(all_infinity),
        std::span<Point::Field_type>(prefix_products));
    for (const Affine& point : actual) {
        EXPECT_TRUE(Affine::is_inf(point));
    }

    Point::to_affine_batch({}, {}, {});
}

TEST(MiddleProductMultipair, AvxPortableMatchesFft)
{
    const std::vector<std::size_t> sizes = { 3, 5, 9, 16 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_multipair_matches_fft<blst::MiddleProductAvxContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b));
    }
}

TEST(MiddleProductAlgorithm, ScalarAutoPolicyMatchesBenchmarkRanges)
{
    const std::pair<std::size_t, bool> cases[] = {
        { 1, true },
        { 31, true },
        { 32, false },
        { 33, true },
        { 49, true },
        { 50, false },
        { 64, false },
        { 65, true },
        { 76, true },
        { 77, false },
        { 128, false },
        { 129, false },
    };

    for (const auto& [n, expected] : cases) {
        EXPECT_EQ(blst::middle_product_auto_uses_multipair_scalar(n), expected)
            << "unexpected scalar Auto policy at n=" << n;
    }
}

TEST(MiddleProductMultipair, AvxWorkspaceReusesFinalExpScratch)
{
    constexpr std::size_t n = 3;
    const std::vector<G1Point> a = make_g1_data(n);
    const std::vector<G2Point> b = make_g2_data(n);
    const auto precomputed = blst::precompute_middle_product_b_multipair_avx(b, n);
    blst::MiddleProductWorkspaceMultipairAvx workspace;
    std::vector<GT> first(n, GT::identity());
    std::vector<GT> second(n, GT::identity());

    blst::middle_product_multipair_avx(a, precomputed, first, workspace);
    EXPECT_TRUE(workspace.final_exp_scratch.has_capacity(n));

    blst::middle_product_multipair_avx(a, precomputed, second, workspace);

    EXPECT_TRUE(workspace.final_exp_scratch.has_capacity(n));
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(GT::equal(first[i], second[i]));
    }
}

TEST(MiddleProductMultipair, AvxPortableUpdateAndCombineReconstructsFullResult)
{
    const std::vector<std::size_t> sizes = { 4, 6, 10, 18 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_multipair_update_and_combine_matches_full<blst::MiddleProductAvxContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b));
    }
}

TEST(MiddleProductAlgorithm, AvxPortableFftMatchesExplicitFft)
{
    const std::vector<std::size_t> sizes = {
        3,
        16,
        32,
    };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_algorithm_dispatch_matches_explicit<blst::MiddleProductAvxContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b),
            MiddleProductAlgorithm::Fft);
    }
}

TEST(MiddleProductAlgorithm, AvxPreparedLinesPreserveInfinitySemantics)
{
    constexpr std::size_t n = 8;
    const std::vector<G1Point> finite_a = make_g1_data(n);
    const std::vector<G2Point> zero_b(2 * n, G2Point::infinity());
    const auto zero_b_precomputed =
        fft_precompute_for_context<blst::MiddleProductAvxContext>(zero_b, n, nullptr);

    ASSERT_EQ(zero_b_precomputed.prepared_lines.size(), zero_b_precomputed.fft_size);
    ASSERT_EQ(zero_b_precomputed.b_fft_is_zero.size(), zero_b_precomputed.fft_size);
    for (uint8_t is_zero : zero_b_precomputed.b_fft_is_zero) {
        EXPECT_NE(is_zero, 0);
    }

    std::vector<GT> out(n, GT::identity());
    fft_middle_product_for_context<blst::MiddleProductAvxContext>(
        finite_a,
        zero_b_precomputed,
        out);
    for (const GT& value : out) {
        EXPECT_TRUE(GT::equal(value, GT::identity()));
    }

    const std::vector<G1Point> zero_a(n, G1Point::infinity());
    const std::vector<G2Point> finite_b = make_g2_data(n);
    const auto finite_b_precomputed =
        fft_precompute_for_context<blst::MiddleProductAvxContext>(finite_b, n, nullptr);
    fft_middle_product_for_context<blst::MiddleProductAvxContext>(
        zero_a,
        finite_b_precomputed,
        out);
    for (const GT& value : out) {
        EXPECT_TRUE(GT::equal(value, GT::identity()));
    }
}

TEST(MiddleProductAlgorithm, AvxPortableMultipairMatchesExplicitMultipair)
{
    const std::vector<std::size_t> sizes = {
        3,
        16,
        24,
    };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_algorithm_dispatch_matches_explicit<blst::MiddleProductAvxContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b),
            MiddleProductAlgorithm::Multipair);
    }
}

TEST(MiddleProductAlgorithm, AvxPortableAutoUsesExpectedAlgorithm)
{
    const std::vector<std::size_t> sizes = { 3, 32 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_auto_dispatch_uses_expected_algorithm<blst::MiddleProductAvxContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b));
    }
}

TEST(MiddleProductAlgorithm, AvxAutoPolicyMatchesBenchmarkRanges)
{
    const std::pair<std::size_t, bool> cases[] = {
        { 1, true },
        { 31, true },
        { 32, false },
        { 33, true },
        { 50, true },
        { 51, false },
        { 64, false },
        { 65, true },
        { 76, true },
        { 77, false },
        { 128, false },
        { 129, false },
    };

    for (const auto& [n, expected] : cases) {
        EXPECT_EQ(blst::middle_product_auto_uses_multipair_avx(n), expected)
            << "unexpected AVX Auto policy at n=" << n;
    }
}
#endif

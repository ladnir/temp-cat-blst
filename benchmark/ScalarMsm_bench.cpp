#include <benchmark/benchmark.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#include "bindings/blst.h"

namespace {

blst_fr fr_from_u64(std::uint64_t value)
{
    const std::array<std::uint64_t, 4> words = {value, 0, 0, 0};
    blst_fr out = {};
    blst_fr_from_uint64(&out, words.data());
    return out;
}

std::vector<blst_fr> make_scalars(std::size_t count)
{
    const blst_fr multiplier = fr_from_u64(6364136223846793005ULL);
    const blst_fr increment = fr_from_u64(1442695040888963407ULL);
    blst_fr state = fr_from_u64(17);
    std::vector<blst_fr> scalars(count);
    for (std::size_t i = 0; i < count; ++i) {
        scalars[i] = state;
        blst_fr_mul(&state, &state, &multiplier);
        blst_fr_add(&state, &state, &increment);
    }
    return scalars;
}

std::vector<blst_p1_affine> make_g1_points(std::size_t count)
{
    std::vector<blst_p1_affine> affine(count);
    if (count == 0) {
        return affine;
    }

    std::vector<blst_p1> projective(count);
    const blst_p1 generator = *blst_p1_generator();
    blst_p1 current = generator;
    for (std::size_t i = 0; i < count; ++i) {
        projective[i] = current;
        blst_p1 next = {};
        blst_p1_add_or_double(&next, &current, &generator);
        current = next;
    }
    const blst_p1* descriptors[] = {projective.data(), nullptr};
    blst_p1s_to_affine(affine.data(), descriptors, count);
    return affine;
}

std::vector<blst_p2_affine> make_g2_points(std::size_t count)
{
    std::vector<blst_p2_affine> affine(count);
    if (count == 0) {
        return affine;
    }

    std::vector<blst_p2> projective(count);
    const blst_p2 generator = *blst_p2_generator();
    blst_p2 current = generator;
    for (std::size_t i = 0; i < count; ++i) {
        projective[i] = current;
        blst_p2 next = {};
        blst_p2_add_or_double(&next, &current, &generator);
        current = next;
    }
    const blst_p2* descriptors[] = {projective.data(), nullptr};
    blst_p2s_to_affine(affine.data(), descriptors, count);
    return affine;
}

std::vector<limb_t> scratch_words(std::size_t bytes)
{
    if (bytes == 0 || bytes % sizeof(limb_t) != 0) {
        throw std::runtime_error("invalid scalar MSM scratch size");
    }
    return std::vector<limb_t>(bytes / sizeof(limb_t));
}

struct G1Inputs
{
    explicit G1Inputs(std::size_t count)
        : points(make_g1_points(count)),
          scalars(make_scalars(count)),
          canonical(count),
          canonical_scratch(scratch_words(
              blst_p1s_mult_pippenger_scratch_sizeof(count))),
          native_scratch(scratch_words(
              blst_p1s_mult_pippenger_fr_scratch_sizeof(count))),
          selected_scratch(scratch_words(
              blst_p1s_mult_fr_scratch_sizeof(count)))
    {
        if (!blst_scalars_from_fr(canonical.data(), scalars.data(), count)) {
            throw std::runtime_error("G1 scalar preparation failed");
        }
        const blst_p1_affine* point_descriptors[] = {points.data(), nullptr};
        const byte* scalar_descriptors[] = {
            reinterpret_cast<const byte*>(canonical.data()), nullptr};
        blst_p1 canonical_result = {};
        blst_p1 native_result = {};
        blst_p1s_mult_pippenger(
            &canonical_result, point_descriptors, count, scalar_descriptors,
            255, canonical_scratch.data());
        if (!blst_p1s_mult_pippenger_fr(
                &native_result, points.data(), scalars.data(), count,
                native_scratch.data(), native_scratch.size() * sizeof(limb_t)) ||
            !blst_p1_is_equal(&canonical_result, &native_result)) {
            throw std::runtime_error("G1 native-Fr and canonical MSM differ");
        }
    }

    std::vector<blst_p1_affine> points;
    std::vector<blst_fr> scalars;
    std::vector<blst_scalar> canonical;
    std::vector<limb_t> canonical_scratch;
    std::vector<limb_t> native_scratch;
    std::vector<limb_t> selected_scratch;
};

struct G2Inputs
{
    explicit G2Inputs(std::size_t count)
        : points(make_g2_points(count)),
          scalars(make_scalars(count)),
          canonical(count),
          canonical_scratch(scratch_words(
              blst_p2s_mult_pippenger_scratch_sizeof(count))),
          native_scratch(scratch_words(
              blst_p2s_mult_pippenger_fr_scratch_sizeof(count))),
          selected_scratch(scratch_words(
              blst_p2s_mult_fr_scratch_sizeof(count)))
    {
        if (!blst_scalars_from_fr(canonical.data(), scalars.data(), count)) {
            throw std::runtime_error("G2 scalar preparation failed");
        }
        const blst_p2_affine* point_descriptors[] = {points.data(), nullptr};
        const byte* scalar_descriptors[] = {
            reinterpret_cast<const byte*>(canonical.data()), nullptr};
        blst_p2 canonical_result = {};
        blst_p2 native_result = {};
        blst_p2s_mult_pippenger(
            &canonical_result, point_descriptors, count, scalar_descriptors,
            255, canonical_scratch.data());
        if (!blst_p2s_mult_pippenger_fr(
                &native_result, points.data(), scalars.data(), count,
                native_scratch.data(), native_scratch.size() * sizeof(limb_t)) ||
            !blst_p2_is_equal(&canonical_result, &native_result)) {
            throw std::runtime_error("G2 native-Fr and canonical MSM differ");
        }
    }

    std::vector<blst_p2_affine> points;
    std::vector<blst_fr> scalars;
    std::vector<blst_scalar> canonical;
    std::vector<limb_t> canonical_scratch;
    std::vector<limb_t> native_scratch;
    std::vector<limb_t> selected_scratch;
};

void set_counters(
    benchmark::State& state,
    std::size_t count,
    std::size_t scratch,
    std::size_t prepared_scalars)
{
    state.SetItemsProcessed(
        static_cast<std::int64_t>(state.iterations() * count));
    state.counters["points"] = static_cast<double>(count);
    state.counters["prepared_scalar_bytes"] =
        static_cast<double>(prepared_scalars);
    state.counters["scratch_bytes"] = static_cast<double>(scratch);
}

void BM_ScalarMsm_FrToCanonical(benchmark::State& state)
{
    const auto count = static_cast<std::size_t>(state.range(0));
    const auto scalars = make_scalars(count);
    std::vector<blst_scalar> canonical(count);
    for (auto _ : state) {
        const bool success =
            blst_scalars_from_fr(canonical.data(), scalars.data(), count);
        benchmark::DoNotOptimize(success);
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(
        static_cast<std::int64_t>(state.iterations() * count));
}

void BM_ScalarMsm_G1Canonical(benchmark::State& state)
{
    const auto count = static_cast<std::size_t>(state.range(0));
    G1Inputs input(count);
    const blst_p1_affine* point_descriptors[] = {input.points.data(), nullptr};
    const byte* scalar_descriptors[] = {
        reinterpret_cast<const byte*>(input.canonical.data()), nullptr};
    blst_p1 output = {};
    for (auto _ : state) {
        blst_p1s_mult_pippenger(
            &output, point_descriptors, count, scalar_descriptors, 255,
            input.canonical_scratch.data());
        benchmark::DoNotOptimize(output);
        benchmark::ClobberMemory();
    }
    set_counters(
        state, count, input.canonical_scratch.size() * sizeof(limb_t),
        input.canonical.size() * sizeof(blst_scalar));
}

void BM_ScalarMsm_G1NativeFr(benchmark::State& state)
{
    const auto count = static_cast<std::size_t>(state.range(0));
    G1Inputs input(count);
    blst_p1 output = {};
    for (auto _ : state) {
        const bool success = blst_p1s_mult_pippenger_fr(
            &output, input.points.data(), input.scalars.data(), count,
            input.native_scratch.data(),
            input.native_scratch.size() * sizeof(limb_t));
        benchmark::DoNotOptimize(success);
        benchmark::DoNotOptimize(output);
        benchmark::ClobberMemory();
    }
    set_counters(
        state, count, input.native_scratch.size() * sizeof(limb_t), 0);
}

void BM_ScalarMsm_G1SelectedFr(benchmark::State& state)
{
    const auto count = static_cast<std::size_t>(state.range(0));
    G1Inputs input(count);
    blst_p1 output = {};
    for (auto _ : state) {
        const bool success = blst_p1s_mult_fr(
            &output, input.points.data(), input.scalars.data(), count,
            input.selected_scratch.data(),
            input.selected_scratch.size() * sizeof(limb_t));
        benchmark::DoNotOptimize(success);
        benchmark::DoNotOptimize(output);
        benchmark::ClobberMemory();
    }
    set_counters(
        state, count, input.selected_scratch.size() * sizeof(limb_t), 0);
}

void BM_ScalarMsm_G1ExplicitWindow(benchmark::State& state)
{
    const auto count = static_cast<std::size_t>(state.range(0));
    const auto window = static_cast<std::size_t>(state.range(1));
    G1Inputs input(count);
    auto window_scratch = scratch_words(
        blst_p1s_mult_pippenger_window_scratch_sizeof(window));
    const blst_p1_affine* point_descriptors[] = {input.points.data(), nullptr};
    const byte* scalar_descriptors[] = {
        reinterpret_cast<const byte*>(input.canonical.data()), nullptr};
    blst_p1 expected = {};
    blst_p1 output = {};
    blst_p1s_mult_pippenger(
        &expected, point_descriptors, count, scalar_descriptors, 255,
        input.canonical_scratch.data());
    blst_p1s_mult_pippenger_window(
        &output, point_descriptors, count, scalar_descriptors, 255,
        window_scratch.data(), window);
    if (!blst_p1_is_equal(&expected, &output)) {
        throw std::runtime_error("G1 explicit-window MSM differs");
    }

    for (auto _ : state) {
        blst_p1s_mult_pippenger_window(
            &output, point_descriptors, count, scalar_descriptors, 255,
            window_scratch.data(), window);
        benchmark::DoNotOptimize(output);
        benchmark::ClobberMemory();
    }
    set_counters(
        state, count, window_scratch.size() * sizeof(limb_t),
        input.canonical.size() * sizeof(blst_scalar));
    state.counters["window"] = static_cast<double>(window);
}

void BM_ScalarMsm_G2Canonical(benchmark::State& state)
{
    const auto count = static_cast<std::size_t>(state.range(0));
    G2Inputs input(count);
    const blst_p2_affine* point_descriptors[] = {input.points.data(), nullptr};
    const byte* scalar_descriptors[] = {
        reinterpret_cast<const byte*>(input.canonical.data()), nullptr};
    blst_p2 output = {};
    for (auto _ : state) {
        blst_p2s_mult_pippenger(
            &output, point_descriptors, count, scalar_descriptors, 255,
            input.canonical_scratch.data());
        benchmark::DoNotOptimize(output);
        benchmark::ClobberMemory();
    }
    set_counters(
        state, count, input.canonical_scratch.size() * sizeof(limb_t),
        input.canonical.size() * sizeof(blst_scalar));
}

void BM_ScalarMsm_G2NativeFr(benchmark::State& state)
{
    const auto count = static_cast<std::size_t>(state.range(0));
    G2Inputs input(count);
    blst_p2 output = {};
    for (auto _ : state) {
        const bool success = blst_p2s_mult_pippenger_fr(
            &output, input.points.data(), input.scalars.data(), count,
            input.native_scratch.data(),
            input.native_scratch.size() * sizeof(limb_t));
        benchmark::DoNotOptimize(success);
        benchmark::DoNotOptimize(output);
        benchmark::ClobberMemory();
    }
    set_counters(
        state, count, input.native_scratch.size() * sizeof(limb_t), 0);
}

void BM_ScalarMsm_G2SelectedFr(benchmark::State& state)
{
    const auto count = static_cast<std::size_t>(state.range(0));
    G2Inputs input(count);
    blst_p2 output = {};
    for (auto _ : state) {
        const bool success = blst_p2s_mult_fr(
            &output, input.points.data(), input.scalars.data(), count,
            input.selected_scratch.data(),
            input.selected_scratch.size() * sizeof(limb_t));
        benchmark::DoNotOptimize(success);
        benchmark::DoNotOptimize(output);
        benchmark::ClobberMemory();
    }
    set_counters(
        state, count, input.selected_scratch.size() * sizeof(limb_t), 0);
}

void BM_ScalarMsm_G2ExplicitWindow(benchmark::State& state)
{
    const auto count = static_cast<std::size_t>(state.range(0));
    const auto window = static_cast<std::size_t>(state.range(1));
    G2Inputs input(count);
    auto window_scratch = scratch_words(
        blst_p2s_mult_pippenger_window_scratch_sizeof(window));
    const blst_p2_affine* point_descriptors[] = {input.points.data(), nullptr};
    const byte* scalar_descriptors[] = {
        reinterpret_cast<const byte*>(input.canonical.data()), nullptr};
    blst_p2 expected = {};
    blst_p2 output = {};
    blst_p2s_mult_pippenger(
        &expected, point_descriptors, count, scalar_descriptors, 255,
        input.canonical_scratch.data());
    blst_p2s_mult_pippenger_window(
        &output, point_descriptors, count, scalar_descriptors, 255,
        window_scratch.data(), window);
    if (!blst_p2_is_equal(&expected, &output)) {
        throw std::runtime_error("G2 explicit-window MSM differs");
    }

    for (auto _ : state) {
        blst_p2s_mult_pippenger_window(
            &output, point_descriptors, count, scalar_descriptors, 255,
            window_scratch.data(), window);
        benchmark::DoNotOptimize(output);
        benchmark::ClobberMemory();
    }
    set_counters(
        state, count, window_scratch.size() * sizeof(limb_t),
        input.canonical.size() * sizeof(blst_scalar));
    state.counters["window"] = static_cast<double>(window);
}

bool register_scalar_msm_benchmarks()
{
    const auto environment_log_size = [](const char* name, int fallback) {
        const char* value = std::getenv(name);
        if (value == nullptr) {
            return fallback;
        }
        const int parsed = std::stoi(value);
        if (parsed < 0 || parsed > 20) {
            throw std::runtime_error(std::string(name) + " must lie in [0, 20]");
        }
        return parsed;
    };
    const int min_log_size =
        environment_log_size("HYPERCAT_MSM_MIN_LOG_SIZE", 0);
    const int max_log_size =
        environment_log_size("HYPERCAT_MSM_MAX_LOG_SIZE", 18);
    if (min_log_size > max_log_size) {
        throw std::runtime_error("MSM minimum log size exceeds maximum");
    }
    const char* group_value = std::getenv("HYPERCAT_MSM_GROUP");
    const std::string group = group_value == nullptr ? "all" : group_value;
    if (group != "all" && group != "g1" && group != "g2") {
        throw std::runtime_error("HYPERCAT_MSM_GROUP must be all, g1, or g2");
    }
    const int min_window =
        environment_log_size("HYPERCAT_MSM_WINDOW_MIN", 0);
    const int max_window =
        environment_log_size("HYPERCAT_MSM_WINDOW_MAX", 0);
    if ((min_window == 0) != (max_window == 0) ||
        min_window > max_window) {
        throw std::runtime_error(
            "MSM window bounds must both be zero or form an ordered range");
    }

    for (int log_size = min_log_size; log_size <= max_log_size; ++log_size) {
        const std::int64_t count = std::int64_t{1} << log_size;
        benchmark::RegisterBenchmark(
            "ScalarMsm/FrToCanonical", BM_ScalarMsm_FrToCanonical)
            ->Arg(count)
            ->UseRealTime();
        if (group == "all" || group == "g1") {
            benchmark::RegisterBenchmark(
                "ScalarMsm/G1Canonical", BM_ScalarMsm_G1Canonical)
                ->Arg(count)
                ->UseRealTime();
            benchmark::RegisterBenchmark(
                "ScalarMsm/G1NativeFr", BM_ScalarMsm_G1NativeFr)
                ->Arg(count)
                ->UseRealTime();
            benchmark::RegisterBenchmark(
                "ScalarMsm/G1SelectedFr", BM_ScalarMsm_G1SelectedFr)
                ->Arg(count)
                ->UseRealTime();
            for (int window = min_window; window <= max_window && window != 0;
                 ++window) {
                benchmark::RegisterBenchmark(
                    "ScalarMsm/G1ExplicitWindow",
                    BM_ScalarMsm_G1ExplicitWindow)
                    ->Args({count, window})
                    ->UseRealTime();
            }
        }
        if (group == "all" || group == "g2") {
            benchmark::RegisterBenchmark(
                "ScalarMsm/G2Canonical", BM_ScalarMsm_G2Canonical)
                ->Arg(count)
                ->UseRealTime();
            benchmark::RegisterBenchmark(
                "ScalarMsm/G2NativeFr", BM_ScalarMsm_G2NativeFr)
                ->Arg(count)
                ->UseRealTime();
            benchmark::RegisterBenchmark(
                "ScalarMsm/G2SelectedFr", BM_ScalarMsm_G2SelectedFr)
                ->Arg(count)
                ->UseRealTime();
            for (int window = min_window; window <= max_window && window != 0;
                 ++window) {
                benchmark::RegisterBenchmark(
                    "ScalarMsm/G2ExplicitWindow",
                    BM_ScalarMsm_G2ExplicitWindow)
                    ->Args({count, window})
                    ->UseRealTime();
            }
        }
    }
    return true;
}

const bool scalar_msm_benchmarks_registered = register_scalar_msm_benchmarks();

} // namespace

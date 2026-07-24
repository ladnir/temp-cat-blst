#include <benchmark/benchmark.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#include "bindings/blst.h"

namespace {

blst_fr batch_affine_fr_from_u64(std::uint64_t value)
{
    const std::array<std::uint64_t, 4> words = {value, 0, 0, 0};
    blst_fr out = {};
    blst_fr_from_uint64(&out, words.data());
    return out;
}

std::vector<blst_fr> batch_affine_make_scalars(std::size_t count)
{
    const blst_fr multiplier =
        batch_affine_fr_from_u64(6364136223846793005ULL);
    const blst_fr increment =
        batch_affine_fr_from_u64(1442695040888963407ULL);
    blst_fr state = batch_affine_fr_from_u64(17);
    std::vector<blst_fr> scalars(count);
    for (auto& scalar : scalars) {
        scalar = state;
        blst_fr_mul(&state, &state, &multiplier);
        blst_fr_add(&state, &state, &increment);
    }
    return scalars;
}

std::vector<blst_p1_affine> batch_affine_make_points(std::size_t count)
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

std::vector<blst_p2_affine> batch_affine_make_g2_points(std::size_t count)
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

std::vector<limb_t> batch_affine_scratch_words(std::size_t bytes)
{
    if (bytes == 0 || bytes % sizeof(limb_t) != 0) {
        throw std::runtime_error("invalid batch-affine scratch size");
    }
    return std::vector<limb_t>(bytes / sizeof(limb_t));
}

struct BatchAffineInputs
{
    BatchAffineInputs(
        std::size_t count,
        std::size_t window,
        std::size_t batch_size)
        : points(batch_affine_make_points(count)),
          scalars(batch_affine_make_scalars(count)),
          pippenger_scratch(batch_affine_scratch_words(
              blst_p1s_mult_pippenger_fr_scratch_sizeof(count))),
          batch_scratch(batch_affine_scratch_words(
              blst_p1s_mult_batch_affine_fr_scratch_sizeof(
                  count, window, batch_size)))
    {
        blst_p1 expected = {};
        blst_p1 actual = {};
        if (!blst_p1s_mult_pippenger_fr(
                &expected, points.data(), scalars.data(), count,
                pippenger_scratch.data(),
                pippenger_scratch.size() * sizeof(limb_t)) ||
            !blst_p1s_mult_batch_affine_fr(
                &actual, points.data(), scalars.data(), count,
                batch_scratch.data(),
                batch_scratch.size() * sizeof(limb_t),
                window, batch_size) ||
            !blst_p1_is_equal(&expected, &actual)) {
            throw std::runtime_error(
                "production batch-affine MSM differs from Pippenger");
        }
    }

    std::vector<blst_p1_affine> points;
    std::vector<blst_fr> scalars;
    std::vector<limb_t> pippenger_scratch;
    std::vector<limb_t> batch_scratch;
};

struct G2BatchAffineInputs
{
    G2BatchAffineInputs(
        std::size_t count,
        std::size_t window,
        std::size_t batch_size)
        : points(batch_affine_make_g2_points(count)),
          scalars(batch_affine_make_scalars(count)),
          pippenger_scratch(batch_affine_scratch_words(
              blst_p2s_mult_pippenger_fr_scratch_sizeof(count))),
          batch_scratch(batch_affine_scratch_words(
              blst_p2s_mult_batch_affine_fr_scratch_sizeof(
                  count, window, batch_size)))
    {
        blst_p2 expected = {};
        blst_p2 actual = {};
        if (!blst_p2s_mult_pippenger_fr(
                &expected, points.data(), scalars.data(), count,
                pippenger_scratch.data(),
                pippenger_scratch.size() * sizeof(limb_t)) ||
            !blst_p2s_mult_batch_affine_fr(
                &actual, points.data(), scalars.data(), count,
                batch_scratch.data(),
                batch_scratch.size() * sizeof(limb_t),
                window, batch_size) ||
            !blst_p2_is_equal(&expected, &actual)) {
            throw std::runtime_error(
                "production G2 batch-affine MSM differs from Pippenger");
        }
    }

    std::vector<blst_p2_affine> points;
    std::vector<blst_fr> scalars;
    std::vector<limb_t> pippenger_scratch;
    std::vector<limb_t> batch_scratch;
};

template<std::size_t BatchSize>
void BM_G1BatchAffineProduction(benchmark::State& state)
{
    const auto count = static_cast<std::size_t>(state.range(0));
    const auto window = static_cast<std::size_t>(state.range(1));
    BatchAffineInputs input(count, window, BatchSize);
    blst_p1 output = {};
    for (auto _ : state) {
        const bool success = blst_p1s_mult_batch_affine_fr(
            &output, input.points.data(), input.scalars.data(), count,
            input.batch_scratch.data(),
            input.batch_scratch.size() * sizeof(limb_t),
            window, BatchSize);
        benchmark::DoNotOptimize(success);
        benchmark::DoNotOptimize(output);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(
        static_cast<std::int64_t>(state.iterations() * count));
    state.counters["points"] = static_cast<double>(count);
    state.counters["window"] = static_cast<double>(window);
    state.counters["batch_size"] = static_cast<double>(BatchSize);
    state.counters["scratch_bytes"] = static_cast<double>(
        input.batch_scratch.size() * sizeof(limb_t));
}

template<std::size_t BatchSize>
void BM_G2BatchAffineProduction(benchmark::State& state)
{
    const auto count = static_cast<std::size_t>(state.range(0));
    const auto window = static_cast<std::size_t>(state.range(1));
    G2BatchAffineInputs input(count, window, BatchSize);
    blst_p2 output = {};
    for (auto _ : state) {
        const bool success = blst_p2s_mult_batch_affine_fr(
            &output, input.points.data(), input.scalars.data(), count,
            input.batch_scratch.data(),
            input.batch_scratch.size() * sizeof(limb_t),
            window, BatchSize);
        benchmark::DoNotOptimize(success);
        benchmark::DoNotOptimize(output);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(
        static_cast<std::int64_t>(state.iterations() * count));
    state.counters["points"] = static_cast<double>(count);
    state.counters["window"] = static_cast<double>(window);
    state.counters["batch_size"] = static_cast<double>(BatchSize);
    state.counters["scratch_bytes"] = static_cast<double>(
        input.batch_scratch.size() * sizeof(limb_t));
}

int batch_affine_environment_int(
    const char* name,
    int fallback,
    int minimum,
    int maximum)
{
    const char* value = std::getenv(name);
    if (value == nullptr) {
        return fallback;
    }
    const int parsed = std::stoi(value);
    if (parsed < minimum || parsed > maximum) {
        throw std::runtime_error(
            std::string(name) + " is outside its supported range");
    }
    return parsed;
}

std::size_t batch_affine_default_window(std::size_t count)
{
    if (count < 4) {
        return 1;
    }
    if (count < 32) {
        return 3;
    }
    return static_cast<std::size_t>(
        std::ceil(std::log(static_cast<double>(count))));
}

bool register_batch_affine_benchmarks()
{
    const char* enabled = std::getenv("HYPERCAT_MSM_BATCH_AFFINE");
    if (enabled == nullptr || std::string(enabled) != "1") {
        return true;
    }

    const int min_log_size = batch_affine_environment_int(
        "HYPERCAT_MSM_MIN_LOG_SIZE", 10, 1, 20);
    const int max_log_size = batch_affine_environment_int(
        "HYPERCAT_MSM_MAX_LOG_SIZE", 16, 1, 20);
    if (min_log_size > max_log_size) {
        throw std::runtime_error(
            "batch-affine minimum log size exceeds maximum");
    }
    const int explicit_min_window = batch_affine_environment_int(
        "HYPERCAT_BATCH_AFFINE_WINDOW_MIN", 0, 0, 20);
    const int explicit_max_window = batch_affine_environment_int(
        "HYPERCAT_BATCH_AFFINE_WINDOW_MAX", 0, 0, 20);
    if ((explicit_min_window == 0) != (explicit_max_window == 0) ||
        explicit_min_window > explicit_max_window) {
        throw std::runtime_error(
            "batch-affine window bounds must both be zero or ordered");
    }

    const char* batch_value =
        std::getenv("HYPERCAT_BATCH_AFFINE_BATCH");
    const std::string batch =
        batch_value == nullptr ? "all" : batch_value;
    if (batch != "all" && batch != "128" && batch != "384") {
        throw std::runtime_error(
            "batch-affine batch must be all, 128, or 384");
    }
    const char* group_value =
        std::getenv("HYPERCAT_BATCH_AFFINE_GROUP");
    const std::string group =
        group_value == nullptr ? "g1" : group_value;
    if (group != "all" && group != "g1" && group != "g2") {
        throw std::runtime_error(
            "batch-affine group must be all, g1, or g2");
    }

    for (int log_size = min_log_size;
         log_size <= max_log_size; ++log_size) {
        const std::int64_t count = std::int64_t{1} << log_size;
        const int min_window = explicit_min_window == 0
            ? static_cast<int>(batch_affine_default_window(
                  static_cast<std::size_t>(count)))
            : explicit_min_window;
        const int max_window = explicit_max_window == 0
            ? min_window
            : explicit_max_window;
        for (int window = min_window; window <= max_window; ++window) {
            if ((group == "all" || group == "g1") &&
                (batch == "all" || batch == "128")) {
                benchmark::RegisterBenchmark(
                    "ScalarMsm/G1BatchAffineProduction128",
                    BM_G1BatchAffineProduction<128>)
                    ->Args({count, window})
                    ->UseRealTime();
            }
            if ((group == "all" || group == "g1") &&
                (batch == "all" || batch == "384")) {
                benchmark::RegisterBenchmark(
                    "ScalarMsm/G1BatchAffineProduction384",
                    BM_G1BatchAffineProduction<384>)
                    ->Args({count, window})
                    ->UseRealTime();
            }
            if ((group == "all" || group == "g2") &&
                (batch == "all" || batch == "128")) {
                benchmark::RegisterBenchmark(
                    "ScalarMsm/G2BatchAffineProduction128",
                    BM_G2BatchAffineProduction<128>)
                    ->Args({count, window})
                    ->UseRealTime();
            }
            if ((group == "all" || group == "g2") &&
                (batch == "all" || batch == "384")) {
                benchmark::RegisterBenchmark(
                    "ScalarMsm/G2BatchAffineProduction384",
                    BM_G2BatchAffineProduction<384>)
                    ->Args({count, window})
                    ->UseRealTime();
            }
        }
    }
    return true;
}

const bool batch_affine_benchmarks_registered =
    register_batch_affine_benchmarks();

} // namespace

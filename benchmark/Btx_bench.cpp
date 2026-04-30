#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "BTE/btx/core.h"

using namespace bte::btx;

namespace
{
    class DeterministicRng
    {
      public:
        explicit DeterministicRng(std::uint64_t seed)
            : state_(seed != 0 ? seed : 0x9e3779b97f4a7c15ULL)
        {
        }

        std::uint64_t next_u64()
        {
            std::uint64_t x = state_;
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            state_ = x;
            return x * 2685821657736338717ULL;
        }

        void fill(const std::span<std::uint8_t> out)
        {
            std::uint64_t word = 0;
            std::size_t remaining = 0;
            for (std::uint8_t& byte : out) {
                if (remaining == 0) {
                    word = next_u64();
                    remaining = sizeof(word);
                }
                byte = static_cast<std::uint8_t>(word & 0xffU);
                word >>= 8;
                --remaining;
            }
        }

      private:
        std::uint64_t state_;
    };

    std::vector<Ciphertext> extract_ciphertexts(const std::span<const EncryptionResult> encrypted)
    {
        std::vector<Ciphertext> ciphertexts;
        ciphertexts.reserve(encrypted.size());
        for (const EncryptionResult& item : encrypted) {
            ciphertexts.push_back(item.ciphertext);
        }
        return ciphertexts;
    }

    void set_batch_counter(benchmark::State& state, const std::size_t batch_size)
    {
        state.counters["batch"] = static_cast<double>(batch_size);
        state.SetItemsProcessed(
            static_cast<std::int64_t>(state.iterations()) *
            static_cast<std::int64_t>(batch_size));
    }

    void set_threshold_counter(benchmark::State& state, const std::size_t threshold)
    {
        state.counters["threshold"] = static_cast<double>(threshold);
        state.SetItemsProcessed(
            static_cast<std::int64_t>(state.iterations()) *
            static_cast<std::int64_t>(threshold));
    }

    void set_batch_and_threshold_counters(
        benchmark::State& state,
        const std::size_t batch_size,
        const std::size_t threshold)
    {
        state.counters["batch"] = static_cast<double>(batch_size);
        state.counters["threshold"] = static_cast<double>(threshold);
        state.SetItemsProcessed(
            static_cast<std::int64_t>(state.iterations()) *
            static_cast<std::int64_t>(batch_size));
    }

    void apply_batch_args(benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int batch_size : { 1, 8, 32, 64, 128, 256, 512 }) {
            benchmark_handle->Arg(batch_size);
        }
    }

    void apply_threshold_args(benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int threshold : { 2, 4, 8, 16, 512 }) {
            benchmark_handle->Arg(threshold);
        }
    }

    struct BtxBenchmarkContext
    {
        BtxParameters params{};
        DeterministicRng rng;
        SetupResult setup;
        PreparedPublicParameters prepared_public_parameters;
        std::vector<EncryptionResult> encrypted;
        std::vector<Ciphertext> ciphertexts;
        ThresholdRequest request;
        std::vector<PartialDecryptionShare> threshold_shares;
        blst::G1Point aggregate_predecryption_key = blst::G1Point::infinity();
        BatchDecryptionPrecomputation precomputed_naive;
        BatchDecryptionPrecomputation precomputed_fft_public;
        BatchDecryptionPrecomputation precomputed_fft_prepared;
        BatchDecryptionPrecomputation precomputed_multipair;
        BatchDecryptionPrecomputation precomputed_auto;

        BtxBenchmarkContext(
            const std::size_t batch_size,
            const std::size_t threshold,
            const std::size_t num_servers,
            const std::uint64_t seed)
            : params{ std::max<std::size_t>(batch_size, 1), num_servers, threshold, 32,
                      "BTE_BTX_KDF_SALT_V1", "BTE_BTX_KDF_INFO_V1" },
              rng(seed),
              setup(Setup(params, [this](std::span<std::uint8_t> out) { rng.fill(out); })),
              prepared_public_parameters(PreparePublicParameters(setup.public_parameters))
        {
            encrypted.reserve(batch_size);
            for (std::size_t i = 0; i < batch_size; ++i) {
                encrypted.push_back(Encrypt(
                    params,
                    setup.public_parameters,
                    [this](std::span<std::uint8_t> out) { rng.fill(out); }));
            }

            ciphertexts = extract_ciphertexts(
                std::span<const EncryptionResult>(encrypted.data(), encrypted.size()));
            request = BuildThresholdRequest(setup.public_parameters, ciphertexts);

            threshold_shares.reserve(threshold);
            for (std::size_t i = 0; i < threshold; ++i) {
                threshold_shares.push_back(CreatePartialDecryptionShare(
                    setup.secret_shares[i],
                    request));
            }

            aggregate_predecryption_key = CombinePartialDecryptionShares(
                std::span<const PartialDecryptionShare>(
                    threshold_shares.data(),
                    threshold_shares.size()),
                setup.public_parameters);

            precomputed_naive = PrecomputeBatchDecryption(
                setup.public_parameters,
                ciphertexts,
                CrossTermAlgorithm::Naive);
            precomputed_fft_public = PrecomputeBatchDecryption(
                setup.public_parameters,
                ciphertexts,
                CrossTermAlgorithm::MiddleProductFft);
            precomputed_fft_prepared = PrecomputeBatchDecryption(
                prepared_public_parameters,
                ciphertexts,
                CrossTermAlgorithm::MiddleProductFft);
            precomputed_multipair = PrecomputeBatchDecryption(
                setup.public_parameters,
                ciphertexts,
                CrossTermAlgorithm::MiddleProductMultipair);
            precomputed_auto = PrecomputeBatchDecryption(
                setup.public_parameters,
                ciphertexts,
                CrossTermAlgorithm::MiddleProductAuto);

            // Warm the lazily materialized alpha cache so prepared finalize
            // benchmarks measure the steady-state online path.
            benchmark::DoNotOptimize(FinalizeBatchDecryption(
                params,
                prepared_public_parameters,
                precomputed_fft_prepared,
                aggregate_predecryption_key));
        }
    };

    constexpr std::size_t kDefaultBatchSize = 32;
    constexpr std::size_t kDefaultThreshold = 8;
    constexpr std::size_t kDefaultServerCount = 16;

    static void BM_Btx_Setup(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        const BtxParameters params{
            batch_size,
            kDefaultServerCount,
            kDefaultThreshold,
            32,
            "BTE_BTX_KDF_SALT_V1",
            "BTE_BTX_KDF_INFO_V1"
        };
        DeterministicRng rng(0xb7000000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(Setup(
                params,
                [&rng](std::span<std::uint8_t> out) { rng.fill(out); }));
            benchmark::ClobberMemory();
        }

        set_batch_and_threshold_counters(state, batch_size, kDefaultThreshold);
    }

    static void BM_Btx_PreparePublicParameters(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb7010000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(PreparePublicParameters(ctx.setup.public_parameters));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_Encrypt(benchmark::State& state)
    {
        BtxBenchmarkContext ctx(kDefaultBatchSize, kDefaultThreshold, kDefaultServerCount, 0xb7020000ULL);
        DeterministicRng rng(0xb7021000ULL);

        for (auto _ : state) {
            benchmark::DoNotOptimize(Encrypt(
                ctx.params,
                ctx.setup.public_parameters,
                [&rng](std::span<std::uint8_t> out) { rng.fill(out); }));
            benchmark::ClobberMemory();
        }

        state.counters["batch"] = 1.0;
        state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()));
    }

    static void BM_Btx_BuildThresholdRequest(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb7030000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(BuildThresholdRequest(ctx.setup.public_parameters, ctx.ciphertexts));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_CreatePartialDecryptionShare(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb7040000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(CreatePartialDecryptionShare(
                ctx.setup.secret_shares.front(),
                ctx.request));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_CombinePartialDecryptionShares(benchmark::State& state)
    {
        const std::size_t threshold = static_cast<std::size_t>(state.range(0));
        const std::size_t server_count = std::max(kDefaultServerCount, threshold);
        BtxBenchmarkContext ctx(kDefaultBatchSize, threshold, server_count, 0xb7050000ULL + threshold);

        for (auto _ : state) {
            benchmark::DoNotOptimize(CombinePartialDecryptionShares(
                std::span<const PartialDecryptionShare>(
                    ctx.threshold_shares.data(),
                    ctx.threshold_shares.size()),
                ctx.setup.public_parameters));
            benchmark::ClobberMemory();
        }

        set_threshold_counter(state, threshold);
    }

    static void BM_Btx_ComputeCrossTerms_Naive(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb7060000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(ComputeCrossTerms(
                ctx.setup.public_parameters,
                ctx.ciphertexts,
                CrossTermAlgorithm::Naive));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_ComputeCrossTerms_FftPublic(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb7070000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(ComputeCrossTerms(
                ctx.setup.public_parameters,
                ctx.ciphertexts,
                CrossTermAlgorithm::MiddleProductFft));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_ComputeCrossTerms_FftPrepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb7080000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(ComputeCrossTerms(
                ctx.prepared_public_parameters,
                ctx.ciphertexts,
                CrossTermAlgorithm::MiddleProductFft));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_ComputeCrossTerms_Multipair(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb7090000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(ComputeCrossTerms(
                ctx.setup.public_parameters,
                ctx.ciphertexts,
                CrossTermAlgorithm::MiddleProductMultipair));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_ComputeCrossTerms_Auto(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb70a0000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(ComputeCrossTerms(
                ctx.setup.public_parameters,
                ctx.ciphertexts,
                CrossTermAlgorithm::MiddleProductAuto));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_FinalizeBatchDecryption_DirectAlpha(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb70b0000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(FinalizeBatchDecryption(
                ctx.params,
                ctx.setup.public_parameters,
                ctx.precomputed_fft_prepared,
                ctx.aggregate_predecryption_key));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_FinalizeBatchDecryption_PreparedAlpha(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb70c0000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(FinalizeBatchDecryption(
                ctx.params,
                ctx.prepared_public_parameters,
                ctx.precomputed_fft_prepared,
                ctx.aggregate_predecryption_key));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_RecoverBatchKeys_Naive(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb70d0000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(RecoverBatchKeys(
                ctx.params,
                ctx.setup.public_parameters,
                ctx.ciphertexts,
                ctx.aggregate_predecryption_key,
                CrossTermAlgorithm::Naive));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_RecoverBatchKeys_FftPublic(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb70e0000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(RecoverBatchKeys(
                ctx.params,
                ctx.setup.public_parameters,
                ctx.ciphertexts,
                ctx.aggregate_predecryption_key,
                CrossTermAlgorithm::MiddleProductFft));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_RecoverBatchKeys_FftPrepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb70f0000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(RecoverBatchKeys(
                ctx.params,
                ctx.prepared_public_parameters,
                ctx.ciphertexts,
                ctx.aggregate_predecryption_key,
                CrossTermAlgorithm::MiddleProductFft));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_RecoverBatchKeys_Multipair(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb7100000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(RecoverBatchKeys(
                ctx.params,
                ctx.setup.public_parameters,
                ctx.ciphertexts,
                ctx.aggregate_predecryption_key,
                CrossTermAlgorithm::MiddleProductMultipair));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_RecoverBatchKeys_Auto(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xb7110000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(RecoverBatchKeys(
                ctx.params,
                ctx.setup.public_parameters,
                ctx.ciphertexts,
                ctx.aggregate_predecryption_key,
                CrossTermAlgorithm::MiddleProductAuto));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }
}

BENCHMARK(BM_Btx_Setup)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_PreparePublicParameters)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_Encrypt);
BENCHMARK(BM_Btx_BuildThresholdRequest)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_CreatePartialDecryptionShare)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_CombinePartialDecryptionShares)->Apply(apply_threshold_args);
BENCHMARK(BM_Btx_ComputeCrossTerms_Naive)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_ComputeCrossTerms_FftPublic)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_ComputeCrossTerms_FftPrepared)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_ComputeCrossTerms_Multipair)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_ComputeCrossTerms_Auto)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_FinalizeBatchDecryption_DirectAlpha)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_FinalizeBatchDecryption_PreparedAlpha)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_RecoverBatchKeys_Naive)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_RecoverBatchKeys_FftPublic)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_RecoverBatchKeys_FftPrepared)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_RecoverBatchKeys_Multipair)->Apply(apply_batch_args);
BENCHMARK(BM_Btx_RecoverBatchKeys_Auto)->Apply(apply_batch_args);

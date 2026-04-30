#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "BTE/pfe/core.h"

using namespace bte::pfe;

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
        for (const int batch_size : { 1, 8, 32, 64, 128, 512 }) {
            benchmark_handle->Arg(batch_size);
        }
    }

    void apply_threshold_args(benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int threshold : { 2, 4, 8, 16, 512 }) {
            benchmark_handle->Arg(threshold);
        }
    }

    struct PfeBenchmarkContext
    {
        PfeParameters params{};
        DeterministicRng rng;
        SetupResult setup;
        PreparedPublicParameters prepared_public_parameters;
        std::vector<EncryptionResult> encrypted;
        std::vector<Ciphertext> ciphertexts;
        ThresholdRequest request;
        std::vector<PartialDecryptionShare> threshold_shares;
        blst::G1Point aggregate_predecryption_key = blst::G1Point::infinity();
        BatchDecryptionPrecomputation precomputed_public;
        BatchDecryptionPrecomputation precomputed_prepared;

        PfeBenchmarkContext(
            const std::size_t batch_size,
            const std::size_t threshold,
            const std::size_t num_servers,
            const std::uint64_t seed)
            : params{ std::max<std::size_t>(batch_size, 1), num_servers, threshold, 32,
                      "BTE_PFE_KDF_SALT_V1", "BTE_PFE_KDF_INFO_V1" },
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

            precomputed_public = PrecomputeBatchDecryption(setup.public_parameters, ciphertexts);
            precomputed_prepared = PrecomputeBatchDecryption(prepared_public_parameters, ciphertexts);
        }
    };

    constexpr std::size_t kDefaultBatchSize = 32;
    constexpr std::size_t kDefaultThreshold = 8;
    constexpr std::size_t kDefaultServerCount = 16;

    static void BM_Pfe_Setup(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        const PfeParameters params{
            batch_size,
            kDefaultServerCount,
            kDefaultThreshold,
            32,
            "BTE_PFE_KDF_SALT_V1",
            "BTE_PFE_KDF_INFO_V1"
        };
        DeterministicRng rng(0xcf500000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(Setup(
                params,
                [&rng](std::span<std::uint8_t> out) { rng.fill(out); }));
            benchmark::ClobberMemory();
        }

        set_batch_and_threshold_counters(state, batch_size, kDefaultThreshold);
    }

    static void BM_Pfe_PreparePublicParameters(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PfeBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xcf510000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(PreparePublicParameters(ctx.setup.public_parameters));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Pfe_Encrypt(benchmark::State& state)
    {
        PfeBenchmarkContext ctx(kDefaultBatchSize, kDefaultThreshold, kDefaultServerCount, 0xcf520000ULL);
        DeterministicRng rng(0xcf520100ULL);

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

    static void BM_Pfe_BuildThresholdRequest(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PfeBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xcf530000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(BuildThresholdRequest(ctx.setup.public_parameters, ctx.ciphertexts));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Pfe_CreatePartialDecryptionShare(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PfeBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xcf540000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(CreatePartialDecryptionShare(
                ctx.setup.secret_shares.front(),
                ctx.request));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Pfe_CombinePartialDecryptionShares(benchmark::State& state)
    {
        const std::size_t threshold = static_cast<std::size_t>(state.range(0));
        const std::size_t server_count = std::max(kDefaultServerCount, threshold);
        PfeBenchmarkContext ctx(kDefaultBatchSize, threshold, server_count, 0xcf550000ULL + threshold);

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

    static void BM_Pfe_PrecomputeBatchDecryption_Public(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PfeBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xcf560000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(PrecomputeBatchDecryption(
                ctx.setup.public_parameters,
                ctx.ciphertexts));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Pfe_PrecomputeBatchDecryption_Prepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PfeBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xcf570000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(PrecomputeBatchDecryption(
                ctx.prepared_public_parameters,
                ctx.ciphertexts));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Pfe_FinalizeBatchDecryption_Public(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PfeBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xcf580000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(FinalizeBatchDecryption(
                ctx.params,
                ctx.setup.public_parameters,
                ctx.precomputed_public,
                ctx.aggregate_predecryption_key));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Pfe_FinalizeBatchDecryption_Prepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PfeBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xcf590000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(FinalizeBatchDecryption(
                ctx.params,
                ctx.prepared_public_parameters,
                ctx.precomputed_prepared,
                ctx.aggregate_predecryption_key));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Pfe_RecoverBatchKeys_Public(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PfeBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xcf5a0000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(RecoverBatchKeys(
                ctx.params,
                ctx.setup.public_parameters,
                ctx.ciphertexts,
                ctx.aggregate_predecryption_key));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Pfe_RecoverBatchKeys_Prepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PfeBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xcf5b0000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(RecoverBatchKeys(
                ctx.params,
                ctx.prepared_public_parameters,
                ctx.ciphertexts,
                ctx.aggregate_predecryption_key));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }
}

BENCHMARK(BM_Pfe_Setup)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Pfe_PreparePublicParameters)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Pfe_Encrypt)->UseRealTime();
BENCHMARK(BM_Pfe_BuildThresholdRequest)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Pfe_CreatePartialDecryptionShare)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Pfe_CombinePartialDecryptionShares)->Apply(apply_threshold_args)->UseRealTime();
BENCHMARK(BM_Pfe_PrecomputeBatchDecryption_Public)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Pfe_PrecomputeBatchDecryption_Prepared)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Pfe_FinalizeBatchDecryption_Public)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Pfe_FinalizeBatchDecryption_Prepared)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Pfe_RecoverBatchKeys_Public)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Pfe_RecoverBatchKeys_Prepared)->Apply(apply_batch_args)->UseRealTime();

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <vector>

#include "BTE/beat/core.h"
#include "BTE/btx/core.h"
#include "BTE/pfe/core.h"

namespace
{
    namespace beat = bte::beat;
    namespace pfe = bte::pfe;
    namespace btx = bte::btx;

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

    std::vector<pfe::Ciphertext> extract_pfe_ciphertexts(
        const std::span<const pfe::EncryptionResult> encrypted)
    {
        std::vector<pfe::Ciphertext> ciphertexts;
        ciphertexts.reserve(encrypted.size());
        for (const pfe::EncryptionResult& item : encrypted) {
            ciphertexts.push_back(item.ciphertext);
        }
        return ciphertexts;
    }

    std::vector<btx::Ciphertext> extract_btx_ciphertexts(
        const std::span<const btx::EncryptionResult> encrypted)
    {
        std::vector<btx::Ciphertext> ciphertexts;
        ciphertexts.reserve(encrypted.size());
        for (const btx::EncryptionResult& item : encrypted) {
            ciphertexts.push_back(item.ciphertext);
        }
        return ciphertexts;
    }

    std::vector<blst::G1Point> extract_btx_c1_points(const std::span<const btx::Ciphertext> ciphertexts)
    {
        std::vector<blst::G1Point> points;
        points.reserve(ciphertexts.size());
        for (const btx::Ciphertext& ciphertext : ciphertexts) {
            points.push_back(ciphertext.c1);
        }
        return points;
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

    void apply_batch_args(benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int batch_size : { 1, 8, 32, 64, 128, 256, 512 }) {
            benchmark_handle->Arg(batch_size);
        }
    }

    void apply_threshold_args(benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int threshold : { 2, 4, 8, 16 }) {
            benchmark_handle->Arg(threshold);
        }
    }

    struct ComparisonContext;

    std::unique_ptr<ComparisonContext> build_context(
        benchmark::State& state,
        const std::size_t batch_size,
        const std::size_t threshold,
        const std::size_t num_servers,
        const std::uint64_t seed);

    struct ComparisonContext
    {
        std::size_t batch_size = 0;

        beat::BeatParameters beat_params{};
        DeterministicRng beat_rng;
        beat::SetupResult beat_setup;
        beat::PreparedPublicParameters beat_prepared;
        std::vector<beat::EncryptionResult> beat_encrypted;
        std::vector<beat::Ciphertext> beat_ciphertexts;
        beat::ThresholdRequest beat_request;
        std::vector<beat::PartialDecryptionShare> beat_shares;
        blst::G1Point beat_aggregate_predecryption_key = blst::G1Point::infinity();
        beat::BatchDecryptionPrecomputation beat_precomputed;

        pfe::PfeParameters pfe_params{};
        DeterministicRng pfe_rng;
        pfe::SetupResult pfe_setup;
        pfe::PreparedPublicParameters pfe_prepared;
        std::vector<pfe::EncryptionResult> pfe_encrypted;
        std::vector<pfe::Ciphertext> pfe_ciphertexts;
        pfe::ThresholdRequest pfe_request;
        std::vector<pfe::PartialDecryptionShare> pfe_shares;
        blst::G1Point pfe_aggregate_predecryption_key = blst::G1Point::infinity();
        pfe::BatchDecryptionPrecomputation pfe_precomputed;

        btx::BtxParameters btx_params{};
        DeterministicRng btx_rng;
        btx::SetupResult btx_setup;
        btx::PreparedPublicParameters btx_prepared;
        std::vector<btx::EncryptionResult> btx_encrypted;
        std::vector<btx::Ciphertext> btx_ciphertexts;
        btx::ThresholdRequest btx_request;
        std::vector<btx::PartialDecryptionShare> btx_shares;
        blst::G1Point btx_aggregate_predecryption_key = blst::G1Point::infinity();
        btx::BatchDecryptionPrecomputation btx_precomputed;

        ComparisonContext(
            const std::size_t in_batch_size,
            const std::size_t threshold,
            const std::size_t num_servers,
            const std::uint64_t seed)
            : batch_size(std::max<std::size_t>(in_batch_size, 1)),
              beat_params{ batch_size, num_servers, threshold, 32,
                           "BTE_BEAT_KDF_SALT_V1", "BTE_BEAT_KDF_INFO_V1" },
              beat_rng(seed ^ 0x0101010101010101ULL),
              beat_setup(beat::Setup(
                  beat_params,
                  [this](std::span<std::uint8_t> out) { beat_rng.fill(out); })),
              beat_prepared(beat::PreparePublicParameters(beat_setup.public_parameters)),
              pfe_params{ batch_size, num_servers, threshold, 32,
                          "BTE_PFE_KDF_SALT_V1", "BTE_PFE_KDF_INFO_V1" },
              pfe_rng(seed ^ 0x1111111111111111ULL),
              pfe_setup(pfe::Setup(
                  pfe_params,
                  [this](std::span<std::uint8_t> out) { pfe_rng.fill(out); })),
              pfe_prepared(pfe::PreparePublicParameters(pfe_setup.public_parameters)),
              btx_params{ batch_size, num_servers, threshold, 32,
                          "BTE_BTX_KDF_SALT_V1", "BTE_BTX_KDF_INFO_V1" },
              btx_rng(seed ^ 0x2222222222222222ULL),
              btx_setup(btx::Setup(
                  btx_params,
                  [this](std::span<std::uint8_t> out) { btx_rng.fill(out); })),
              btx_prepared(btx::PreparePublicParameters(btx_setup.public_parameters))
        {
            beat_encrypted.reserve(batch_size);
            pfe_encrypted.reserve(batch_size);
            btx_encrypted.reserve(batch_size);
            for (std::size_t i = 0; i < batch_size; ++i) {
                beat_encrypted.push_back(beat::Encrypt(
                    beat_params,
                    beat_setup.public_parameters,
                    [this](std::span<std::uint8_t> out) { beat_rng.fill(out); }));
                pfe_encrypted.push_back(pfe::Encrypt(
                    pfe_params,
                    pfe_setup.public_parameters,
                    [this](std::span<std::uint8_t> out) { pfe_rng.fill(out); }));
                btx_encrypted.push_back(btx::Encrypt(
                    btx_params,
                    btx_setup.public_parameters,
                    [this](std::span<std::uint8_t> out) { btx_rng.fill(out); }));
            }

            beat_ciphertexts.reserve(batch_size);
            for (const beat::EncryptionResult& item : beat_encrypted) {
                beat_ciphertexts.push_back(item.ciphertext);
            }
            pfe_ciphertexts = extract_pfe_ciphertexts(
                std::span<const pfe::EncryptionResult>(pfe_encrypted.data(), pfe_encrypted.size()));
            btx_ciphertexts = extract_btx_ciphertexts(
                std::span<const btx::EncryptionResult>(btx_encrypted.data(), btx_encrypted.size()));

            beat_request = beat::BuildThresholdRequest(beat_setup.public_parameters, beat_ciphertexts);
            pfe_request = pfe::BuildThresholdRequest(pfe_setup.public_parameters, pfe_ciphertexts);
            btx_request = btx::BuildThresholdRequest(btx_setup.public_parameters, btx_ciphertexts);

            beat_shares.reserve(threshold);
            pfe_shares.reserve(threshold);
            btx_shares.reserve(threshold);
            for (std::size_t i = 0; i < threshold; ++i) {
                beat_shares.push_back(beat::CreatePartialDecryptionShare(
                    beat_setup.secret_shares[i],
                    beat_request));
                pfe_shares.push_back(pfe::CreatePartialDecryptionShare(
                    pfe_setup.secret_shares[i],
                    pfe_request));
                btx_shares.push_back(btx::CreatePartialDecryptionShare(
                    btx_setup.secret_shares[i],
                    btx_request));
            }

            beat_aggregate_predecryption_key = beat::CombinePartialDecryptionShares(
                std::span<const beat::PartialDecryptionShare>(beat_shares.data(), beat_shares.size()),
                beat_setup.public_parameters);
            pfe_aggregate_predecryption_key = pfe::CombinePartialDecryptionShares(
                std::span<const pfe::PartialDecryptionShare>(pfe_shares.data(), pfe_shares.size()),
                pfe_setup.public_parameters);
            btx_aggregate_predecryption_key = btx::CombinePartialDecryptionShares(
                std::span<const btx::PartialDecryptionShare>(btx_shares.data(), btx_shares.size()),
                btx_setup.public_parameters);

            beat_precomputed = beat::PrecomputeBatchDecryption(beat_prepared, beat_ciphertexts);
            pfe_precomputed = pfe::PrecomputeBatchDecryption(pfe_prepared, pfe_ciphertexts);
            btx_precomputed = btx::PrecomputeBatchDecryption(
                btx_prepared,
                btx_ciphertexts,
                btx::CrossTermAlgorithm::MiddleProductFft);

            benchmark::DoNotOptimize(beat::FinalizeBatchDecryption(
                beat_params,
                beat_prepared,
                beat_precomputed,
                beat_aggregate_predecryption_key));

            // Warm BTX's lazily populated alpha cache so the timed open path is
            // only the actual finalize step.
            benchmark::DoNotOptimize(btx::FinalizeBatchDecryption(
                btx_params,
                btx_prepared,
                btx_precomputed,
                btx_aggregate_predecryption_key));
        }
    };

    constexpr std::size_t kDefaultBatchSize = 32;
    constexpr std::size_t kDefaultThreshold = 8;
    constexpr std::size_t kDefaultServerCount = 16;

    std::unique_ptr<ComparisonContext> build_context(
        benchmark::State& state,
        const std::size_t batch_size,
        const std::size_t threshold,
        const std::size_t num_servers,
        const std::uint64_t seed)
    {
        try {
            return std::make_unique<ComparisonContext>(
                batch_size,
                threshold,
                num_servers,
                seed);
        } catch (const std::exception& ex) {
            std::fprintf(stderr, "three-way benchmark context error: %s\n", ex.what());
            state.SkipWithError("three-way benchmark context construction failed");
        } catch (...) {
            std::fprintf(stderr, "three-way benchmark context error: unknown exception\n");
            state.SkipWithError("unknown exception while constructing three-way benchmark context");
        }
        return nullptr;
    }

    static void BM_BtxPfeBeat_PartialDecrypt_Pfe(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        const std::unique_ptr<ComparisonContext> ctx =
            build_context(state, batch_size, kDefaultThreshold, kDefaultServerCount, 0xf7050000ULL + batch_size);
        if (!ctx) {
            return;
        }

        for (auto _ : state) {
            benchmark::DoNotOptimize(pfe::CreatePartialDecryptionShare(
                ctx->pfe_setup.secret_shares.front(),
                ctx->pfe_request));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_BtxPfeBeat_PartialDecrypt_Btx(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        const std::unique_ptr<ComparisonContext> ctx =
            build_context(state, batch_size, kDefaultThreshold, kDefaultServerCount, 0xf7060000ULL + batch_size);
        if (!ctx) {
            return;
        }

        for (auto _ : state) {
            benchmark::DoNotOptimize(btx::CreatePartialDecryptionShare(
                ctx->btx_setup.secret_shares.front(),
                ctx->btx_request));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_BtxPfeBeat_PartialDecrypt_Beat(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        const std::unique_ptr<ComparisonContext> ctx =
            build_context(state, batch_size, kDefaultThreshold, kDefaultServerCount, 0xf7070000ULL + batch_size);
        if (!ctx) {
            return;
        }

        for (auto _ : state) {
            benchmark::DoNotOptimize(beat::CreatePartialDecryptionShare(
                ctx->beat_setup.secret_shares.front(),
                ctx->beat_request));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_BtxPfeBeat_Combine_Pfe(benchmark::State& state)
    {
        const std::size_t threshold = static_cast<std::size_t>(state.range(0));
        const std::size_t server_count = std::max(kDefaultServerCount, threshold);
        const std::unique_ptr<ComparisonContext> ctx =
            build_context(state, kDefaultBatchSize, threshold, server_count, 0xf7080000ULL + threshold);
        if (!ctx) {
            return;
        }

        for (auto _ : state) {
            benchmark::DoNotOptimize(pfe::CombinePartialDecryptionShares(
                std::span<const pfe::PartialDecryptionShare>(ctx->pfe_shares.data(), threshold),
                ctx->pfe_setup.public_parameters));
            benchmark::ClobberMemory();
        }

        set_threshold_counter(state, threshold);
    }

    static void BM_BtxPfeBeat_Combine_Btx(benchmark::State& state)
    {
        const std::size_t threshold = static_cast<std::size_t>(state.range(0));
        const std::size_t server_count = std::max(kDefaultServerCount, threshold);
        const std::unique_ptr<ComparisonContext> ctx =
            build_context(state, kDefaultBatchSize, threshold, server_count, 0xf7090000ULL + threshold);
        if (!ctx) {
            return;
        }

        for (auto _ : state) {
            benchmark::DoNotOptimize(btx::CombinePartialDecryptionShares(
                std::span<const btx::PartialDecryptionShare>(ctx->btx_shares.data(), threshold),
                ctx->btx_setup.public_parameters));
            benchmark::ClobberMemory();
        }

        set_threshold_counter(state, threshold);
    }

    static void BM_BtxPfeBeat_Combine_Beat(benchmark::State& state)
    {
        const std::size_t threshold = static_cast<std::size_t>(state.range(0));
        const std::size_t server_count = std::max(kDefaultServerCount, threshold);
        const std::unique_ptr<ComparisonContext> ctx =
            build_context(state, kDefaultBatchSize, threshold, server_count, 0xf70a0000ULL + threshold);
        if (!ctx) {
            return;
        }

        for (auto _ : state) {
            benchmark::DoNotOptimize(beat::CombinePartialDecryptionShares(
                std::span<const beat::PartialDecryptionShare>(ctx->beat_shares.data(), threshold),
                ctx->beat_setup.public_parameters));
            benchmark::ClobberMemory();
        }

        set_threshold_counter(state, threshold);
    }

    static void BM_BtxPfeBeat_Precompute_Pfe(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        const std::unique_ptr<ComparisonContext> ctx =
            build_context(state, batch_size, kDefaultThreshold, kDefaultServerCount, 0xf7100000ULL + batch_size);
        if (!ctx) {
            return;
        }

        for (auto _ : state) {
            benchmark::DoNotOptimize(pfe::PrecomputeBatchDecryption(
                ctx->pfe_prepared,
                ctx->pfe_ciphertexts));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_BtxPfeBeat_Precompute_Btx(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        const std::unique_ptr<ComparisonContext> ctx =
            build_context(state, batch_size, kDefaultThreshold, kDefaultServerCount, 0xf7200000ULL + batch_size);
        if (!ctx) {
            return;
        }

        for (auto _ : state) {
            benchmark::DoNotOptimize(btx::PrecomputeBatchDecryption(
                ctx->btx_prepared,
                ctx->btx_ciphertexts,
                btx::CrossTermAlgorithm::MiddleProductFft));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_BtxPfeBeat_Precompute_Beat(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        const std::unique_ptr<ComparisonContext> ctx =
            build_context(state, batch_size, kDefaultThreshold, kDefaultServerCount, 0xf7300000ULL + batch_size);
        if (!ctx) {
            return;
        }

        for (auto _ : state) {
            benchmark::DoNotOptimize(beat::PrecomputeBatchDecryption(
                ctx->beat_prepared,
                ctx->beat_ciphertexts));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_BtxPfeBeat_Open_Pfe(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        const std::unique_ptr<ComparisonContext> ctx =
            build_context(state, batch_size, kDefaultThreshold, kDefaultServerCount, 0xf7400000ULL + batch_size);
        if (!ctx) {
            return;
        }

        for (auto _ : state) {
            benchmark::DoNotOptimize(pfe::FinalizeBatchDecryption(
                ctx->pfe_params,
                ctx->pfe_prepared,
                ctx->pfe_precomputed,
                ctx->pfe_aggregate_predecryption_key));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_BtxPfeBeat_Open_Btx(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        const std::unique_ptr<ComparisonContext> ctx =
            build_context(state, batch_size, kDefaultThreshold, kDefaultServerCount, 0xf7500000ULL + batch_size);
        if (!ctx) {
            return;
        }

        for (auto _ : state) {
            benchmark::DoNotOptimize(btx::FinalizeBatchDecryption(
                ctx->btx_params,
                ctx->btx_prepared,
                ctx->btx_precomputed,
                ctx->btx_aggregate_predecryption_key));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_BtxPfeBeat_Open_Beat(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        const std::unique_ptr<ComparisonContext> ctx =
            build_context(state, batch_size, kDefaultThreshold, kDefaultServerCount, 0xf7600000ULL + batch_size);
        if (!ctx) {
            return;
        }

        for (auto _ : state) {
            benchmark::DoNotOptimize(beat::FinalizeBatchDecryption(
                ctx->beat_params,
                ctx->beat_prepared,
                ctx->beat_precomputed,
                ctx->beat_aggregate_predecryption_key));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }
}

BENCHMARK(BM_BtxPfeBeat_PartialDecrypt_Pfe)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_BtxPfeBeat_PartialDecrypt_Btx)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_BtxPfeBeat_PartialDecrypt_Beat)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_BtxPfeBeat_Combine_Pfe)->Apply(apply_threshold_args)->UseRealTime();
BENCHMARK(BM_BtxPfeBeat_Combine_Btx)->Apply(apply_threshold_args)->UseRealTime();
BENCHMARK(BM_BtxPfeBeat_Combine_Beat)->Apply(apply_threshold_args)->UseRealTime();
BENCHMARK(BM_BtxPfeBeat_Precompute_Pfe)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_BtxPfeBeat_Precompute_Btx)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_BtxPfeBeat_Precompute_Beat)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_BtxPfeBeat_Open_Pfe)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_BtxPfeBeat_Open_Btx)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_BtxPfeBeat_Open_Beat)->Apply(apply_batch_args)->UseRealTime();

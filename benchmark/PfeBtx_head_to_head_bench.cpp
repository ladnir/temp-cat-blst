#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstdint>
#include <span>
#include <vector>

#include "BTE/btx/core.h"
#include "BTE/pfe/core.h"

namespace
{
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

    void set_batch_counter(benchmark::State& state, const std::size_t batch_size)
    {
        state.counters["batch"] = static_cast<double>(batch_size);
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

    struct ComparisonContext
    {
        std::size_t batch_size = 0;

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
            pfe_encrypted.reserve(batch_size);
            btx_encrypted.reserve(batch_size);
            for (std::size_t i = 0; i < batch_size; ++i) {
                pfe_encrypted.push_back(pfe::Encrypt(
                    pfe_params,
                    pfe_setup.public_parameters,
                    [this](std::span<std::uint8_t> out) { pfe_rng.fill(out); }));
                btx_encrypted.push_back(btx::Encrypt(
                    btx_params,
                    btx_setup.public_parameters,
                    [this](std::span<std::uint8_t> out) { btx_rng.fill(out); }));
            }

            pfe_ciphertexts = extract_pfe_ciphertexts(
                std::span<const pfe::EncryptionResult>(pfe_encrypted.data(), pfe_encrypted.size()));
            btx_ciphertexts = extract_btx_ciphertexts(
                std::span<const btx::EncryptionResult>(btx_encrypted.data(), btx_encrypted.size()));

            pfe_request = pfe::BuildThresholdRequest(pfe_setup.public_parameters, pfe_ciphertexts);
            btx_request = btx::BuildThresholdRequest(btx_setup.public_parameters, btx_ciphertexts);

            pfe_shares.reserve(threshold);
            btx_shares.reserve(threshold);
            for (std::size_t i = 0; i < threshold; ++i) {
                pfe_shares.push_back(pfe::CreatePartialDecryptionShare(
                    pfe_setup.secret_shares[i],
                    pfe_request));
                btx_shares.push_back(btx::CreatePartialDecryptionShare(
                    btx_setup.secret_shares[i],
                    btx_request));
            }

            pfe_aggregate_predecryption_key = pfe::CombinePartialDecryptionShares(
                std::span<const pfe::PartialDecryptionShare>(pfe_shares.data(), pfe_shares.size()),
                pfe_setup.public_parameters);
            btx_aggregate_predecryption_key = btx::CombinePartialDecryptionShares(
                std::span<const btx::PartialDecryptionShare>(btx_shares.data(), btx_shares.size()),
                btx_setup.public_parameters);

            pfe_precomputed = pfe::PrecomputeBatchDecryption(pfe_prepared, pfe_ciphertexts);
            btx_precomputed = btx::PrecomputeBatchDecryption(
                btx_prepared,
                btx_ciphertexts,
                btx::CrossTermAlgorithm::MiddleProductFft);

            // Warm BTX's lazily materialized alpha cache so the timed "online"
            // benchmarks measure only the prepared finalize path.
            benchmark::DoNotOptimize(btx::FinalizeBatchDecryption(
                btx_params,
                btx_prepared,
                btx_precomputed,
                btx_aggregate_predecryption_key));
        }
    };

    constexpr std::size_t kDefaultThreshold = 8;
    constexpr std::size_t kDefaultServerCount = 16;

    static void BM_PfeVsBtx_Offline_PfePrepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        ComparisonContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xd7100000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(pfe::PrecomputeBatchDecryption(
                ctx.pfe_prepared,
                ctx.pfe_ciphertexts));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_PfeVsBtx_Offline_BtxPrepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        ComparisonContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xd7200000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(btx::PrecomputeBatchDecryption(
                ctx.btx_prepared,
                ctx.btx_ciphertexts,
                btx::CrossTermAlgorithm::MiddleProductFft));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_PfeVsBtx_Online_PfePrepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        ComparisonContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xd7300000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(pfe::FinalizeBatchDecryption(
                ctx.pfe_params,
                ctx.pfe_prepared,
                ctx.pfe_precomputed,
                ctx.pfe_aggregate_predecryption_key));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_PfeVsBtx_Online_BtxPrepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        ComparisonContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xd7400000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(btx::FinalizeBatchDecryption(
                ctx.btx_params,
                ctx.btx_prepared,
                ctx.btx_precomputed,
                ctx.btx_aggregate_predecryption_key));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_PfeVsBtx_EndToEnd_PfePrepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        ComparisonContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xd7500000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(pfe::RecoverBatchKeys(
                ctx.pfe_params,
                ctx.pfe_prepared,
                ctx.pfe_ciphertexts,
                ctx.pfe_aggregate_predecryption_key));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_PfeVsBtx_EndToEnd_BtxPrepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        ComparisonContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xd7600000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(btx::RecoverBatchKeys(
                ctx.btx_params,
                ctx.btx_prepared,
                ctx.btx_ciphertexts,
                ctx.btx_aggregate_predecryption_key,
                btx::CrossTermAlgorithm::MiddleProductFft));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }
}

BENCHMARK(BM_PfeVsBtx_Offline_PfePrepared)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_PfeVsBtx_Offline_BtxPrepared)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_PfeVsBtx_Online_PfePrepared)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_PfeVsBtx_Online_BtxPrepared)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_PfeVsBtx_EndToEnd_PfePrepared)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_PfeVsBtx_EndToEnd_BtxPrepared)->Apply(apply_batch_args)->UseRealTime();

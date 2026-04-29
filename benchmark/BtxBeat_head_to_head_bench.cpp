#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstdint>
#include <span>
#include <vector>

#include "BTE/btx/core.h"
#include "blst/fft/middle_product.h"

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

    std::vector<blst::G1Point> extract_c1_points(const std::span<const Ciphertext> ciphertexts)
    {
        std::vector<blst::G1Point> points;
        points.reserve(ciphertexts.size());
        for (const Ciphertext& ciphertext : ciphertexts) {
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

    void apply_batch_args(benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int batch_size : { 1, 8, 32, 64, 128, 256, 512 }) {
            benchmark_handle->Arg(batch_size);
        }
    }

    struct HeadToHeadContext
    {
        BtxParameters params{};
        DeterministicRng rng;
        SetupResult setup;
        PreparedPublicParameters prepared_public_parameters;
        std::vector<EncryptionResult> encrypted;
        std::vector<Ciphertext> ciphertexts;
        std::vector<blst::G1Point> c1_points;
        ThresholdRequest request;
        std::vector<PartialDecryptionShare> threshold_shares;
        blst::G1Point aggregate_predecryption_key = blst::G1Point::infinity();
        BatchDecryptionPrecomputation btx_precomputation;

        HeadToHeadContext(
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
            c1_points = extract_c1_points(std::span<const Ciphertext>(ciphertexts.data(), ciphertexts.size()));
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

            btx_precomputation = PrecomputeBatchDecryption(
                prepared_public_parameters,
                ciphertexts,
                CrossTermAlgorithm::MiddleProductFft);
        }
    };

    constexpr std::size_t kDefaultThreshold = 8;
    constexpr std::size_t kDefaultServerCount = 16;

    static void BM_BtxVsBeat_Offline_BtxPrepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        HeadToHeadContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xbea70000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(ComputeCrossTerms(
                ctx.prepared_public_parameters,
                ctx.ciphertexts,
                CrossTermAlgorithm::MiddleProductFft));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_BtxVsBeat_Offline_BeatPreparedMiddleProduct(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        HeadToHeadContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xbea80000ULL + batch_size);

        std::vector<blst::GT> out(batch_size, blst::GT::identity());
        blst::MiddleProductWorkspace workspace;

        for (auto _ : state) {
            blst::middle_product(
                std::span<const blst::G1Point>(ctx.c1_points.data(), ctx.c1_points.size()),
                ctx.prepared_public_parameters.fft_precomputed_b_by_batch_size[batch_size],
                std::span<blst::GT>(out.data(), out.size()),
                workspace);
            benchmark::DoNotOptimize(out.data());
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_BtxVsBeat_Online_BtxPrepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        HeadToHeadContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xbea90000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(FinalizeBatchDecryption(
                ctx.params,
                ctx.prepared_public_parameters,
                ctx.btx_precomputation,
                ctx.aggregate_predecryption_key));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_BtxVsBeat_Online_BeatPreparedUpdate(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        HeadToHeadContext ctx(batch_size, kDefaultThreshold, kDefaultServerCount, 0xbeaa0000ULL + batch_size);

        std::vector<blst::GT> out(batch_size, blst::GT::identity());

        for (auto _ : state) {
            blst::update(
                ctx.aggregate_predecryption_key,
                ctx.prepared_public_parameters.alpha_precomputed_b_by_batch_size[batch_size],
                std::span<blst::GT>(out.data(), out.size()));
            benchmark::DoNotOptimize(out.data());
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }
}

BENCHMARK(BM_BtxVsBeat_Offline_BtxPrepared)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_Offline_BeatPreparedMiddleProduct)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_Online_BtxPrepared)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_Online_BeatPreparedUpdate)->Apply(apply_batch_args)->UseRealTime();

#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "BTE/building_blocks/blst_cpp.h"
#include "BTE/building_blocks/hkdf_sha256.h"
#include "BTE/btx/core.h"
#include "BTE/transparent_bte/core.h"

namespace
{
    namespace paper = bte::transparent_bte;
    namespace btx = bte::btx;

    constexpr std::size_t kThreshold = 8;
    constexpr std::size_t kServers = 16;

    class DeterministicRng
    {
      public:
        explicit DeterministicRng(std::uint64_t seed) : state_(seed) {}

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
            unsigned int available = 0;
            for (std::uint8_t& byte : out) {
                if (available == 0) {
                    word = next_u64();
                    available = 8;
                }
                byte = static_cast<std::uint8_t>(word);
                word >>= 8;
                --available;
            }
        }

      private:
        std::uint64_t state_;
    };

    blst::Fr fr(const std::uint64_t value)
    {
        return blst::Fr::from_uint64(value);
    }

    blst::Fr full_width_fr(DeterministicRng& rng)
    {
        std::array<std::uint64_t, 4> limbs = {
            rng.next_u64(), rng.next_u64(), rng.next_u64(), rng.next_u64()
        };
        // A 254-bit nonzero value is below the BLS12-381 scalar modulus and
        // avoids the small-coefficient fixtures that distorted earlier work.
        limbs[3] &= 0x3fffffffffffffffULL;
        limbs[3] |= 0x2000000000000000ULL;
        blst::Fr out{};
        blst::Fr::from_uint64(out, limbs.data());
        return out;
    }

    std::vector<btx::Ciphertext> extract_btx_ciphertexts(
        const std::span<const btx::EncryptionResult> encrypted)
    {
        std::vector<btx::Ciphertext> out;
        out.reserve(encrypted.size());
        for (const btx::EncryptionResult& item : encrypted) {
            out.push_back(item.ciphertext);
        }
        return out;
    }

    void apply_batch_args(benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int batch_size : { 1, 2, 4, 8, 16, 24, 32, 40, 48, 49, 50, 64, 128, 160, 192, 208, 224, 256, 512 }) {
            benchmark_handle->Arg(batch_size);
        }
    }

    void apply_naive_batch_args(benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int batch_size : { 1, 2, 4, 8, 16, 24, 32, 40, 48, 49, 50, 64, 128, 160, 192, 208, 224, 256 }) {
            benchmark_handle->Arg(batch_size);
        }
    }

    void apply_quadratic_pairing_args(
        benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int batch_size : { 1, 2, 4, 8, 16, 24, 32, 40, 48, 49, 50, 64, 96, 128 }) {
            benchmark_handle->Arg(batch_size);
        }
    }

    void apply_msm_batch_args(benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int batch_size : { 1, 2, 4, 8, 16, 24, 32, 40, 48, 49, 50, 64, 128, 160, 192, 208, 224, 256, 384, 512 }) {
            benchmark_handle->Arg(batch_size);
        }
    }

    void apply_msm_profile_args(benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int batch_size : { 16, 32, 50 }) {
            benchmark_handle->Arg(batch_size);
        }
    }

    void set_batch_counter(benchmark::State& state, const std::size_t batch_size)
    {
        state.counters["batch"] = static_cast<double>(batch_size);
        state.SetItemsProcessed(
            static_cast<std::int64_t>(state.iterations()) *
            static_cast<std::int64_t>(batch_size));
    }

    std::vector<std::array<std::uint8_t, 32>> derive_paper_keys(
        const std::span<const blst::GT> plaintexts)
    {
        constexpr std::string_view salt_string = "BTE_2026_2205_KDF_SALT_V1";
        constexpr std::string_view info_string = "BTE_2026_2205_KDF_INFO_V1";
        const auto salt = std::span<const std::uint8_t>(
            reinterpret_cast<const std::uint8_t*>(salt_string.data()), salt_string.size());
        const auto info = std::span<const std::uint8_t>(
            reinterpret_cast<const std::uint8_t*>(info_string.data()), info_string.size());

        std::vector<std::array<std::uint8_t, 32>> keys(plaintexts.size());
        std::array<std::uint8_t, 12 * 48> serialized{};
        for (std::size_t i = 0; i < plaintexts.size(); ++i) {
            blst::GT::serialize(serialized.data(), plaintexts[i]);
            const std::vector<std::uint8_t> key = bte::building_blocks::hkdf_sha256(
                salt, serialized, info, keys[i].size());
            std::copy(key.begin(), key.end(), keys[i].begin());
        }
        return keys;
    }

    paper::BatchPrecomputation precompute_paper_best_measured(
        const std::span<const paper::Ciphertext> ciphertexts)
    {
        // Native row MSM beats both the repeated-multiply oracle and the
        // Lemma 1 implementation at every measured size from B=4 through 512.
        // The one-element batch has no cross terms, while direct multiplication
        // is narrowly faster for the single pair at B=2. Above the measured
        // range, retain the asymptotically faster path.
        constexpr std::size_t largest_measured_msm_win = 512;
        if (ciphertexts.size() == 2) {
            return paper::PrecomputeBatchDecryptionNaive(ciphertexts);
        }
        if (ciphertexts.size() > 2 &&
            ciphertexts.size() <= largest_measured_msm_win) {
            return paper::PrecomputeBatchDecryptionMsm(ciphertexts);
        }
        return paper::PrecomputeBatchDecryptionFast(ciphertexts);
    }

    struct PaperContext
    {
        paper::SetupResult setup;
        std::vector<paper::Ciphertext> ciphertexts;
        std::vector<paper::PartialDecryptionShare> shares;
        blst::G1Point combined = blst::G1Point::infinity();
        paper::BatchPrecomputation precomputation;

        explicit PaperContext(
            const std::size_t batch_size,
            const bool prepare_decryption = true)
        {
            DeterministicRng paper_rng(0x22050000ULL + batch_size);
            std::vector<blst::Fr> sharing;
            sharing.reserve(kThreshold);
            for (std::size_t i = 0; i < kThreshold; ++i) {
                sharing.push_back(full_width_fr(paper_rng));
            }
            constexpr std::array<std::uint8_t, 25> label = {
                '2','0','2','6','/','2','2','0','5','/','g','-','p','r','i','m','e','/','b','e','n','c','h','/','1'
            };
            const blst::G1Point g_prime = bte::building_blocks::hash_to_g1(
                label, "BTE_TRANSPARENT_2026_2205_G_PRIME_V1");
            setup = paper::SetupFromSharingPolynomial(sharing, kServers, g_prime);

            ciphertexts.reserve(batch_size);
            for (std::size_t i = 0; i < batch_size; ++i) {
                ciphertexts.push_back(paper::Encrypt(
                    setup.public_parameters,
                    full_width_fr(paper_rng),
                    blst::GT::identity()));
            }
            shares.reserve(kThreshold);
            for (std::size_t i = 0; i < kThreshold; ++i) {
                shares.push_back(paper::PartialDecrypt(setup.secret_shares[i], ciphertexts));
            }
            combined = paper::CombinePartialDecryptionShares(
                setup.public_parameters, shares);
            if (prepare_decryption) {
                precomputation = precompute_paper_best_measured(ciphertexts);
            }
        }
    };

    struct BtxContext
    {
        btx::BtxParameters params;
        DeterministicRng rng;
        btx::SetupResult setup;
        btx::PreparedPublicParameters prepared;
        std::vector<btx::EncryptionResult> encrypted;
        std::vector<btx::Ciphertext> ciphertexts;
        btx::ThresholdRequest request;
        std::vector<btx::PartialDecryptionShare> shares;
        blst::G1Point combined = blst::G1Point::infinity();
        btx::BatchDecryptionPrecomputation precomputation;

        explicit BtxContext(
            const std::size_t batch_size,
            const bool prepare_decryption = true)
            : params{ batch_size, kServers, kThreshold, 32,
                      "BTE_BTX_KDF_SALT_V1", "BTE_BTX_KDF_INFO_V1" },
              rng(0x22050000ULL + batch_size),
              setup(btx::Setup(
                  params,
                  [this](std::span<std::uint8_t> out) { rng.fill(out); })),
              prepared(btx::PreparePublicParameters(setup.public_parameters))
        {
            encrypted.reserve(batch_size);
            for (std::size_t i = 0; i < batch_size; ++i) {
                encrypted.push_back(btx::Encrypt(
                    params,
                    setup.public_parameters,
                    [this](std::span<std::uint8_t> out) { rng.fill(out); }));
            }
            ciphertexts = extract_btx_ciphertexts(encrypted);
            request = btx::BuildThresholdRequest(setup.public_parameters, ciphertexts);
            shares.reserve(kThreshold);
            for (std::size_t i = 0; i < kThreshold; ++i) {
                shares.push_back(btx::CreatePartialDecryptionShare(
                    setup.secret_shares[i], request));
            }
            combined = btx::CombinePartialDecryptionShares(shares, setup.public_parameters);
            if (prepare_decryption) {
                precomputation = btx::PrecomputeBatchDecryption(
                    prepared, ciphertexts, btx::CrossTermAlgorithm::MiddleProductFft);
                benchmark::DoNotOptimize(btx::FinalizeBatchDecryption(
                    params, prepared, precomputation, combined));
            }
        }
    };

    static void BM_Paper_Encrypt(benchmark::State& state)
    {
        PaperContext ctx(1, false);
        std::uint64_t counter = 500001;
        for (auto _ : state) {
            benchmark::DoNotOptimize(paper::Encrypt(
                ctx.setup.public_parameters,
                fr(counter++),
                blst::GT::identity()));
            benchmark::ClobberMemory();
        }
        state.SetItemsProcessed(state.iterations());
    }

    static void BM_Btx_Encrypt_Comparison(benchmark::State& state)
    {
        BtxContext ctx(1, false);
        DeterministicRng rng(0xb7e00000ULL);
        for (auto _ : state) {
            benchmark::DoNotOptimize(btx::Encrypt(
                ctx.params,
                ctx.setup.public_parameters,
                [&rng](std::span<std::uint8_t> out) { rng.fill(out); }));
            benchmark::ClobberMemory();
        }
        state.SetItemsProcessed(state.iterations());
    }

    static void BM_Paper_PartialDecrypt(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PaperContext ctx(batch_size, false);
        for (auto _ : state) {
            benchmark::DoNotOptimize(paper::PartialDecrypt(
                ctx.setup.secret_shares.front(), ctx.ciphertexts));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_PartialDecrypt_Comparison(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxContext ctx(batch_size, false);
        for (auto _ : state) {
            benchmark::DoNotOptimize(btx::CreatePartialDecryptionShare(
                ctx.setup.secret_shares.front(), ctx.request));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Paper_VerifyShare(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PaperContext ctx(batch_size, false);
        for (auto _ : state) {
            benchmark::DoNotOptimize(paper::VerifyPartialDecryptionShare(
                ctx.setup.public_parameters, ctx.shares.front(), ctx.ciphertexts));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Paper_CombineShares(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PaperContext ctx(batch_size, false);
        for (auto _ : state) {
            benchmark::DoNotOptimize(paper::CombinePartialDecryptionShares(
                ctx.setup.public_parameters, ctx.shares));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_VerifyCombined_Comparison(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxContext ctx(batch_size, false);
        btx::AggregatePredecryptionKeyVerificationWorkspace workspace;
        for (auto _ : state) {
            benchmark::DoNotOptimize(btx::VerifyAggregatePredecryptionKey(
                ctx.request, ctx.combined, ctx.setup.public_parameters, workspace));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_CombineShares_Comparison(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxContext ctx(batch_size, false);
        for (auto _ : state) {
            benchmark::DoNotOptimize(btx::CombinePartialDecryptionShares(
                ctx.shares, ctx.setup.public_parameters));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Paper_Precompute_Naive(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PaperContext ctx(batch_size);
        for (auto _ : state) {
            benchmark::DoNotOptimize(
                paper::PrecomputeBatchDecryptionNaive(ctx.ciphertexts));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Paper_Precompute_Fast(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PaperContext ctx(batch_size);
        for (auto _ : state) {
            benchmark::DoNotOptimize(
                paper::PrecomputeBatchDecryptionFast(ctx.ciphertexts));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Paper_Precompute_QuadraticPairing(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PaperContext ctx(batch_size);
        for (auto _ : state) {
            benchmark::DoNotOptimize(
                paper::PrecomputeBatchDecryptionQuadraticPairing(
                    ctx.ciphertexts));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Paper_Precompute_Msm(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PaperContext ctx(batch_size, false);
        for (auto _ : state) {
            benchmark::DoNotOptimize(
                paper::PrecomputeBatchDecryptionMsm(ctx.ciphertexts));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Paper_Precompute_Msm_Profile(
        benchmark::State& state,
        const paper::MsmWeightMode weight_mode)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PaperContext ctx(batch_size, false);
        paper::MsmPrecomputeTimings sum{};
        for (auto _ : state) {
            paper::MsmPrecomputeTimings current{};
            benchmark::DoNotOptimize(
                paper::BenchmarkPrecomputeBatchDecryptionMsm(
                    ctx.ciphertexts, weight_mode, &current));
            sum.point_hash_ns += current.point_hash_ns;
            sum.output_shell_ns += current.output_shell_ns;
            sum.cauchy_weights_ns += current.cauchy_weights_ns;
            sum.buffer_setup_ns += current.buffer_setup_ns;
            sum.g1_msm_ns += current.g1_msm_ns;
            sum.g2_msm_ns += current.g2_msm_ns;
            sum.affine_ns += current.affine_ns;
            sum.pairing_ns += current.pairing_ns;
            sum.total_ns += current.total_ns;
            benchmark::ClobberMemory();
        }

        const double ns_to_average_ms =
            1.0 / (static_cast<double>(state.iterations()) * 1'000'000.0);
        state.counters["hash_ms"] = sum.point_hash_ns * ns_to_average_ms;
        state.counters["shell_ms"] = sum.output_shell_ns * ns_to_average_ms;
        state.counters["weights_ms"] = sum.cauchy_weights_ns * ns_to_average_ms;
        state.counters["buffers_ms"] = sum.buffer_setup_ns * ns_to_average_ms;
        state.counters["g1_msm_ms"] = sum.g1_msm_ns * ns_to_average_ms;
        state.counters["g2_msm_ms"] = sum.g2_msm_ns * ns_to_average_ms;
        state.counters["affine_ms"] = sum.affine_ns * ns_to_average_ms;
        state.counters["pairing_ms"] = sum.pairing_ns * ns_to_average_ms;
        state.counters["profile_total_ms"] = sum.total_ns * ns_to_average_ms;
        set_batch_counter(state, batch_size);
    }

    static void BM_Paper_Precompute_Msm_Profile_Cauchy(
        benchmark::State& state)
    {
        BM_Paper_Precompute_Msm_Profile(state, paper::MsmWeightMode::Cauchy);
    }

    static void BM_Paper_Precompute_Msm_Profile_Unit(
        benchmark::State& state)
    {
        BM_Paper_Precompute_Msm_Profile(state, paper::MsmWeightMode::Unit);
    }

    static void BM_Paper_Precompute_BestMeasured(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PaperContext ctx(batch_size);
        for (auto _ : state) {
            benchmark::DoNotOptimize(
                precompute_paper_best_measured(ctx.ciphertexts));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_Precompute_Prepared(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxContext ctx(batch_size);
        for (auto _ : state) {
            benchmark::DoNotOptimize(btx::PrecomputeBatchDecryption(
                ctx.prepared,
                ctx.ciphertexts,
                btx::CrossTermAlgorithm::MiddleProductFft));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Paper_Finalize_Algebra(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PaperContext ctx(batch_size);
        for (auto _ : state) {
            benchmark::DoNotOptimize(paper::FinalizeBatchDecryption(
                ctx.precomputation, ctx.combined, ctx.ciphertexts));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Paper_Finalize_WithKdf(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PaperContext ctx(batch_size);
        for (auto _ : state) {
            const auto plaintexts = paper::FinalizeBatchDecryption(
                ctx.precomputation, ctx.combined, ctx.ciphertexts);
            benchmark::DoNotOptimize(derive_paper_keys(plaintexts));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_Finalize_WithKdf(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxContext ctx(batch_size);
        for (auto _ : state) {
            benchmark::DoNotOptimize(btx::FinalizeBatchDecryption(
                ctx.params, ctx.prepared, ctx.precomputation, ctx.combined));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Paper_DecryptTotal_WithKdf(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PaperContext ctx(batch_size);
        for (auto _ : state) {
            const auto precomputation =
                precompute_paper_best_measured(ctx.ciphertexts);
            const auto plaintexts = paper::FinalizeBatchDecryption(
                precomputation, ctx.combined, ctx.ciphertexts);
            benchmark::DoNotOptimize(derive_paper_keys(plaintexts));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Paper_DecryptTotal_QuadraticPairing_WithKdf(
        benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        PaperContext ctx(batch_size);
        for (auto _ : state) {
            const auto precomputation =
                paper::PrecomputeBatchDecryptionQuadraticPairing(
                    ctx.ciphertexts);
            const auto plaintexts = paper::FinalizeBatchDecryption(
                precomputation, ctx.combined, ctx.ciphertexts);
            benchmark::DoNotOptimize(derive_paper_keys(plaintexts));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }

    static void BM_Btx_DecryptTotal_WithKdf(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        BtxContext ctx(batch_size);
        for (auto _ : state) {
            const auto precomputation = btx::PrecomputeBatchDecryption(
                ctx.prepared,
                ctx.ciphertexts,
                btx::CrossTermAlgorithm::MiddleProductFft);
            benchmark::DoNotOptimize(btx::FinalizeBatchDecryption(
                ctx.params, ctx.prepared, precomputation, ctx.combined));
            benchmark::ClobberMemory();
        }
        set_batch_counter(state, batch_size);
    }
}

BENCHMARK(BM_Paper_Encrypt)->UseRealTime();
BENCHMARK(BM_Btx_Encrypt_Comparison)->UseRealTime();
BENCHMARK(BM_Paper_PartialDecrypt)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Btx_PartialDecrypt_Comparison)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Paper_VerifyShare)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Paper_CombineShares)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Btx_VerifyCombined_Comparison)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Btx_CombineShares_Comparison)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Paper_Precompute_Naive)->Apply(apply_naive_batch_args)->UseRealTime();
BENCHMARK(BM_Paper_Precompute_Fast)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Paper_Precompute_QuadraticPairing)
    ->Apply(apply_quadratic_pairing_args)
    ->UseRealTime();
BENCHMARK(BM_Paper_Precompute_Msm)->Apply(apply_msm_batch_args)->UseRealTime();
BENCHMARK(BM_Paper_Precompute_Msm_Profile_Cauchy)
    ->Apply(apply_msm_profile_args)
    ->UseRealTime();
BENCHMARK(BM_Paper_Precompute_Msm_Profile_Unit)
    ->Apply(apply_msm_profile_args)
    ->UseRealTime();
BENCHMARK(BM_Paper_Precompute_BestMeasured)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Btx_Precompute_Prepared)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Paper_Finalize_Algebra)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Paper_Finalize_WithKdf)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Btx_Finalize_WithKdf)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Paper_DecryptTotal_WithKdf)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_Paper_DecryptTotal_QuadraticPairing_WithKdf)
    ->Apply(apply_quadratic_pairing_args)
    ->UseRealTime();
BENCHMARK(BM_Btx_DecryptTotal_WithKdf)->Apply(apply_batch_args)->UseRealTime();

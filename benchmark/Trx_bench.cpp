#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "BTE/building_blocks/blst_cpp.h"
#include "BTE/trx/core.h"

using namespace bte::trx;

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

    blst::Fr fr_from_u64(const std::uint64_t value)
    {
        return blst::Fr::from_uint64(value);
    }

    blst::Fr fr_add(const blst::Fr& lhs, const blst::Fr& rhs)
    {
        blst::Fr out{};
        blst::Fr::add(out, lhs, rhs);
        return out;
    }

    blst::Fr fr_mul(const blst::Fr& lhs, const blst::Fr& rhs)
    {
        blst::Fr out{};
        blst::Fr::mul(out, lhs, rhs);
        return out;
    }

    blst::Fr random_fr(DeterministicRng& rng)
    {
        std::array<std::uint8_t, 48> seed{};
        rng.fill(seed);
        const std::array<std::uint8_t, 32> scalar_bytes =
            bte::building_blocks::hash_to_scalar_le(
                std::span<const std::uint8_t>(seed.data(), seed.size()),
                "BTE_TRX_BENCH_RNG_FR_V1");
        blst::Fr out{};
        blst::Fr::from_scalar(out, scalar_bytes.data());
        return out;
    }

    std::vector<ValidatorSecretShare> make_shamir_secret_shares(
        const std::uint32_t validator_count,
        const std::uint32_t threshold,
        DeterministicRng& rng)
    {
        if (threshold == 0 || threshold > validator_count) {
            throw std::invalid_argument("Invalid threshold configuration");
        }

        std::vector<blst::Fr> coefficients;
        coefficients.reserve(threshold);
        for (std::uint32_t i = 0; i < threshold; ++i) {
            coefficients.push_back(random_fr(rng));
        }

        std::vector<ValidatorSecretShare> shares;
        shares.reserve(validator_count);
        for (std::uint32_t index = 1; index <= validator_count; ++index) {
            const blst::Fr x = fr_from_u64(index);
            blst::Fr power = blst::Fr::one();
            blst::Fr value = blst::Fr::zero();
            for (const blst::Fr& coefficient : coefficients) {
                value = fr_add(value, fr_mul(coefficient, power));
                power = fr_mul(power, x);
            }

            ValidatorSecretShare share{};
            share.index = index;
            share.sk_share = value;
            shares.push_back(share);
        }
        return shares;
    }

    Bytes make_associated_data(const std::size_t index)
    {
        return Bytes{
            static_cast<std::uint8_t>('a'),
            static_cast<std::uint8_t>('d'),
            static_cast<std::uint8_t>('-'),
            static_cast<std::uint8_t>(index & 0xffU),
            static_cast<std::uint8_t>((index >> 8U) & 0xffU)
        };
    }

    std::vector<CiphertextEnvelope> to_envelopes(const std::vector<EncryptedCiphertext>& encrypted)
    {
        std::vector<CiphertextEnvelope> envelopes;
        envelopes.reserve(encrypted.size());
        for (const EncryptedCiphertext& item : encrypted) {
            envelopes.push_back({ item.ciphertext, item.associated_data });
        }
        return envelopes;
    }

    struct TrxBenchmarkContext
    {
        TrxParameters params{};
        DeterministicRng rng;
        Crs crs;
        std::vector<ValidatorSecretShare> secret_shares;
        std::vector<ValidatorPublicShare> public_shares;
        PublicKeyMaterial public_key_material;
        EncryptionKey encryption_key;
        std::vector<EncryptedCiphertext> encrypted;
        std::vector<CiphertextEnvelope> envelopes;
        std::vector<VerifiedCiphertext> verified_ciphertexts;
        BatchDigest digest;
        std::vector<blst::G1Point> proofs;
        std::vector<PartialDecryptionShare> threshold_shares;
        blst::G1Point sigma = blst::G1Point::infinity();

        TrxBenchmarkContext(
            const std::size_t batch_size,
            const std::uint32_t threshold,
            const std::uint32_t validator_count,
            const std::size_t num_contexts,
            const std::uint64_t seed)
            : params{ std::max<std::size_t>(batch_size, 1), num_contexts, 32, 1,
                      "BTE_TRX_H1_V1",
                      "BTE_TRX_HF_V1",
                      "BTE_TRX_KDF_SALT_V1",
                      "BTE_TRX_KDF_INFO_V1",
                      "BTE_TRX_OTS_BLS_SIG_V1" },
              rng(seed),
              crs(GenerateCRS(params, [this](std::span<std::uint8_t> out) { rng.fill(out); }))
        {
            secret_shares = make_shamir_secret_shares(validator_count, threshold, rng);
            public_shares.reserve(secret_shares.size());
            for (const ValidatorSecretShare& share : secret_shares) {
                public_shares.push_back(DeriveValidatorPublicShare(share, crs.h_tau));
            }

            std::vector<std::uint32_t> threshold_subset;
            threshold_subset.reserve(threshold);
            for (std::uint32_t index = 1; index <= threshold; ++index) {
                threshold_subset.push_back(index);
            }

            public_key_material = DerivePublicKeyMaterial(
                std::span<const ValidatorPublicShare>(public_shares.data(), public_shares.size()),
                std::span<const std::uint32_t>(threshold_subset.data(), threshold_subset.size()));
            encryption_key = ToEncryptionKey(public_key_material);

            encrypted.reserve(batch_size);
            for (std::size_t i = 0; i < batch_size; ++i) {
                encrypted.push_back(Encrypt(
                    params,
                    encryption_key,
                    make_associated_data(i),
                    [this](std::span<std::uint8_t> out) { rng.fill(out); }));
            }
            envelopes = to_envelopes(encrypted);

            verified_ciphertexts.reserve(encrypted.size());
            for (const EncryptedCiphertext& item : encrypted) {
                const std::optional<VerifiedCiphertext> verified = VerifyCiphertext(
                    params,
                    item.ciphertext,
                    std::span<const std::uint8_t>(
                        item.associated_data.data(),
                        item.associated_data.size()));
                if (!verified.has_value()) {
                    throw std::runtime_error("TRX benchmark generated an invalid ciphertext");
                }
                verified_ciphertexts.push_back(*verified);
            }

            digest = BuildBatchDigestFromVerified(
                params,
                std::span<const VerifiedCiphertext>(
                    verified_ciphertexts.data(),
                    verified_ciphertexts.size()),
                0,
                crs);
            proofs = ComputeEvalProofs(digest, crs);

            threshold_shares.reserve(threshold);
            for (std::uint32_t i = 0; i < threshold; ++i) {
                threshold_shares.push_back(CreatePartialDecryption(
                    params,
                    secret_shares[i],
                    digest,
                    public_key_material));
            }
            sigma = CombineShares(
                std::span<const PartialDecryptionShare>(
                    threshold_shares.data(),
                    threshold_shares.size()),
                digest,
                public_key_material);
        }
    };

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
        for (const int threshold : { 2, 4, 8, 16, 32, 512 }) {
            benchmark_handle->Arg(threshold);
        }
    }

    constexpr std::uint32_t kDefaultValidatorCount = 32;
    constexpr std::uint32_t kDefaultThreshold = 8;
    constexpr std::size_t kDefaultContexts = 4;

    static void BM_Trx_GenerateCRS(benchmark::State& state)
    {
        const std::size_t max_batch_size = static_cast<std::size_t>(state.range(0));
        TrxParameters params{
            max_batch_size, kDefaultContexts, 32, 1,
            "BTE_TRX_H1_V1",
            "BTE_TRX_HF_V1",
            "BTE_TRX_KDF_SALT_V1",
            "BTE_TRX_KDF_INFO_V1",
            "BTE_TRX_OTS_BLS_SIG_V1"
        };
        DeterministicRng rng(0xabc00000ULL + max_batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(GenerateCRS(
                params,
                [&rng](std::span<std::uint8_t> out) { rng.fill(out); }));
            benchmark::ClobberMemory();
        }

        state.counters["batch"] = static_cast<double>(max_batch_size);
        state.counters["contexts"] = static_cast<double>(kDefaultContexts);
    }

    static void BM_Trx_DerivePublicKeyMaterial(benchmark::State& state)
    {
        const std::uint32_t threshold = static_cast<std::uint32_t>(state.range(0));
        const std::uint32_t validator_count = std::max(kDefaultValidatorCount, threshold);
        TrxBenchmarkContext ctx(8, threshold, validator_count, 2, 0x11110000ULL + threshold);

        std::vector<std::uint32_t> subset;
        subset.reserve(threshold);
        for (std::uint32_t i = 1; i <= threshold; ++i) {
            subset.push_back(i);
        }

        for (auto _ : state) {
            benchmark::DoNotOptimize(DerivePublicKeyMaterial(
                std::span<const ValidatorPublicShare>(ctx.public_shares.data(), ctx.public_shares.size()),
                std::span<const std::uint32_t>(subset.data(), subset.size())));
            benchmark::ClobberMemory();
        }

        set_threshold_counter(state, threshold);
    }

    static void BM_Trx_Encrypt(benchmark::State& state)
    {
        TrxBenchmarkContext ctx(8, kDefaultThreshold, kDefaultValidatorCount, 2, 0x22220000ULL);
        DeterministicRng rng(0x2222abcdULL);
        const Bytes associated_data = make_associated_data(0);

        for (auto _ : state) {
            benchmark::DoNotOptimize(Encrypt(
                ctx.params,
                ctx.encryption_key,
                associated_data,
                [&rng](std::span<std::uint8_t> out) { rng.fill(out); }));
            benchmark::ClobberMemory();
        }

        state.counters["batch"] = 1.0;
        state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()));
    }

    static void BM_Trx_VerifyCiphertext(benchmark::State& state)
    {
        TrxBenchmarkContext ctx(8, kDefaultThreshold, kDefaultValidatorCount, 2, 0x33330000ULL);
        const EncryptedCiphertext& encrypted = ctx.encrypted.front();

        for (auto _ : state) {
            benchmark::DoNotOptimize(VerifyCiphertext(
                ctx.params,
                encrypted.ciphertext,
                encrypted.associated_data));
            benchmark::ClobberMemory();
        }

        state.counters["batch"] = 1.0;
        state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()));
    }

    static void BM_Trx_BuildBatchDigest(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        TrxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultValidatorCount, 2, 0x44440000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(BuildBatchDigest(
                ctx.params,
                std::span<const CiphertextEnvelope>(ctx.envelopes.data(), ctx.envelopes.size()),
                0,
                ctx.crs));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Trx_ComputeEvalProofs(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        TrxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultValidatorCount, 2, 0x55550000ULL + batch_size);

        for (auto _ : state) {
            BatchDigest digest = ctx.digest;
            digest.has_cached_eval_proofs = false;
            digest.cached_eval_proofs.clear();
            benchmark::DoNotOptimize(ComputeEvalProofs(digest, ctx.crs));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Trx_CreatePartialDecryption(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        TrxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultValidatorCount, 2, 0x66660000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(CreatePartialDecryption(
                ctx.params,
                ctx.secret_shares.front(),
                ctx.digest,
                ctx.public_key_material));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Trx_VerifyPartialDecryption(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        TrxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultValidatorCount, 2, 0x77770000ULL + batch_size);
        const PartialDecryptionShare& share = ctx.threshold_shares.front();
        const ValidatorPublicShare& public_share = ctx.public_shares.front();

        for (auto _ : state) {
            benchmark::DoNotOptimize(VerifyPartialDecryption(
                ctx.params,
                share,
                public_share,
                ctx.digest,
                ctx.public_key_material));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Trx_CombineShares(benchmark::State& state)
    {
        const std::uint32_t threshold = static_cast<std::uint32_t>(state.range(0));
        const std::uint32_t validator_count = std::max(kDefaultValidatorCount, threshold);
        TrxBenchmarkContext ctx(32, threshold, validator_count, 2, 0x88880000ULL + threshold);

        for (auto _ : state) {
            benchmark::DoNotOptimize(CombineShares(
                std::span<const PartialDecryptionShare>(
                    ctx.threshold_shares.data(),
                    ctx.threshold_shares.size()),
                ctx.digest,
                ctx.public_key_material));
            benchmark::ClobberMemory();
        }

        set_threshold_counter(state, threshold);
    }

    static void BM_Trx_DecapsulateOne(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        TrxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultValidatorCount, 2, 0x99990000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(DecapsulateOne(
                ctx.params,
                ctx.digest,
                0,
                ctx.proofs.front(),
                ctx.sigma,
                ctx.crs,
                ctx.public_key_material));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Trx_Phase_CtxtCheck(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        TrxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultValidatorCount, 2, 0x33340000ULL + batch_size);

        for (auto _ : state) {
            for (const EncryptedCiphertext& encrypted : ctx.encrypted) {
                benchmark::DoNotOptimize(VerifyCiphertext(
                    ctx.params,
                    encrypted.ciphertext,
                    encrypted.associated_data));
            }
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Trx_Phase_Precompute(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        TrxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultValidatorCount, 2, 0x44450000ULL + batch_size);

        for (auto _ : state) {
            BatchDigest digest = BuildBatchDigestFromVerified(
                ctx.params,
                std::span<const VerifiedCiphertext>(
                    ctx.verified_ciphertexts.data(),
                    ctx.verified_ciphertexts.size()),
                0,
                ctx.crs);
            benchmark::DoNotOptimize(ComputeEvalProofs(digest, ctx.crs));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Trx_Phase_PartialDecrypt(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        TrxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultValidatorCount, 2, 0x66670000ULL + batch_size);

        for (auto _ : state) {
            benchmark::DoNotOptimize(CreatePartialDecryption(
                ctx.params,
                ctx.secret_shares.front(),
                ctx.digest,
                ctx.public_key_material));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Trx_Phase_ServerCheck(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        TrxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultValidatorCount, 2, 0x77780000ULL + batch_size);
        const PartialDecryptionShare& share = ctx.threshold_shares.front();
        const ValidatorPublicShare& public_share = ctx.public_shares.front();

        for (auto _ : state) {
            benchmark::DoNotOptimize(VerifyPartialDecryption(
                ctx.params,
                share,
                public_share,
                ctx.digest,
                ctx.public_key_material));
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Trx_Phase_Combine(benchmark::State& state)
    {
        const std::uint32_t threshold = static_cast<std::uint32_t>(state.range(0));
        const std::uint32_t validator_count = std::max(kDefaultValidatorCount, threshold);
        TrxBenchmarkContext ctx(32, threshold, validator_count, 2, 0x88890000ULL + threshold);

        for (auto _ : state) {
            benchmark::DoNotOptimize(CombineShares(
                std::span<const PartialDecryptionShare>(
                    ctx.threshold_shares.data(),
                    ctx.threshold_shares.size()),
                ctx.digest,
                ctx.public_key_material));
            benchmark::ClobberMemory();
        }

        set_threshold_counter(state, threshold);
    }

    static void BM_Trx_Phase_Open(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        TrxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultValidatorCount, 2, 0x999a0000ULL + batch_size);

        for (auto _ : state) {
            for (std::size_t i = 0; i < ctx.proofs.size(); ++i) {
                benchmark::DoNotOptimize(DecapsulateOne(
                    ctx.params,
                    ctx.digest,
                    i,
                    ctx.proofs[i],
                    ctx.sigma,
                    ctx.crs,
                    ctx.public_key_material));
            }
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }

    static void BM_Trx_EndToEndBatch(benchmark::State& state)
    {
        const std::size_t batch_size = static_cast<std::size_t>(state.range(0));
        TrxBenchmarkContext ctx(batch_size, kDefaultThreshold, kDefaultValidatorCount, 2, 0xabcd0000ULL + batch_size);

        for (auto _ : state) {
            const BatchDigest digest = BuildBatchDigest(
                ctx.params,
                std::span<const CiphertextEnvelope>(ctx.envelopes.data(), ctx.envelopes.size()),
                0,
                ctx.crs);
            const std::vector<blst::G1Point> proofs = ComputeEvalProofs(digest, ctx.crs);

            std::vector<PartialDecryptionShare> shares;
            shares.reserve(ctx.threshold_shares.size());
            for (std::size_t i = 0; i < ctx.threshold_shares.size(); ++i) {
                shares.push_back(CreatePartialDecryption(
                    ctx.params,
                    ctx.secret_shares[i],
                    digest,
                    ctx.public_key_material));
            }

            const blst::G1Point sigma = CombineShares(
                std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
                digest,
                ctx.public_key_material);

            for (std::size_t i = 0; i < proofs.size(); ++i) {
                benchmark::DoNotOptimize(DecapsulateOne(
                    ctx.params,
                    digest,
                    i,
                    proofs[i],
                    sigma,
                    ctx.crs,
                    ctx.public_key_material));
            }
            benchmark::ClobberMemory();
        }

        set_batch_counter(state, batch_size);
    }
}

BENCHMARK(BM_Trx_GenerateCRS)->Apply(apply_batch_args);
BENCHMARK(BM_Trx_DerivePublicKeyMaterial)->Apply(apply_threshold_args);
BENCHMARK(BM_Trx_Encrypt);
BENCHMARK(BM_Trx_VerifyCiphertext);
BENCHMARK(BM_Trx_BuildBatchDigest)->Apply(apply_batch_args);
BENCHMARK(BM_Trx_ComputeEvalProofs)->Apply(apply_batch_args);
BENCHMARK(BM_Trx_CreatePartialDecryption)->Apply(apply_batch_args);
BENCHMARK(BM_Trx_VerifyPartialDecryption)->Apply(apply_batch_args);
BENCHMARK(BM_Trx_CombineShares)->Apply(apply_threshold_args);
BENCHMARK(BM_Trx_DecapsulateOne)->Apply(apply_batch_args);
BENCHMARK(BM_Trx_EndToEndBatch)->Apply(apply_batch_args);
BENCHMARK(BM_Trx_Phase_CtxtCheck)->Apply(apply_batch_args);
BENCHMARK(BM_Trx_Phase_Precompute)->Apply(apply_batch_args);
BENCHMARK(BM_Trx_Phase_PartialDecrypt)->Apply(apply_batch_args);
BENCHMARK(BM_Trx_Phase_ServerCheck)->Apply(apply_batch_args);
BENCHMARK(BM_Trx_Phase_Combine)->Apply(apply_threshold_args);
BENCHMARK(BM_Trx_Phase_Open)->Apply(apply_batch_args);

#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "BTE/building_blocks/blst_cpp.h"
#include "BTE/btx/core.h"
#if BLST_DEFAULT_AVX
#include "blst/avx/curve/map_to_G.h"
#endif
#include "blst/curve/multi_scalar.h"
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

    blst::Fr fr_zero()
    {
        blst::Fr zero{};
        std::memset(&zero, 0, sizeof(zero));
        return zero;
    }

    bool fr_is_zero(const blst::Fr& value)
    {
        return blst::Fr::equal(value, fr_zero());
    }

    std::array<std::uint8_t, 32> fr_to_scalar_bytes(const blst::Fr& scalar)
    {
        std::array<std::uint8_t, 32> bytes{};
        blst::Fr::to_scaler(bytes.data(), scalar);
        return bytes;
    }

    blst::Fr bytes_to_fr(const std::array<std::uint8_t, 32>& bytes)
    {
        blst::Fr out{};
        blst::Fr::from_scalar(out, bytes.data());
        return out;
    }

    blst::Fr random_nonzero_fr(DeterministicRng& rng, const std::string& dst)
    {
        std::array<std::uint8_t, 48> seed{};
        std::uint32_t counter = 0;
        for (;;) {
            rng.fill(seed);
            seed[0] ^= static_cast<std::uint8_t>(counter);
            seed[1] ^= static_cast<std::uint8_t>(counter >> 8);
            seed[2] ^= static_cast<std::uint8_t>(counter >> 16);
            seed[3] ^= static_cast<std::uint8_t>(counter >> 24);
            ++counter;

            const blst::Fr out = bytes_to_fr(bte::building_blocks::hash_to_scalar_le(seed, dst));
            if (!fr_is_zero(out)) {
                return out;
            }
        }
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

    blst::Fr fr_neg(const blst::Fr& value)
    {
        blst::Fr out{};
        blst::Fr::neg(out, value);
        return out;
    }

    blst::G1Point add_g1(const blst::G1Point& lhs, const blst::G1Point& rhs)
    {
        blst::G1Point out{};
        blst::G1Point::dadd(out, lhs, rhs);
        return out;
    }

    blst::G1Point scale_g1(const blst::G1Point& point, const blst::Fr& scalar)
    {
        if (fr_is_zero(scalar)) {
            return blst::G1Point::infinity();
        }

        const std::array<std::uint8_t, 32> scalar_bytes = fr_to_scalar_bytes(scalar);
        blst::G1Point out{};
        blst::G1Point::mult(out, point, scalar_bytes.data(), 256);
        return out;
    }

    bool is_identity_g1(const blst::G1Point& point)
    {
        return blst::G1Point::equal(point, blst::G1Point::infinity()) != 0;
    }

    struct PreparedG1Msm
    {
        std::vector<blst::G1Point> direct_points;
        std::vector<blst::Fr> direct_scalars;
        std::vector<blst::scalar::G1AffinePoint> affine_points;
        std::vector<std::uint8_t> scalar_blob;
        mutable std::vector<blst::Pointxyzz<blst::Fp>> scratch;
        bool use_pippenger = false;

        void prepare(
            const std::span<const blst::G1Point> points,
            const std::span<const blst::Fr> scalars)
        {
            if (points.size() != scalars.size()) {
                throw std::invalid_argument("MSM inputs must have matching lengths");
            }
            use_pippenger = points.size() >= 8;
            if (!use_pippenger) {
                direct_points.assign(points.begin(), points.end());
                direct_scalars.assign(scalars.begin(), scalars.end());
                affine_points.clear();
                scalar_blob.clear();
                scratch.clear();
                return;
            }

            direct_points.clear();
            direct_scalars.clear();
            affine_points.resize(points.size());
            scalar_blob.assign(points.size() * 32, 0);
            for (std::size_t i = 0; i < points.size(); ++i) {
#if BLST_DEFAULT_AVX
                const blst::scalar::G1Point scalar_point =
                    blst::avx::avx_map_to_G_detail::to_scalar(points[i]);
                blst::scalar::G1Point::to_affine(affine_points[i], scalar_point);
#else
                blst::scalar::G1Point::to_affine(affine_points[i], points[i]);
#endif
                const std::array<std::uint8_t, 32> scalar_bytes = fr_to_scalar_bytes(scalars[i]);
                std::memcpy(scalar_blob.data() + (i * 32), scalar_bytes.data(), scalar_bytes.size());
            }

            const std::size_t scratch_bytes =
                blst::msm_mult_pippenger_scratch_sizeof_<blst::Pointxyzz<blst::Fp>>(points.size());
            scratch.resize(scratch_bytes / sizeof(blst::Pointxyzz<blst::Fp>));
        }

        blst::G1Point run() const
        {
            if (!use_pippenger) {
                blst::G1Point out = blst::G1Point::infinity();
                for (std::size_t i = 0; i < direct_points.size(); ++i) {
                    out = add_g1(out, scale_g1(direct_points[i], direct_scalars[i]));
                }
                return out;
            }

            const blst::scalar::G1AffinePoint* point_params[2] = { affine_points.data(), nullptr };
            const blst::byte* scalar_params[2] = { scalar_blob.data(), nullptr };
            blst::scalar::G1Point scalar_out{};
            blst::msm_mult_pippenger<
                blst::scalar::Point<blst::scalar::Fp>,
                blst::scalar::AffinePoint<blst::scalar::Fp>,
                blst::scalar::Pointxyzz<blst::scalar::Fp>>(
                scalar_out,
                point_params,
                affine_points.size(),
                scalar_params,
                256,
                scratch.data(),
                blst::pippenger_window_size(affine_points.size()));
#if BLST_DEFAULT_AVX
            return blst::avx::avx_map_to_G_detail::from_scalar(scalar_out);
#else
            return scalar_out;
#endif
        }
    };

    void append_g1(std::vector<std::uint8_t>& out, const blst::G1Point& point)
    {
        std::array<std::uint8_t, 48> compressed{};
        blst::G1Point::compress(compressed.data(), point);
        out.insert(out.end(), compressed.begin(), compressed.end());
    }

    blst::Fr beat_dleq_challenge(
        const blst::G1Point& g,
        const blst::G1Point& h,
        const blst::G1Point& x,
        const blst::G1Point& y,
        const blst::G1Point& a,
        const blst::G1Point& b)
    {
        std::vector<std::uint8_t> transcript;
        transcript.reserve(6 * 48);
        append_g1(transcript, g);
        append_g1(transcript, h);
        append_g1(transcript, x);
        append_g1(transcript, y);
        append_g1(transcript, a);
        append_g1(transcript, b);
        return bytes_to_fr(bte::building_blocks::hash_to_scalar_le(
            transcript,
            "BTE_BEAT_CLIENT_DLEQ_V1"));
    }

    struct BeatDleqProof
    {
        blst::Fr challenge{};
        blst::Fr response{};
    };

    struct BeatClientCiphertextCheck
    {
        blst::G1Point h = blst::G1Point::infinity();
        blst::G1Point x = blst::G1Point::infinity();
        blst::G1Point y = blst::G1Point::infinity();
        blst::G1Point a = blst::G1Point::infinity();
        blst::G1Point b = blst::G1Point::infinity();
        BeatDleqProof proof{};
    };

    BeatClientCiphertextCheck make_beat_client_ciphertext_check(
        DeterministicRng& rng,
        const std::size_t index)
    {
        const std::string suffix = std::to_string(index);
        const blst::G1Point g = blst::G1Point::generator();
        const blst::Fr delta = random_nonzero_fr(rng, "BTE_BEAT_DLEQ_DELTA_" + suffix);
        const blst::Fr r = random_nonzero_fr(rng, "BTE_BEAT_DLEQ_WITNESS_" + suffix);
        const blst::Fr nonce = random_nonzero_fr(rng, "BTE_BEAT_DLEQ_NONCE_" + suffix);

        BeatClientCiphertextCheck out{};
        out.h = scale_g1(g, delta);
        out.x = scale_g1(g, r);
        out.y = scale_g1(out.h, r);

        out.a = scale_g1(g, nonce);
        out.b = scale_g1(out.h, nonce);
        out.proof.challenge = beat_dleq_challenge(g, out.h, out.x, out.y, out.a, out.b);
        out.proof.response = fr_add(nonce, fr_mul(out.proof.challenge, r));
        return out;
    }

    bool verify_beat_client_ciphertext_check(const BeatClientCiphertextCheck& check)
    {
        const blst::G1Point g = blst::G1Point::generator();
        const blst::Fr challenge =
            beat_dleq_challenge(g, check.h, check.x, check.y, check.a, check.b);
        if (!blst::Fr::equal(challenge, check.proof.challenge)) {
            return false;
        }

        const blst::G1Point lhs_g = scale_g1(g, check.proof.response);
        const blst::G1Point rhs_g = add_g1(check.a, scale_g1(check.x, check.proof.challenge));
        const blst::G1Point lhs_h = scale_g1(check.h, check.proof.response);
        const blst::G1Point rhs_h = add_g1(check.b, scale_g1(check.y, check.proof.challenge));
        return blst::G1Point::equal(lhs_g, rhs_g) != 0 &&
               blst::G1Point::equal(lhs_h, rhs_h) != 0;
    }

    struct BeatClientBatchCheck
    {
        std::vector<BeatClientCiphertextCheck> checks;
        std::vector<blst::Fr> weights;
        PreparedG1Msm first_equation;
        PreparedG1Msm second_equation;
    };

    void prepare_beat_client_batch_check(BeatClientBatchCheck& batch)
    {
        const blst::G1Point g = blst::G1Point::generator();
        blst::Fr g_scalar = fr_zero();
        std::vector<blst::G1Point> first_points;
        std::vector<blst::Fr> first_scalars;
        std::vector<blst::G1Point> second_points;
        std::vector<blst::Fr> second_scalars;
        first_points.reserve(2 * batch.checks.size() + 1);
        first_scalars.reserve(2 * batch.checks.size() + 1);
        second_points.reserve(3 * batch.checks.size());
        second_scalars.reserve(3 * batch.checks.size());

        for (std::size_t i = 0; i < batch.checks.size(); ++i) {
            const BeatClientCiphertextCheck& check = batch.checks[i];
            const blst::Fr weighted_response = fr_mul(batch.weights[i], check.proof.response);
            const blst::Fr weighted_challenge = fr_mul(batch.weights[i], check.proof.challenge);
            g_scalar = fr_add(g_scalar, weighted_response);

            first_points.push_back(check.x);
            first_scalars.push_back(fr_neg(weighted_challenge));
            first_points.push_back(check.a);
            first_scalars.push_back(fr_neg(batch.weights[i]));

            second_points.push_back(check.h);
            second_scalars.push_back(weighted_response);
            second_points.push_back(check.y);
            second_scalars.push_back(fr_neg(weighted_challenge));
            second_points.push_back(check.b);
            second_scalars.push_back(fr_neg(batch.weights[i]));
        }

        first_points.push_back(g);
        first_scalars.push_back(g_scalar);

        batch.first_equation.prepare(
            std::span<const blst::G1Point>(first_points.data(), first_points.size()),
            std::span<const blst::Fr>(first_scalars.data(), first_scalars.size()));
        batch.second_equation.prepare(
            std::span<const blst::G1Point>(second_points.data(), second_points.size()),
            std::span<const blst::Fr>(second_scalars.data(), second_scalars.size()));
    }

    BeatClientBatchCheck make_beat_client_batch_check(
        DeterministicRng& rng,
        const std::size_t client_index,
        const std::size_t copies)
    {
        BeatClientBatchCheck out{};
        out.checks.reserve(copies);
        out.weights.reserve(copies);
        for (std::size_t i = 0; i < copies; ++i) {
            out.checks.push_back(
                make_beat_client_ciphertext_check(rng, client_index * copies + i));
            out.weights.push_back(random_nonzero_fr(
                rng,
                "BTE_BEAT_BATCH_DLEQ_WEIGHT_" + std::to_string(client_index) + "_" +
                    std::to_string(i)));
        }
        prepare_beat_client_batch_check(out);
        return out;
    }

    bool verify_beat_client_batch_check(const BeatClientBatchCheck& batch)
    {
        if (batch.checks.size() != batch.weights.size()) {
            return false;
        }

        const blst::G1Point g = blst::G1Point::generator();
        for (std::size_t i = 0; i < batch.checks.size(); ++i) {
            const BeatClientCiphertextCheck& check = batch.checks[i];
            const blst::Fr challenge =
                beat_dleq_challenge(g, check.h, check.x, check.y, check.a, check.b);
            if (!blst::Fr::equal(challenge, check.proof.challenge)) {
                return false;
            }
        }

        const blst::G1Point first_result = batch.first_equation.run();
        if (!is_identity_g1(first_result)) {
            return false;
        }

        const blst::G1Point second_result = batch.second_equation.run();
        return is_identity_g1(second_result);
    }

    std::vector<BeatClientBatchCheck> make_beat_client_batch_checks(
        const std::size_t client_count,
        const std::size_t copies,
        const std::uint64_t seed)
    {
        DeterministicRng rng(seed);
        std::vector<BeatClientBatchCheck> checks;
        checks.reserve(client_count);
        for (std::size_t i = 0; i < client_count; ++i) {
            checks.push_back(make_beat_client_batch_check(rng, i, copies));
        }
        return checks;
    }

    void set_batch_counter(benchmark::State& state, const std::size_t batch_size)
    {
        state.counters["batch"] = static_cast<double>(batch_size);
        state.SetItemsProcessed(
            static_cast<std::int64_t>(state.iterations()) *
            static_cast<std::int64_t>(batch_size));
    }

    struct BeatCuckooParameters
    {
        std::size_t batch_size;
        std::size_t copies;
        std::size_t subbatch_count;
        std::size_t subbatch_size;
        std::size_t supported_batch_size;
        double security_bits;
    };

    constexpr std::array<BeatCuckooParameters, 5> kBeat40BitParameters = {{
        { 32, 2, 13, 8, 38, 40.0 },
        { 64, 2, 16, 8, 69, 40.0 },
        { 128, 2, 14, 16, 162, 40.0 },
        { 256, 2, 12, 32, 351, 40.0 },
        { 512, 2, 10, 64, 606, 40.0 },
    }};

    constexpr std::array<BeatCuckooParameters, 5> kBeat128BitParameters = {{
        { 32, 20, 8, 8, 32, 512.3 },
        { 64, 20, 9, 8, 64, 243.6 },
        { 128, 20, 9, 16, 128, 475.5 },
        { 256, 20, 9, 32, 256, 936.6 },
        { 512, 20, 9, 64, 512, 1125.9 },
    }};

    BeatCuckooParameters select_parameters(
        const std::span<const BeatCuckooParameters> table,
        const std::size_t batch_size)
    {
        const auto it = std::find_if(
            table.begin(),
            table.end(),
            [batch_size](const BeatCuckooParameters& params) {
                return params.batch_size == batch_size;
            });
        if (it == table.end()) {
            throw std::invalid_argument("Unsupported BEAT++ benchmark batch size");
        }
        return *it;
    }

    void set_beat_cuckoo_counters(
        benchmark::State& state,
        const BeatCuckooParameters& params)
    {
        const std::size_t slots = params.subbatch_count * params.subbatch_size;
        state.counters["batch"] = static_cast<double>(params.batch_size);
        state.counters["copies"] = static_cast<double>(params.copies);
        state.counters["subbatches"] = static_cast<double>(params.subbatch_count);
        state.counters["subbatch"] = static_cast<double>(params.subbatch_size);
        state.counters["slots"] = static_cast<double>(slots);
        state.counters["load"] = static_cast<double>(params.batch_size) / static_cast<double>(slots);
        state.counters["supported"] = static_cast<double>(params.supported_batch_size);
        state.counters["security_bits"] = params.security_bits;
        state.SetItemsProcessed(
            static_cast<std::int64_t>(state.iterations()) *
            static_cast<std::int64_t>(params.batch_size));
    }

    void apply_batch_args(benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int batch_size : { 1, 8, 32, 64, 128, 256, 512 }) {
            benchmark_handle->Arg(batch_size);
        }
    }

    void apply_paper_batch_args(benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int batch_size : { 32, 64, 128, 256, 512 }) {
            benchmark_handle->Arg(batch_size);
        }
    }

    void apply_paper_batch_threshold_args(benchmark::internal::Benchmark* benchmark_handle)
    {
        for (const int batch_size : { 32, 64, 128, 256, 512 }) {
            for (const int threshold : { 2, 4, 8, 16, 512 }) {
                benchmark_handle->Args({ batch_size, threshold });
            }
        }
    }

    void set_beat_cuckoo_threshold_counters(
        benchmark::State& state,
        const BeatCuckooParameters& params,
        const std::size_t threshold)
    {
        set_beat_cuckoo_counters(state, params);
        state.counters["threshold"] = static_cast<double>(threshold);
        state.counters["share_sets"] = static_cast<double>(params.subbatch_count);
        state.counters["threshold_shares"] =
            static_cast<double>(params.subbatch_count * threshold);
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
            benchmark::DoNotOptimize(FinalizeBatchDecryption(
                params,
                prepared_public_parameters,
                btx_precomputation,
                aggregate_predecryption_key));
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

    void run_beat_plus_plus_precompute(
        benchmark::State& state,
        const BeatCuckooParameters& params,
        const std::uint64_t seed)
    {
        HeadToHeadContext ctx(
            params.subbatch_size,
            kDefaultThreshold,
            kDefaultServerCount,
            seed + params.batch_size);

        std::vector<blst::GT> out(params.subbatch_size, blst::GT::identity());
        blst::MiddleProductWorkspace workspace;
        const blst::MiddleProductPrecomputedB& precomputed_b =
            ctx.prepared_public_parameters.fft_precomputed_b_by_batch_size[params.subbatch_size];

        for (auto _ : state) {
            for (std::size_t i = 0; i < params.subbatch_count; ++i) {
                blst::middle_product(
                    std::span<const blst::G1Point>(ctx.c1_points.data(), ctx.c1_points.size()),
                    precomputed_b,
                    std::span<blst::GT>(out.data(), out.size()),
                    workspace);
            }
            benchmark::DoNotOptimize(out.data());
            benchmark::ClobberMemory();
        }

        set_beat_cuckoo_counters(state, params);
    }

    static void BM_BtxVsBeat_Offline_BeatPlusPlus40(benchmark::State& state)
    {
        const BeatCuckooParameters params = select_parameters(
            std::span<const BeatCuckooParameters>(
                kBeat40BitParameters.data(),
                kBeat40BitParameters.size()),
            static_cast<std::size_t>(state.range(0)));
        run_beat_plus_plus_precompute(state, params, 0xbea80000ULL);
    }

    static void BM_BtxVsBeat_Offline_BeatPlusPlus128(benchmark::State& state)
    {
        const BeatCuckooParameters params = select_parameters(
            std::span<const BeatCuckooParameters>(
                kBeat128BitParameters.data(),
                kBeat128BitParameters.size()),
            static_cast<std::size_t>(state.range(0)));
        run_beat_plus_plus_precompute(state, params, 0xbea81000ULL);
    }

    void run_beat_plus_plus_ctxt_check(
        benchmark::State& state,
        const BeatCuckooParameters& params,
        const std::uint64_t seed)
    {
        const std::size_t checked_ciphertexts = params.batch_size * params.copies;
        const std::vector<BeatClientBatchCheck> checks =
            make_beat_client_batch_checks(params.batch_size, params.copies, seed + params.batch_size);
        for (const BeatClientBatchCheck& check : checks) {
            for (const BeatClientCiphertextCheck& copy : check.checks) {
                if (!verify_beat_client_ciphertext_check(copy)) {
                    throw std::runtime_error("Generated invalid BEAT++ client DLEq check fixture");
                }
            }
            if (!verify_beat_client_batch_check(check)) {
                throw std::runtime_error("Generated invalid BEAT++ batched client DLEq check fixture");
            }
        }

        for (auto _ : state) {
            bool ok = true;
            for (const BeatClientBatchCheck& check : checks) {
                ok = verify_beat_client_batch_check(check) && ok;
            }
            benchmark::DoNotOptimize(ok);
            benchmark::ClobberMemory();
        }

        state.counters["checked_ciphertexts"] = static_cast<double>(checked_ciphertexts);
        state.counters["client_batches"] = static_cast<double>(params.batch_size);
        set_beat_cuckoo_counters(state, params);
    }

    static void BM_BtxVsBeat_CtxtCheck_BeatPlusPlus40(benchmark::State& state)
    {
        const BeatCuckooParameters params = select_parameters(
            std::span<const BeatCuckooParameters>(
                kBeat40BitParameters.data(),
                kBeat40BitParameters.size()),
            static_cast<std::size_t>(state.range(0)));
        run_beat_plus_plus_ctxt_check(state, params, 0xbea82000ULL);
    }

    static void BM_BtxVsBeat_CtxtCheck_BeatPlusPlus128(benchmark::State& state)
    {
        const BeatCuckooParameters params = select_parameters(
            std::span<const BeatCuckooParameters>(
                kBeat128BitParameters.data(),
                kBeat128BitParameters.size()),
            static_cast<std::size_t>(state.range(0)));
        run_beat_plus_plus_ctxt_check(state, params, 0xbea83000ULL);
    }

    void run_beat_plus_plus_combine(
        benchmark::State& state,
        const BeatCuckooParameters& params,
        const std::size_t threshold,
        const std::uint64_t seed)
    {
        const std::size_t server_count = std::max(kDefaultServerCount, threshold);
        HeadToHeadContext ctx(
            params.subbatch_size,
            threshold,
            server_count,
            seed + (params.batch_size << 8) + threshold);

        for (auto _ : state) {
            blst::G1Point accumulator = blst::G1Point::infinity();
            for (std::size_t i = 0; i < params.subbatch_count; ++i) {
                const blst::G1Point combined = CombinePartialDecryptionShares(
                    std::span<const PartialDecryptionShare>(
                        ctx.threshold_shares.data(),
                        ctx.threshold_shares.size()),
                    ctx.setup.public_parameters);
                accumulator = add_g1(accumulator, combined);
            }
            benchmark::DoNotOptimize(accumulator);
            benchmark::ClobberMemory();
        }

        set_beat_cuckoo_threshold_counters(state, params, threshold);
    }

    static void BM_BtxVsBeat_Combine_BeatPlusPlus40(benchmark::State& state)
    {
        const BeatCuckooParameters params = select_parameters(
            std::span<const BeatCuckooParameters>(
                kBeat40BitParameters.data(),
                kBeat40BitParameters.size()),
            static_cast<std::size_t>(state.range(0)));
        run_beat_plus_plus_combine(
            state,
            params,
            static_cast<std::size_t>(state.range(1)),
            0xbea84000ULL);
    }

    static void BM_BtxVsBeat_Combine_BeatPlusPlus128(benchmark::State& state)
    {
        const BeatCuckooParameters params = select_parameters(
            std::span<const BeatCuckooParameters>(
                kBeat128BitParameters.data(),
                kBeat128BitParameters.size()),
            static_cast<std::size_t>(state.range(0)));
        run_beat_plus_plus_combine(
            state,
            params,
            static_cast<std::size_t>(state.range(1)),
            0xbea85000ULL);
    }

    struct ServerCheckInputs
    {
        std::vector<blst::G1AffinePoint> g1_affine;
        std::vector<blst::G2AffinePoint> g2_affine;
        std::vector<const blst::G1AffinePoint*> g1_ptrs;
        std::vector<const blst::G2AffinePoint*> g2_ptrs;
    };

    ServerCheckInputs make_server_check_inputs(
        const HeadToHeadContext& ctx,
        const std::size_t repeat_count = 1)
    {
        ServerCheckInputs inputs{};
        const std::size_t pair_count = ctx.c1_points.size();
        inputs.g1_affine.resize(pair_count);
        inputs.g2_affine.resize(pair_count);
        inputs.g1_ptrs.reserve(pair_count * repeat_count);
        inputs.g2_ptrs.reserve(pair_count * repeat_count);

        for (std::size_t i = 0; i < pair_count; ++i) {
            blst::G1Point::to_affine(inputs.g1_affine[i], ctx.c1_points[i]);
            blst::G2Point::to_affine(inputs.g2_affine[i], ctx.setup.public_parameters.g2_powers[i + 1]);
        }

        for (std::size_t batch = 0; batch < repeat_count; ++batch) {
            for (std::size_t i = 0; i < pair_count; ++i) {
                inputs.g1_ptrs.push_back(&inputs.g1_affine[i]);
                inputs.g2_ptrs.push_back(&inputs.g2_affine[i]);
            }
        }

        return inputs;
    }

    blst::GT run_server_check_multipairing(const ServerCheckInputs& inputs)
    {
        blst::GT miller{};
        blst::GT::miller_loop_n(
            miller,
            inputs.g2_ptrs.data(),
            inputs.g1_ptrs.data(),
            inputs.g1_ptrs.size());
        blst::GT out{};
        blst::GT::final_mult(out, miller);
        return out;
    }

    void run_beat_plus_plus_server_check(
        benchmark::State& state,
        const BeatCuckooParameters& params,
        const std::uint64_t seed)
    {
        HeadToHeadContext ctx(
            params.subbatch_size,
            kDefaultThreshold,
            kDefaultServerCount,
            seed + params.batch_size);
        const ServerCheckInputs inputs = make_server_check_inputs(ctx, params.subbatch_count);

        for (auto _ : state) {
            benchmark::DoNotOptimize(run_server_check_multipairing(inputs));
            benchmark::ClobberMemory();
        }

        state.counters["checked_pairs"] =
            static_cast<double>(params.subbatch_count * params.subbatch_size);
        state.counters["server_check_batches"] = 1.0;
        set_beat_cuckoo_counters(state, params);
    }

    static void BM_BtxVsBeat_ServerCheck_BeatPlusPlus40(benchmark::State& state)
    {
        const BeatCuckooParameters params = select_parameters(
            std::span<const BeatCuckooParameters>(
                kBeat40BitParameters.data(),
                kBeat40BitParameters.size()),
            static_cast<std::size_t>(state.range(0)));
        run_beat_plus_plus_server_check(state, params, 0xbea86000ULL);
    }

    static void BM_BtxVsBeat_ServerCheck_BeatPlusPlus128(benchmark::State& state)
    {
        const BeatCuckooParameters params = select_parameters(
            std::span<const BeatCuckooParameters>(
                kBeat128BitParameters.data(),
                kBeat128BitParameters.size()),
            static_cast<std::size_t>(state.range(0)));
        run_beat_plus_plus_server_check(state, params, 0xbea87000ULL);
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

    void run_beat_plus_plus_open(
        benchmark::State& state,
        const BeatCuckooParameters& params,
        const std::uint64_t seed)
    {
        HeadToHeadContext ctx(
            params.subbatch_size,
            kDefaultThreshold,
            kDefaultServerCount,
            seed + params.batch_size);

        std::vector<blst::GT> out(params.subbatch_size, blst::GT::identity());
        const blst::MiddleProductPrecomputedB& precomputed_b =
            ctx.prepared_public_parameters.alpha_precomputed_b_by_batch_size[params.subbatch_size];

        for (auto _ : state) {
            for (std::size_t i = 0; i < params.subbatch_count; ++i) {
                blst::update(
                    ctx.aggregate_predecryption_key,
                    precomputed_b,
                    std::span<blst::GT>(out.data(), out.size()));
            }
            benchmark::DoNotOptimize(out.data());
            benchmark::ClobberMemory();
        }

        set_beat_cuckoo_counters(state, params);
    }

    static void BM_BtxVsBeat_Online_BeatPlusPlus40(benchmark::State& state)
    {
        const BeatCuckooParameters params = select_parameters(
            std::span<const BeatCuckooParameters>(
                kBeat40BitParameters.data(),
                kBeat40BitParameters.size()),
            static_cast<std::size_t>(state.range(0)));
        run_beat_plus_plus_open(state, params, 0xbeaa0000ULL);
    }

    static void BM_BtxVsBeat_Online_BeatPlusPlus128(benchmark::State& state)
    {
        const BeatCuckooParameters params = select_parameters(
            std::span<const BeatCuckooParameters>(
                kBeat128BitParameters.data(),
                kBeat128BitParameters.size()),
            static_cast<std::size_t>(state.range(0)));
        run_beat_plus_plus_open(state, params, 0xbeab0000ULL);
    }
}

BENCHMARK(BM_BtxVsBeat_Offline_BtxPrepared)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_Offline_BeatPlusPlus40)->Apply(apply_paper_batch_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_Offline_BeatPlusPlus128)->Apply(apply_paper_batch_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_CtxtCheck_BeatPlusPlus40)->Apply(apply_paper_batch_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_CtxtCheck_BeatPlusPlus128)->Apply(apply_paper_batch_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_Combine_BeatPlusPlus40)->Apply(apply_paper_batch_threshold_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_Combine_BeatPlusPlus128)->Apply(apply_paper_batch_threshold_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_ServerCheck_BeatPlusPlus40)->Apply(apply_paper_batch_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_ServerCheck_BeatPlusPlus128)->Apply(apply_paper_batch_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_Online_BtxPrepared)->Apply(apply_batch_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_Online_BeatPlusPlus40)->Apply(apply_paper_batch_args)->UseRealTime();
BENCHMARK(BM_BtxVsBeat_Online_BeatPlusPlus128)->Apply(apply_paper_batch_args)->UseRealTime();

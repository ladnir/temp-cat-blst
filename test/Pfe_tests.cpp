#include <array>
#include <cstdint>
#include <cstring>
#include <span>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "BTE/building_blocks/hkdf_sha256.h"
#include "BTE/pfe/core.h"
#include "InsecureRng.h"

using namespace bte::pfe;

namespace
{
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

    blst::Fr fr_sub(const blst::Fr& lhs, const blst::Fr& rhs)
    {
        blst::Fr out{};
        blst::Fr::sub(out, lhs, rhs);
        return out;
    }

    blst::Fr fr_mul(const blst::Fr& lhs, const blst::Fr& rhs)
    {
        blst::Fr out{};
        blst::Fr::mul(out, lhs, rhs);
        return out;
    }

    blst::Fr fr_inverse(const blst::Fr& value)
    {
        if (fr_is_zero(value)) {
            throw std::invalid_argument("Cannot invert zero");
        }

        blst::Fr out{};
        blst::Fr::reciprocal(out, value);
        return out;
    }

    blst::GT add_gt(const blst::GT& lhs, const blst::GT& rhs)
    {
        blst::GT out{};
        blst::GT::add(out, lhs, rhs);
        return out;
    }

    blst::GT sub_gt(const blst::GT& lhs, const blst::GT& rhs)
    {
        blst::GT neg_rhs{};
        blst::GT::neg(neg_rhs, rhs);
        return add_gt(lhs, neg_rhs);
    }

    blst::G1Point add_g1(const blst::G1Point& lhs, const blst::G1Point& rhs)
    {
        blst::G1Point out{};
        blst::G1Point::dadd(out, lhs, rhs);
        return out;
    }

    blst::G1Point sub_g1(const blst::G1Point& lhs, const blst::G1Point& rhs)
    {
        blst::G1Point neg_rhs{};
        blst::G1Point::neg(neg_rhs, rhs);
        return add_g1(lhs, neg_rhs);
    }

    blst::G2Point add_g2(const blst::G2Point& lhs, const blst::G2Point& rhs)
    {
        blst::G2Point out{};
        blst::G2Point::dadd(out, lhs, rhs);
        return out;
    }

    blst::G1Point scale_g1(const blst::G1Point& point, const blst::Fr& scalar)
    {
        if (fr_is_zero(scalar)) {
            return blst::G1Point::infinity();
        }

        std::array<std::uint8_t, 32> scalar_bytes{};
        blst::Fr::to_scaler(scalar_bytes.data(), scalar);
        blst::G1Point out{};
        blst::G1Point::mult(out, point, scalar_bytes.data(), 256);
        return out;
    }

    blst::G2Point scale_g2(const blst::G2Point& point, const blst::Fr& scalar)
    {
        if (fr_is_zero(scalar)) {
            return blst::G2Point::infinity();
        }

        std::array<std::uint8_t, 32> scalar_bytes{};
        blst::Fr::to_scaler(scalar_bytes.data(), scalar);
        blst::G2Point out{};
        blst::G2Point::mult(out, point, scalar_bytes.data(), 256);
        return out;
    }

    blst::GT scale_gt(const blst::GT& value, const blst::Fr& scalar)
    {
        if (fr_is_zero(scalar)) {
            return blst::GT::identity();
        }

        std::array<std::uint8_t, 32> scalar_bytes{};
        blst::Fr::to_scaler(scalar_bytes.data(), scalar);
        blst::GT out{};
        blst::GT::mult(out, value, scalar_bytes.data(), 256);
        return out;
    }

    blst::GT pairing(const blst::G1Point& p, const blst::G2Point& q)
    {
        blst::G1AffinePoint p_aff{};
        blst::G2AffinePoint q_aff{};
        blst::G1Point::to_affine(p_aff, p);
        blst::G2Point::to_affine(q_aff, q);

        blst::MillerResult miller{};
        blst::GT::miller_loop(miller, &q_aff, &p_aff);
        blst::GT out{};
        blst::GT::final_exp(out, miller);
        return out;
    }

    Bytes kdf_from_gt(const PfeParameters& params, const blst::GT& value)
    {
        std::array<std::uint8_t, 12 * 48> bytes{};
        blst::GT::serialize(bytes.data(), value);
        const Bytes ikm(bytes.begin(), bytes.end());
        const std::span<const std::uint8_t> salt(
            reinterpret_cast<const std::uint8_t*>(params.kdf_salt.data()),
            params.kdf_salt.size());
        const std::span<const std::uint8_t> info(
            reinterpret_cast<const std::uint8_t*>(params.kdf_info.data()),
            params.kdf_info.size());
        return bte::building_blocks::hkdf_sha256(salt, ikm, info, params.key_size);
    }

    blst::Fr slot_difference_inverse(const std::size_t slot, const std::size_t other_slot)
    {
        return fr_inverse(fr_sub(fr_from_u64(slot), fr_from_u64(other_slot)));
    }

    std::vector<Bytes> derive_keys_reference(
        const PfeParameters& params,
        const PublicParameters& public_parameters,
        const std::span<const Ciphertext> ciphertexts,
        const blst::G1Point& aggregate_predecryption_key)
    {
        std::vector<Bytes> keys;
        keys.reserve(ciphertexts.size());

        blst::G2Point sum_dk = blst::G2Point::infinity();
        for (const blst::G2Point& dk : public_parameters.dk) {
            sum_dk = add_g2(sum_dk, dk);
        }

        std::vector<blst::GT> tm_self(ciphertexts.size(), blst::GT::identity());
        for (std::size_t i = 0; i < ciphertexts.size(); ++i) {
            tm_self[i] = pairing(ciphertexts[i].c1, public_parameters.dk[i]);
        }

        for (std::size_t j = 0; j < ciphertexts.size(); ++j) {
            const std::size_t slot = j + 1;

            blst::Fr self_coefficient = fr_add(
                fr_inverse(fr_from_u64(slot)),
                fr_inverse(fr_add(fr_from_u64(slot), blst::Fr::one())));

            blst::G1Point weighted_c1 = blst::G1Point::infinity();
            blst::G2Point weighted_dk = blst::G2Point::infinity();
            blst::GT weighted_tm = blst::GT::identity();
            for (std::size_t i = 0; i < ciphertexts.size(); ++i) {
                if (i == j) {
                    continue;
                }

                const blst::Fr coeff = slot_difference_inverse(slot, i + 1);
                self_coefficient = fr_add(self_coefficient, coeff);
                weighted_c1 = add_g1(weighted_c1, scale_g1(ciphertexts[i].c1, coeff));
                weighted_dk = add_g2(weighted_dk, scale_g2(public_parameters.dk[i], coeff));
                weighted_tm = add_gt(weighted_tm, scale_gt(tm_self[i], coeff));
            }

            blst::G1Point l1_input = sub_g1(ciphertexts[j].c2, aggregate_predecryption_key);
            l1_input = add_g1(l1_input, scale_g1(ciphertexts[j].c1, self_coefficient));
            l1_input = sub_g1(l1_input, weighted_c1);

            blst::GT l1 = pairing(l1_input, public_parameters.dk[j]);
            l1 = sub_gt(l1, pairing(ciphertexts[j].c1, weighted_dk));
            l1 = add_gt(l1, weighted_tm);

            const blst::GT l2 = sub_gt(
                pairing(ciphertexts[j].c2, blst::G2Point::generator()),
                pairing(ciphertexts[j].c1, sum_dk));

            const blst::GT rp0 = scale_gt(
                sub_gt(scale_gt(l1, fr_from_u64(slot + 1)), l2),
                fr_from_u64(slot));
            keys.push_back(kdf_from_gt(params, rp0));
        }

        return keys;
    }

    struct PfeFixture
    {
        PfeParameters params;
        InsecureRng rng;
        SetupResult setup;

        explicit PfeFixture(const PfeParameters& custom_params = default_params())
            : params(custom_params),
              rng(0x55aa33cc77ee11ddULL),
              setup(Setup(params, [this](std::span<std::uint8_t> out) { fill_random(out); }))
        {
        }

        static PfeParameters default_params()
        {
            return PfeParameters{4, 4, 2, 32, "BTE_PFE_KDF_SALT_V1", "BTE_PFE_KDF_INFO_V1"};
        }

        void fill_random(const std::span<std::uint8_t> out)
        {
            std::uint64_t word = 0;
            std::size_t remaining = 0;
            for (std::uint8_t& byte : out) {
                if (remaining == 0) {
                    word = rng.next_u64();
                    remaining = sizeof(word);
                }

                byte = static_cast<std::uint8_t>(word & 0xffU);
                word >>= 8;
                --remaining;
            }
        }

        std::vector<EncryptionResult> encrypt_batch()
        {
            std::vector<EncryptionResult> encrypted;
            encrypted.reserve(params.batch_size);
            for (std::size_t i = 0; i < params.batch_size; ++i) {
                encrypted.push_back(Encrypt(
                    params,
                    setup.public_parameters,
                    [this](std::span<std::uint8_t> out) { fill_random(out); }));
            }
            return encrypted;
        }

        std::vector<Ciphertext> extract_ciphertexts(
            const std::span<const EncryptionResult> encrypted) const
        {
            std::vector<Ciphertext> ciphertexts;
            ciphertexts.reserve(encrypted.size());
            for (const EncryptionResult& item : encrypted) {
                ciphertexts.push_back(item.ciphertext);
            }
            return ciphertexts;
        }

        std::vector<PartialDecryptionShare> create_shares(
            const ThresholdRequest& request,
            const std::span<const std::size_t> secret_share_indices) const
        {
            std::vector<PartialDecryptionShare> shares;
            shares.reserve(secret_share_indices.size());
            for (const std::size_t index : secret_share_indices) {
                shares.push_back(CreatePartialDecryptionShare(setup.secret_shares[index], request));
            }
            return shares;
        }
    };
}

TEST(Pfe, SetupPublishesExpectedLayout)
{
    PfeFixture fixture;

    EXPECT_EQ(fixture.setup.public_parameters.dk.size(), fixture.params.batch_size);
    EXPECT_EQ(
        fixture.setup.public_parameters.weighted_dk_excluding_slot.size(),
        fixture.params.batch_size);
    EXPECT_EQ(
        fixture.setup.public_parameters.self_coefficients.size(),
        fixture.params.batch_size);
    EXPECT_EQ(fixture.setup.secret_shares.size(), fixture.params.num_servers);
    EXPECT_EQ(
        fixture.setup.secret_shares.front().p_shares.size(),
        fixture.params.batch_size);
    EXPECT_EQ(
        fixture.setup.public_parameters.toeplitz_fft_size,
        static_cast<std::size_t>(8));
}

TEST(Pfe, EndToEndThresholdFlowRecoversDerivedKeys)
{
    PfeFixture fixture;
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch();
    const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);

    const ThresholdRequest request =
        BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts);
    const std::array<std::size_t, 2> subset = {0, 1};
    const std::vector<PartialDecryptionShare> shares = fixture.create_shares(request, subset);
    const blst::G1Point sbk = CombinePartialDecryptionShares(
        std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
        fixture.setup.public_parameters);

    const BatchDecryptionPrecomputation precomputation = PrecomputeBatchDecryption(
        fixture.setup.public_parameters,
        ciphertexts);
    const std::vector<Bytes> recovered = FinalizeBatchDecryption(
        fixture.params,
        fixture.setup.public_parameters,
        precomputation,
        sbk);

    ASSERT_EQ(recovered.size(), encrypted.size());
    for (std::size_t i = 0; i < encrypted.size(); ++i) {
        EXPECT_EQ(recovered[i], encrypted[i].key) << "key mismatch at slot " << i;
    }
}

TEST(Pfe, RecoverBatchKeysMatchesPrecomputedFinalize)
{
    PfeFixture fixture;
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch();
    const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);

    const ThresholdRequest request =
        BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts);
    const std::array<std::size_t, 2> subset = {1, 3};
    const std::vector<PartialDecryptionShare> shares = fixture.create_shares(request, subset);
    const blst::G1Point sbk = CombinePartialDecryptionShares(
        std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
        fixture.setup.public_parameters);

    const BatchDecryptionPrecomputation precomputation = PrecomputeBatchDecryption(
        fixture.setup.public_parameters,
        ciphertexts);
    const std::vector<Bytes> finalized = FinalizeBatchDecryption(
        fixture.params,
        fixture.setup.public_parameters,
        precomputation,
        sbk);
    const std::vector<Bytes> recovered = RecoverBatchKeys(
        fixture.params,
        fixture.setup.public_parameters,
        ciphertexts,
        sbk);

    EXPECT_EQ(recovered, finalized);
}

TEST(Pfe, PreparedPathMatchesDirectPath)
{
    PfeFixture fixture(PfeParameters{6, 6, 3, 32, "BTE_PFE_KDF_SALT_V1", "BTE_PFE_KDF_INFO_V1"});
    const PreparedPublicParameters prepared =
        PreparePublicParameters(fixture.setup.public_parameters);
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch();
    const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);

    const ThresholdRequest request =
        BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts);
    const std::array<std::size_t, 3> subset = {0, 2, 5};
    const std::vector<PartialDecryptionShare> shares = fixture.create_shares(request, subset);
    const blst::G1Point sbk = CombinePartialDecryptionShares(
        std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
        fixture.setup.public_parameters);

    const BatchDecryptionPrecomputation direct_precomputation = PrecomputeBatchDecryption(
        fixture.setup.public_parameters,
        ciphertexts);
    const BatchDecryptionPrecomputation prepared_precomputation = PrecomputeBatchDecryption(
        prepared,
        ciphertexts);

    EXPECT_EQ(
        FinalizeBatchDecryption(
            fixture.params,
            fixture.setup.public_parameters,
            direct_precomputation,
            sbk),
        FinalizeBatchDecryption(
            fixture.params,
            prepared,
            prepared_precomputation,
            sbk));

    EXPECT_EQ(
        RecoverBatchKeys(
            fixture.params,
            fixture.setup.public_parameters,
            ciphertexts,
            sbk),
        RecoverBatchKeys(
            fixture.params,
            prepared,
            ciphertexts,
            sbk));
}

TEST(Pfe, OptimizedFinalizeMatchesSlowReference)
{
    const std::array<PfeParameters, 4> parameter_sets = {
        PfeParameters{2, 3, 2, 32, "BTE_PFE_KDF_SALT_V1", "BTE_PFE_KDF_INFO_V1"},
        PfeParameters{4, 4, 2, 32, "BTE_PFE_KDF_SALT_V1", "BTE_PFE_KDF_INFO_V1"},
        PfeParameters{5, 5, 3, 32, "BTE_PFE_KDF_SALT_V1", "BTE_PFE_KDF_INFO_V1"},
        PfeParameters{6, 6, 4, 32, "BTE_PFE_KDF_SALT_V1", "BTE_PFE_KDF_INFO_V1"},
    };

    for (const PfeParameters& params : parameter_sets) {
        PfeFixture fixture(params);
        const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch();
        const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);
        const ThresholdRequest request =
            BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts);

        std::vector<std::size_t> subset;
        subset.reserve(params.threshold);
        for (std::size_t i = 0; i < params.threshold; ++i) {
            subset.push_back(i);
        }

        const std::vector<PartialDecryptionShare> shares = fixture.create_shares(request, subset);
        const blst::G1Point sbk = CombinePartialDecryptionShares(
            std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
            fixture.setup.public_parameters);

        const BatchDecryptionPrecomputation precomputation = PrecomputeBatchDecryption(
            fixture.setup.public_parameters,
            ciphertexts);
        const std::vector<Bytes> optimized = FinalizeBatchDecryption(
            fixture.params,
            fixture.setup.public_parameters,
            precomputation,
            sbk);
        const std::vector<Bytes> reference = derive_keys_reference(
            fixture.params,
            fixture.setup.public_parameters,
            ciphertexts,
            sbk);

        EXPECT_EQ(optimized, reference) << "batch_size=" << params.batch_size;
        EXPECT_EQ(optimized.size(), encrypted.size());
        for (std::size_t i = 0; i < optimized.size(); ++i) {
            EXPECT_EQ(optimized[i], encrypted[i].key)
                << "encrypt/decrypt mismatch at slot " << i
                << " for batch_size=" << params.batch_size;
        }
    }
}

TEST(Pfe, DifferentThresholdSubsetsRecoverSameKeys)
{
    PfeFixture fixture(PfeParameters{5, 5, 3, 32, "BTE_PFE_KDF_SALT_V1", "BTE_PFE_KDF_INFO_V1"});
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch();
    const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);
    const ThresholdRequest request =
        BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts);

    const std::array<std::array<std::size_t, 3>, 4> subsets = {{
        {0, 1, 2},
        {0, 2, 4},
        {1, 2, 3},
        {2, 3, 4},
    }};

    std::vector<Bytes> baseline;
    for (std::size_t subset_index = 0; subset_index < subsets.size(); ++subset_index) {
        const std::vector<PartialDecryptionShare> shares =
            fixture.create_shares(request, subsets[subset_index]);
        const blst::G1Point sbk = CombinePartialDecryptionShares(
            std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
            fixture.setup.public_parameters);

        const std::vector<Bytes> recovered = RecoverBatchKeys(
            fixture.params,
            fixture.setup.public_parameters,
            ciphertexts,
            sbk);

        if (subset_index == 0) {
            baseline = recovered;
        } else {
            EXPECT_EQ(recovered, baseline) << "subset index " << subset_index;
        }
    }
}

TEST(Pfe, RejectsNonFullBatch)
{
    PfeFixture fixture;
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch();
    std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);
    ciphertexts.pop_back();

    EXPECT_THROW(
        BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts),
        std::invalid_argument);
}

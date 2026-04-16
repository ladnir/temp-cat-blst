#include <array>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#if __has_include("BTE/beat/core.h")

#include "BTE/beat/core.h"
#include "BTE/building_blocks/lagrange.h"
#include "blst/curve/GT.h"
#include "InsecureRng.h"

using namespace bte::beat;

namespace
{
    std::array<std::uint8_t, 32> fr_to_scalar_bytes(const blst::Fr& scalar)
    {
        std::array<std::uint8_t, 32> bytes{};
        blst::Fr::to_scaler(bytes.data(), scalar);
        return bytes;
    }

    blst::G1Point scale_g1(const blst::G1Point& point, const blst::Fr& scalar)
    {
        const std::array<std::uint8_t, 32> scalar_bytes = fr_to_scalar_bytes(scalar);
        blst::G1Point out{};
        blst::G1Point::mult(out, point, scalar_bytes.data(), 256);
        return out;
    }

    blst::G1Point add_g1(const blst::G1Point& lhs, const blst::G1Point& rhs)
    {
        blst::G1Point out{};
        blst::G1Point::dadd(out, lhs, rhs);
        return out;
    }

    blst::G1Point naive_msm_sum_g1(
        const std::span<const blst::G1Point> points,
        const std::span<const blst::Fr> scalars)
    {
        blst::G1Point out = blst::G1Point::infinity();
        for (std::size_t i = 0; i < points.size(); ++i) {
            out = add_g1(out, scale_g1(points[i], scalars[i]));
        }
        return out;
    }

    struct BeatFixture
    {
        BeatParameters params;
        InsecureRng rng;
        SetupResult setup;
        PreparedPublicParameters prepared_public_parameters;

        explicit BeatFixture(
            BeatParameters params_in = BeatParameters{
                8,
                4,
                2,
                32,
                "BTE_BEAT_KDF_SALT_V1",
                "BTE_BEAT_KDF_INFO_V1"},
            const std::uint64_t seed = 0x1029384756abcdefULL)
            : params(std::move(params_in)),
              rng(seed),
              setup(Setup(params, [this](std::span<std::uint8_t> out) { fill_random(out); })),
              prepared_public_parameters(PreparePublicParameters(setup.public_parameters))
        {
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

        std::vector<EncryptionResult> encrypt_batch(const std::size_t batch_size)
        {
            std::vector<EncryptionResult> encrypted;
            encrypted.reserve(batch_size);
            for (std::size_t i = 0; i < batch_size; ++i) {
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
    };
}

TEST(Beat, SetupPublishesExpectedPowerLayout)
{
    BeatFixture fixture;

    ASSERT_EQ(
        fixture.setup.public_parameters.g2_powers.size(),
        (2 * fixture.params.max_batch_size) + 1);
    EXPECT_TRUE(blst::G2Point::is_inf(
        fixture.setup.public_parameters.g2_powers[fixture.params.max_batch_size + 1]));
    EXPECT_EQ(fixture.setup.secret_shares.size(), fixture.params.num_servers);
    EXPECT_EQ(
        fixture.setup.secret_shares.front().tau_power_shares.size(),
        fixture.params.max_batch_size);
    EXPECT_EQ(
        fixture.prepared_public_parameters.beta_precomputed_b_by_batch_size.size(),
        fixture.params.max_batch_size + 1);
    EXPECT_EQ(
        fixture.prepared_public_parameters.alpha_precomputed_b_by_batch_size.size(),
        fixture.params.max_batch_size + 1);
}

TEST(Beat, EndToEndThresholdFlowRecoversDerivedKeys)
{
    BeatFixture fixture;
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch(4);
    const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);

    const ThresholdRequest request =
        BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts);
    const std::array<PartialDecryptionShare, 2> shares = {
        CreatePartialDecryptionShare(fixture.setup.secret_shares[0], request),
        CreatePartialDecryptionShare(fixture.setup.secret_shares[1], request),
    };
    const blst::G1Point sbk = CombinePartialDecryptionShares(
        std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
        fixture.setup.public_parameters);

    const BatchDecryptionPrecomputation precomputation = PrecomputeBatchDecryption(
        fixture.prepared_public_parameters,
        ciphertexts);
    const std::vector<Bytes> recovered = FinalizeBatchDecryption(
        fixture.params,
        fixture.prepared_public_parameters,
        precomputation,
        sbk);

    ASSERT_EQ(recovered.size(), encrypted.size());
    for (std::size_t i = 0; i < encrypted.size(); ++i) {
        EXPECT_EQ(recovered[i], encrypted[i].key) << "key mismatch at slot " << i;
    }
}

TEST(Beat, PreparedFinalizeUsesCachedPathAndMatchesDirectFinalize)
{
    BeatFixture fixture;
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch(5);
    const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);

    const ThresholdRequest request =
        BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts);
    const std::array<PartialDecryptionShare, 2> shares = {
        CreatePartialDecryptionShare(fixture.setup.secret_shares[0], request),
        CreatePartialDecryptionShare(fixture.setup.secret_shares[2], request),
    };
    const blst::G1Point sbk = CombinePartialDecryptionShares(
        std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
        fixture.setup.public_parameters);

    const BatchDecryptionPrecomputation precomputation = PrecomputeBatchDecryption(
        fixture.prepared_public_parameters,
        ciphertexts);
    const std::vector<Bytes> direct = FinalizeBatchDecryption(
        fixture.params,
        fixture.setup.public_parameters,
        precomputation,
        sbk);
    const std::vector<Bytes> prepared = FinalizeBatchDecryption(
        fixture.params,
        fixture.prepared_public_parameters,
        precomputation,
        sbk);

    ASSERT_EQ(direct.size(), prepared.size());
    for (std::size_t i = 0; i < direct.size(); ++i) {
        EXPECT_EQ(prepared[i], direct[i]) << "prepared path mismatch at slot " << i;
        EXPECT_EQ(prepared[i], encrypted[i].key) << "prepared path key mismatch at slot " << i;
    }
}

TEST(Beat, CreatePartialDecryptionShareMatchesNaiveMsmAtEightSlotCutoff)
{
    BeatFixture fixture(
        BeatParameters{8, 16, 8, 32, "BTE_BEAT_KDF_SALT_V1", "BTE_BEAT_KDF_INFO_V1"},
        0x2233445566778899ULL);
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch(8);
    const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);

    const ThresholdRequest request =
        BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts);
    const PartialDecryptionShare share =
        CreatePartialDecryptionShare(fixture.setup.secret_shares[0], request);

    const blst::G1Point expected = naive_msm_sum_g1(
        std::span<const blst::G1Point>(request.ordered_c1.data(), request.ordered_c1.size()),
        std::span<const blst::Fr>(
            fixture.setup.secret_shares[0].tau_power_shares.data(),
            request.ordered_c1.size()));

    EXPECT_TRUE(blst::G1Point::equal(share.value, expected));
}

TEST(Beat, CombinePartialDecryptionSharesMatchesNaiveMsmAtEightShareCutoff)
{
    BeatFixture fixture(
        BeatParameters{8, 16, 8, 32, "BTE_BEAT_KDF_SALT_V1", "BTE_BEAT_KDF_INFO_V1"},
        0x33445566778899aaULL);
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch(8);
    const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);

    const ThresholdRequest request =
        BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts);

    std::vector<PartialDecryptionShare> shares;
    shares.reserve(fixture.params.threshold);
    for (std::size_t i = 0; i < fixture.params.threshold; ++i) {
        shares.push_back(CreatePartialDecryptionShare(fixture.setup.secret_shares[i], request));
    }

    const blst::G1Point combined = CombinePartialDecryptionShares(
        std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
        fixture.setup.public_parameters);

    std::vector<std::uint32_t> indices;
    indices.reserve(shares.size());
    std::vector<blst::G1Point> share_points;
    share_points.reserve(shares.size());
    for (const PartialDecryptionShare& share : shares) {
        indices.push_back(share.index);
        share_points.push_back(share.value);
    }
    const std::vector<blst::Fr> lambdas = bte::building_blocks::lagrange_at_zero(indices);
    const blst::G1Point expected = naive_msm_sum_g1(
        std::span<const blst::G1Point>(share_points.data(), share_points.size()),
        std::span<const blst::Fr>(lambdas.data(), lambdas.size()));

    EXPECT_TRUE(blst::G1Point::equal(combined, expected));
}

#else

TEST(Beat, ModuleNotAvailableYet)
{
    GTEST_SKIP() << "BTE/beat is not present yet; scaffolded tests will activate when the module lands.";
}

#endif

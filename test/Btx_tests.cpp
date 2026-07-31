#include <array>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "BTE/btx/core.h"
#include "BTE/btx/serialization.h"
#include "blst/curve/GT.h"
#include "InsecureRng.h"

using namespace bte::btx;

namespace
{
    struct BtxFixture
    {
        BtxParameters params;
        InsecureRng rng;
        SetupResult setup;
        PreparedPublicParameters prepared_public_parameters;

        BtxFixture()
            : params{8, 4, 2, 32, "BTE_BTX_KDF_SALT_V1", "BTE_BTX_KDF_INFO_V1"},
              rng(0x1029384756abcdefULL),
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

TEST(Btx, SetupPublishesExpectedPowerLayout)
{
    BtxFixture fixture;

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
        fixture.prepared_public_parameters.fft_precomputed_b_by_batch_size.size(),
        fixture.params.max_batch_size + 1);
    EXPECT_EQ(
        fixture.prepared_public_parameters.alpha_precomputed_b_by_batch_size.size(),
        fixture.params.max_batch_size + 1);
    for (std::size_t batch_size = 1;
         batch_size <= fixture.params.max_batch_size;
         ++batch_size) {
        EXPECT_EQ(
            fixture.prepared_public_parameters.fft_precomputed_b_by_batch_size[batch_size].n,
            0U);
        EXPECT_EQ(
            fixture.prepared_public_parameters.alpha_precomputed_b_by_batch_size[batch_size].n,
            0U);
    }
}

TEST(Btx, EndToEndNaiveThresholdFlowRecoversDerivedKeys)
{
    BtxFixture fixture;
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
        fixture.setup.public_parameters,
        ciphertexts,
        CrossTermAlgorithm::Naive);
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

TEST(Btx, MiddleProductCrossTermsMatchNaive)
{
    BtxFixture fixture;
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch(5);
    const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);

    const std::vector<blst::GT> naive = ComputeCrossTerms(
        fixture.setup.public_parameters,
        ciphertexts,
        CrossTermAlgorithm::Naive);
    const std::vector<blst::GT> fft = ComputeCrossTerms(
        fixture.prepared_public_parameters,
        ciphertexts,
        CrossTermAlgorithm::MiddleProductFft);
    const std::vector<blst::GT> multipair = ComputeCrossTerms(
        fixture.setup.public_parameters,
        ciphertexts,
        CrossTermAlgorithm::MiddleProductMultipair);

    ASSERT_EQ(naive.size(), fft.size());
    ASSERT_EQ(naive.size(), multipair.size());
    for (std::size_t i = 0; i < naive.size(); ++i) {
        EXPECT_TRUE(blst::GT::equal(naive[i], fft[i])) << "fft mismatch at slot " << i;
        EXPECT_TRUE(blst::GT::equal(naive[i], multipair[i])) << "multipair mismatch at slot " << i;
    }
}

TEST(Btx, VerifyAggregatePredecryptionKeyAcceptsCombinedShares)
{
    BtxFixture fixture;
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch(4);
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

    EXPECT_TRUE(VerifyAggregatePredecryptionKey(
        request,
        sbk,
        fixture.setup.public_parameters));
}

TEST(Btx, VerifyAggregatePredecryptionKeyRejectsTampering)
{
    BtxFixture fixture;
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

    blst::G1Point wrong_sbk{};
    blst::G1Point::dadd(wrong_sbk, sbk, blst::G1Point::generator());
    EXPECT_FALSE(VerifyAggregatePredecryptionKey(
        request,
        wrong_sbk,
        fixture.setup.public_parameters));

    ThresholdRequest wrong_request = request;
    blst::G1Point wrong_c1{};
    blst::G1Point::dadd(wrong_c1, wrong_request.ordered_c1[0], blst::G1Point::generator());
    wrong_request.ordered_c1[0] = wrong_c1;
    EXPECT_FALSE(VerifyAggregatePredecryptionKey(
        wrong_request,
        sbk,
        fixture.setup.public_parameters));
}

TEST(Btx, RecoverBatchKeysWithMiddleProductMatchesEncryptorKeys)
{
    BtxFixture fixture;
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch(6);
    const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);

    const ThresholdRequest request =
        BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts);
    const std::array<PartialDecryptionShare, 2> shares = {
        CreatePartialDecryptionShare(fixture.setup.secret_shares[1], request),
        CreatePartialDecryptionShare(fixture.setup.secret_shares[3], request),
    };
    const blst::G1Point sbk = CombinePartialDecryptionShares(
        std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
        fixture.setup.public_parameters);

    const std::vector<Bytes> recovered = RecoverBatchKeys(
        fixture.params,
        fixture.prepared_public_parameters,
        ciphertexts,
        sbk,
        CrossTermAlgorithm::MiddleProductFft);

    ASSERT_EQ(recovered.size(), encrypted.size());
    for (std::size_t i = 0; i < encrypted.size(); ++i) {
        EXPECT_EQ(recovered[i], encrypted[i].key) << "fft middle-product mismatch at slot " << i;
    }
}

TEST(Btx, DefaultCrossTermPathUsesFftAndRecoversKeys)
{
    BtxFixture fixture;
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

    const BatchDecryptionPrecomputation precomputation =
        PrecomputeBatchDecryption(fixture.prepared_public_parameters, ciphertexts);
    EXPECT_EQ(precomputation.algorithm, CrossTermAlgorithm::MiddleProductFft);

    const std::vector<Bytes> recovered = RecoverBatchKeys(
        fixture.params,
        fixture.prepared_public_parameters,
        ciphertexts,
        sbk);

    ASSERT_EQ(recovered.size(), encrypted.size());
    for (std::size_t i = 0; i < encrypted.size(); ++i) {
        EXPECT_EQ(recovered[i], encrypted[i].key) << "default fft mismatch at slot " << i;
    }
}

TEST(Btx, PreparedFinalizeUsesCachedAlphaPathAndMatchesDirectFinalize)
{
    BtxFixture fixture;
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
        ciphertexts,
        CrossTermAlgorithm::MiddleProductFft);
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
        EXPECT_EQ(prepared[i], direct[i]) << "prepared alpha mismatch at slot " << i;
        EXPECT_EQ(prepared[i], encrypted[i].key) << "prepared alpha key mismatch at slot " << i;
    }
}

TEST(Btx, PreparedCachesPopulateOnlyForRequestedBatchSize)
{
    BtxFixture fixture;
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch(3);
    const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);

    const BatchDecryptionPrecomputation precomputation = PrecomputeBatchDecryption(
        fixture.prepared_public_parameters,
        ciphertexts,
        CrossTermAlgorithm::MiddleProductFft);
    EXPECT_EQ(
        fixture.prepared_public_parameters.fft_precomputed_b_by_batch_size[3].n,
        3U);
    EXPECT_EQ(
        fixture.prepared_public_parameters.fft_precomputed_b_by_batch_size[4].n,
        0U);

    const ThresholdRequest request =
        BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts);
    const std::array<PartialDecryptionShare, 2> shares = {
        CreatePartialDecryptionShare(fixture.setup.secret_shares[0], request),
        CreatePartialDecryptionShare(fixture.setup.secret_shares[1], request),
    };
    const blst::G1Point sbk = CombinePartialDecryptionShares(
        std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
        fixture.setup.public_parameters);

    const std::vector<Bytes> recovered = FinalizeBatchDecryption(
        fixture.params,
        fixture.prepared_public_parameters,
        precomputation,
        sbk);
    EXPECT_EQ(
        fixture.prepared_public_parameters.alpha_precomputed_b_by_batch_size[3].n,
        3U);
    EXPECT_EQ(
        fixture.prepared_public_parameters.alpha_precomputed_b_by_batch_size[4].n,
        0U);

    ASSERT_EQ(recovered.size(), encrypted.size());
    for (std::size_t i = 0; i < recovered.size(); ++i) {
        EXPECT_EQ(recovered[i], encrypted[i].key) << "prepared lazy-cache mismatch at slot " << i;
    }
}

TEST(Btx, CombineSharesRejectsInsufficientThreshold)
{
    BtxFixture fixture;
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch(3);
    const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);

    const ThresholdRequest request =
        BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts);
    const PartialDecryptionShare share =
        CreatePartialDecryptionShare(fixture.setup.secret_shares[0], request);
    const std::array<PartialDecryptionShare, 1> insufficient = { share };

    EXPECT_THROW(
        CombinePartialDecryptionShares(
            std::span<const PartialDecryptionShare>(insufficient.data(), insufficient.size()),
            fixture.setup.public_parameters),
        std::invalid_argument);
}

TEST(Btx, SerializationRoundTripPreservesObjects)
{
    BtxFixture fixture;
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch(3);
    const Ciphertext ciphertext = encrypted.front().ciphertext;

    const bte::building_blocks::Bytes encoded_public =
        SerializePublicParameters(fixture.setup.public_parameters);
    const PublicParameters decoded_public = DeserializePublicParameters(encoded_public);
    EXPECT_EQ(decoded_public.max_batch_size, fixture.setup.public_parameters.max_batch_size);
    EXPECT_EQ(decoded_public.num_servers, fixture.setup.public_parameters.num_servers);
    EXPECT_EQ(decoded_public.threshold, fixture.setup.public_parameters.threshold);
    ASSERT_EQ(decoded_public.g2_powers.size(), fixture.setup.public_parameters.g2_powers.size());
    EXPECT_TRUE(blst::GT::equal(
        decoded_public.gt_tau_to_bmax_plus_one,
        fixture.setup.public_parameters.gt_tau_to_bmax_plus_one));

    const bte::building_blocks::Bytes encoded_ciphertext = SerializeCiphertext(ciphertext);
    const Ciphertext decoded_ciphertext = DeserializeCiphertext(encoded_ciphertext);
    EXPECT_TRUE(blst::G1Point::equal(decoded_ciphertext.c1, ciphertext.c1));

    const std::vector<Ciphertext> ciphertexts = fixture.extract_ciphertexts(encrypted);
    const ThresholdRequest request =
        BuildThresholdRequest(fixture.setup.public_parameters, ciphertexts);
    const bte::building_blocks::Bytes encoded_request = SerializeThresholdRequest(request);
    const ThresholdRequest decoded_request = DeserializeThresholdRequest(encoded_request);
    ASSERT_EQ(decoded_request.ordered_c1.size(), request.ordered_c1.size());
    for (std::size_t i = 0; i < request.ordered_c1.size(); ++i) {
        EXPECT_TRUE(blst::G1Point::equal(decoded_request.ordered_c1[i], request.ordered_c1[i]));
    }

    const PartialDecryptionShare share =
        CreatePartialDecryptionShare(fixture.setup.secret_shares[0], request);
    const bte::building_blocks::Bytes encoded_share = SerializePartialDecryptionShare(share);
    const PartialDecryptionShare decoded_share = DeserializePartialDecryptionShare(encoded_share);
    EXPECT_EQ(decoded_share.index, share.index);
    EXPECT_TRUE(blst::G1Point::equal(decoded_share.value, share.value));
}

TEST(Btx, EndToEndFlowWorksAcrossSerializationBoundaries)
{
    BtxFixture fixture;
    const std::vector<EncryptionResult> encrypted = fixture.encrypt_batch(4);

    const PublicParameters public_parameters =
        DeserializePublicParameters(SerializePublicParameters(fixture.setup.public_parameters));
    const PreparedPublicParameters prepared_public_parameters =
        PreparePublicParameters(public_parameters);

    std::vector<Ciphertext> ciphertexts;
    ciphertexts.reserve(encrypted.size());
    for (const EncryptionResult& item : encrypted) {
        ciphertexts.push_back(DeserializeCiphertext(SerializeCiphertext(item.ciphertext)));
    }

    const ThresholdRequest request =
        DeserializeThresholdRequest(SerializeThresholdRequest(
            BuildThresholdRequest(public_parameters, ciphertexts)));

    const std::array<PartialDecryptionShare, 2> shares = {
        DeserializePartialDecryptionShare(SerializePartialDecryptionShare(
            CreatePartialDecryptionShare(fixture.setup.secret_shares[0], request))),
        DeserializePartialDecryptionShare(SerializePartialDecryptionShare(
            CreatePartialDecryptionShare(fixture.setup.secret_shares[2], request))),
    };

    const blst::G1Point sbk = CombinePartialDecryptionShares(
        std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
        public_parameters);

    const std::vector<Bytes> recovered = RecoverBatchKeys(
        fixture.params,
        prepared_public_parameters,
        ciphertexts,
        sbk,
        CrossTermAlgorithm::MiddleProductAuto);

    ASSERT_EQ(recovered.size(), encrypted.size());
    for (std::size_t i = 0; i < encrypted.size(); ++i) {
        EXPECT_EQ(recovered[i], encrypted[i].key) << "serialized flow mismatch at slot " << i;
    }
}

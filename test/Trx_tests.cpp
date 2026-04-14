#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "BTE/building_blocks/lagrange.h"
#include "BTE/building_blocks/blst_cpp.h"
#include "BTE/building_blocks/kzg.h"
#include "BTE/building_blocks/polynomial.h"
#include "BTE/trx/core.h"
#include "InsecureRng.h"

using namespace bte::trx;

namespace
{
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

    blst::Fr fr_neg(const blst::Fr& value)
    {
        blst::Fr out{};
        blst::Fr::neg(out, value);
        return out;
    }

    std::vector<blst::Fr> lagrange_at_zero_naive(const std::span<const std::uint32_t> indices)
    {
        std::vector<blst::Fr> lambdas(indices.size(), blst::Fr::one());
        for (std::size_t i = 0; i < indices.size(); ++i) {
            blst::Fr numerator = blst::Fr::one();
            blst::Fr denominator = blst::Fr::one();
            const blst::Fr i_fr = fr_from_u64(indices[i]);
            for (std::size_t j = 0; j < indices.size(); ++j) {
                if (i == j) {
                    continue;
                }

                numerator = fr_mul(numerator, fr_neg(fr_from_u64(indices[j])));

                blst::Fr diff{};
                blst::Fr::sub(diff, i_fr, fr_from_u64(indices[j]));
                denominator = fr_mul(denominator, diff);
            }

            blst::Fr denominator_inverse{};
            blst::Fr::reciprocal(denominator_inverse, denominator);
            lambdas[i] = fr_mul(numerator, denominator_inverse);
        }
        return lambdas;
    }

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

    blst::G2Point scale_g2(const blst::G2Point& point, const blst::Fr& scalar)
    {
        const std::array<std::uint8_t, 32> scalar_bytes = fr_to_scalar_bytes(scalar);
        blst::G2Point out{};
        blst::G2Point::mult(out, point, scalar_bytes.data(), 256);
        return out;
    }

    blst::G1Point sub_g1(const blst::G1Point& lhs, const blst::G1Point& rhs)
    {
        blst::G1Point neg_rhs{};
        blst::G1Point::neg(neg_rhs, rhs);
        blst::G1Point out{};
        blst::G1Point::dadd(out, lhs, neg_rhs);
        return out;
    }

    blst::G1Point hash_to_g1(const Bytes& message, const std::string& dst)
    {
        return bte::building_blocks::hash_to_g1(message, dst);
    }

    Bytes g2_to_bytes(const blst::G2Point& point)
    {
        std::array<std::uint8_t, 96> compressed{};
        blst::G2Point::compress(compressed.data(), point);
        return Bytes(compressed.begin(), compressed.end());
    }

    std::vector<ValidatorSecretShare> make_linear_secret_shares(
        const std::uint32_t count,
        const blst::Fr& secret,
        const blst::Fr& slope)
    {
        std::vector<ValidatorSecretShare> shares;
        shares.reserve(count);
        for (std::uint32_t i = 1; i <= count; ++i) {
            const blst::Fr x = fr_from_u64(i);
            ValidatorSecretShare share{};
            share.index = i;
            share.sk_share = fr_add(secret, fr_mul(slope, x));
            shares.push_back(share);
        }
        return shares;
    }

    Bytes make_associated_data(const std::uint8_t label)
    {
        return Bytes{ std::uint8_t('a'), std::uint8_t('d'), std::uint8_t('-'), label };
    }

    std::string hex_bytes(const std::span<const std::uint8_t> bytes)
    {
        std::ostringstream out;
        out << std::hex << std::setfill('0');
        for (const std::uint8_t byte : bytes) {
            out << std::setw(2) << static_cast<unsigned>(byte);
        }
        return out.str();
    }

    std::string hex_fr(const blst::Fr& value)
    {
        const std::array<std::uint8_t, 32> bytes = fr_to_scalar_bytes(value);
        return hex_bytes(bytes);
    }

    std::string hex_g1(const blst::G1Point& point)
    {
        std::array<std::uint8_t, 48> bytes{};
        blst::G1Point::compress(bytes.data(), point);
        return hex_bytes(bytes);
    }

    std::string hex_g2(const blst::G2Point& point)
    {
        std::array<std::uint8_t, 96> bytes{};
        blst::G2Point::compress(bytes.data(), point);
        return hex_bytes(bytes);
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

    struct TrxFixture
    {
        TrxParameters params;
        InsecureRng rng;
        Crs crs;
        std::vector<ValidatorSecretShare> secret_shares;
        std::vector<ValidatorPublicShare> public_shares;
        PublicKeyMaterial public_key_material;
        EncryptionKey encryption_key;

        TrxFixture()
            : params{8, 2, 32, 1,
                     "BTE_TRX_H1_V1",
                     "BTE_TRX_HF_V1",
                     "BTE_TRX_KDF_SALT_V1",
                     "BTE_TRX_KDF_INFO_V1",
                     "BTE_TRX_OTS_BLS_SIG_V1"},
              rng(0x123456789abcdefULL),
              crs(GenerateCRS(params, [this](std::span<std::uint8_t> out) { fill_random(out); }))
        {
            const blst::Fr secret = fr_from_u64(17);
            const blst::Fr slope = fr_from_u64(29);
            secret_shares = make_linear_secret_shares(4, secret, slope);

            for (const ValidatorSecretShare& share : secret_shares) {
                public_shares.push_back(DeriveValidatorPublicShare(share, crs.h_tau));
            }

            const std::array<std::uint32_t, 2> threshold_subset = { 1, 2 };
            public_key_material = DerivePublicKeyMaterial(
                std::span<const ValidatorPublicShare>(public_shares.data(), public_shares.size()),
                std::span<const std::uint32_t>(threshold_subset.data(), threshold_subset.size()));
            encryption_key = ToEncryptionKey(public_key_material);
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
    };
}

TEST(Trx, EndToEndRoundTrip)
{
    TrxFixture fixture;

    std::vector<EncryptedCiphertext> encrypted;
    for (std::size_t i = 0; i < 3; ++i) {
        encrypted.push_back(Encrypt(
            fixture.params,
            fixture.encryption_key,
            make_associated_data(static_cast<std::uint8_t>('0' + i)),
            [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); }));
    }

    const std::vector<CiphertextEnvelope> ciphertexts = to_envelopes(encrypted);

    const BatchDigest digest = BuildBatchDigest(
        fixture.params,
        std::span<const CiphertextEnvelope>(ciphertexts.data(), ciphertexts.size()),
        0,
        fixture.crs);
    ASSERT_EQ(digest.valid_ciphertexts.size(), ciphertexts.size());

    const std::vector<blst::G1Point> proofs = ComputeEvalProofs(digest, fixture.crs);
    ASSERT_EQ(proofs.size(), ciphertexts.size());
    for (std::size_t i = 0; i < proofs.size(); ++i) {
        EXPECT_TRUE(VerifyEvaluationProof(digest, i, proofs[i], fixture.crs));
    }

    const PartialDecryptionShare share1 = CreatePartialDecryption(
        fixture.params,
        fixture.secret_shares[0],
        digest,
        fixture.public_key_material);
    const PartialDecryptionShare share2 = CreatePartialDecryption(
        fixture.params,
        fixture.secret_shares[1],
        digest,
        fixture.public_key_material);

    EXPECT_TRUE(VerifyPartialDecryption(
        fixture.params,
        share1,
        fixture.public_shares[0],
        digest,
        fixture.public_key_material));
    EXPECT_TRUE(VerifyPartialDecryption(
        fixture.params,
        share2,
        fixture.public_shares[1],
        digest,
        fixture.public_key_material));

    const std::array<PartialDecryptionShare, 2> shares = { share1, share2 };
    const blst::G1Point sigma = CombineShares(
        std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
        digest,
        fixture.public_key_material);
    const blst::G1Point h1 = hash_to_g1(g2_to_bytes(fixture.public_key_material.pk), fixture.params.h1_dst);
    const blst::G1Point base = sub_g1(h1, digest.commitment);
    const blst::G1Point sigma_expected = scale_g1(base, fr_from_u64(17));
    EXPECT_TRUE(blst::G1Point::equal(sigma, sigma_expected));

    for (std::size_t i = 0; i < digest.valid_ciphertexts.size(); ++i) {
        const auto recovered = DecapsulateOne(
            fixture.params,
            digest,
            i,
            proofs[i],
            sigma,
            fixture.crs,
            fixture.public_key_material);
        ASSERT_TRUE(recovered.has_value());
        EXPECT_EQ(*recovered, encrypted[i].session_key);
    }
}

TEST(Trx, ThresholdBoundary)
{
    TrxFixture fixture;
    const EncryptedCiphertext encrypted = Encrypt(
        fixture.params,
        fixture.encryption_key,
        make_associated_data('t'),
        [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); });

    const std::array<CiphertextEnvelope, 1> ciphertexts = { CiphertextEnvelope{ encrypted.ciphertext, encrypted.associated_data } };
    const BatchDigest digest = BuildBatchDigest(
        fixture.params,
        std::span<const CiphertextEnvelope>(ciphertexts.data(), ciphertexts.size()),
        0,
        fixture.crs);

    const PartialDecryptionShare share1 = CreatePartialDecryption(
        fixture.params,
        fixture.secret_shares[0],
        digest,
        fixture.public_key_material);
    const std::array<PartialDecryptionShare, 1> insufficient = { share1 };
    EXPECT_THROW(
        CombineShares(
            std::span<const PartialDecryptionShare>(insufficient.data(), insufficient.size()),
            digest,
            fixture.public_key_material),
        std::invalid_argument);
}

TEST(Trx, ContextIsolation)
{
    TrxFixture fixture;

    std::vector<EncryptedCiphertext> encrypted;
    for (std::size_t i = 0; i < 2; ++i) {
        encrypted.push_back(Encrypt(
            fixture.params,
            fixture.encryption_key,
            make_associated_data(static_cast<std::uint8_t>('A' + i)),
            [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); }));
    }

    const std::vector<CiphertextEnvelope> ciphertexts = to_envelopes(encrypted);

    const BatchDigest digest0 = BuildBatchDigest(
        fixture.params,
        std::span<const CiphertextEnvelope>(ciphertexts.data(), ciphertexts.size()),
        0,
        fixture.crs);
    const BatchDigest digest1 = BuildBatchDigest(
        fixture.params,
        std::span<const CiphertextEnvelope>(ciphertexts.data(), ciphertexts.size()),
        1,
        fixture.crs);

    EXPECT_FALSE(blst::G1Point::equal(digest0.commitment, digest1.commitment));

    const std::vector<blst::G1Point> proofs0 = ComputeEvalProofs(digest0, fixture.crs);
    const std::vector<blst::G1Point> proofs1 = ComputeEvalProofs(digest1, fixture.crs);
    EXPECT_FALSE(blst::G1Point::equal(proofs0.front(), proofs1.front()));

    const std::array<PartialDecryptionShare, 2> shares0 = {
        CreatePartialDecryption(fixture.params, fixture.secret_shares[0], digest0, fixture.public_key_material),
        CreatePartialDecryption(fixture.params, fixture.secret_shares[1], digest0, fixture.public_key_material)
    };
    const std::array<PartialDecryptionShare, 2> shares1 = {
        CreatePartialDecryption(fixture.params, fixture.secret_shares[0], digest1, fixture.public_key_material),
        CreatePartialDecryption(fixture.params, fixture.secret_shares[1], digest1, fixture.public_key_material)
    };

    const blst::G1Point sigma0 = CombineShares(
        std::span<const PartialDecryptionShare>(shares0.data(), shares0.size()),
        digest0,
        fixture.public_key_material);
    const blst::G1Point sigma1 = CombineShares(
        std::span<const PartialDecryptionShare>(shares1.data(), shares1.size()),
        digest1,
        fixture.public_key_material);

    const auto recovered_mixed_proof = DecapsulateOne(
        fixture.params,
        digest0,
        0,
        proofs1.front(),
        sigma0,
        fixture.crs,
        fixture.public_key_material);
    EXPECT_FALSE(recovered_mixed_proof.has_value());

    const auto recovered_mixed_sigma = DecapsulateOne(
        fixture.params,
        digest0,
        0,
        proofs0.front(),
        sigma1,
        fixture.crs,
        fixture.public_key_material);
    ASSERT_TRUE(recovered_mixed_sigma.has_value());
    EXPECT_NE(*recovered_mixed_sigma, encrypted.front().session_key);
}

TEST(Trx, VerifyCiphertextRejectsTampering)
{
    TrxFixture fixture;
    EncryptedCiphertext encrypted = Encrypt(
        fixture.params,
        fixture.encryption_key,
        make_associated_data('x'),
        [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); });

    encrypted.ciphertext.ct3[0] ^= 0x80U;
    EXPECT_FALSE(VerifyCiphertext(
        fixture.params,
        encrypted.ciphertext,
        std::span<const std::uint8_t>(encrypted.associated_data.data(), encrypted.associated_data.size())).has_value());
}

TEST(Trx, VerifyCiphertextRejectsMalformedCt1)
{
    TrxFixture fixture;
    EncryptedCiphertext encrypted = Encrypt(
        fixture.params,
        fixture.encryption_key,
        make_associated_data('y'),
        [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); });

    encrypted.ciphertext.ct1.fill(0xffU);
    EXPECT_FALSE(VerifyCiphertext(
        fixture.params,
        encrypted.ciphertext,
        std::span<const std::uint8_t>(encrypted.associated_data.data(), encrypted.associated_data.size())).has_value());
}

TEST(Trx, BuildBatchDigestRejectsDuplicateTags)
{
    TrxFixture fixture;
    const EncryptedCiphertext encrypted = Encrypt(
        fixture.params,
        fixture.encryption_key,
        make_associated_data('d'),
        [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); });

    const std::array<CiphertextEnvelope, 2> ciphertexts = {
        CiphertextEnvelope{ encrypted.ciphertext, encrypted.associated_data },
        CiphertextEnvelope{ encrypted.ciphertext, encrypted.associated_data }
    };

    EXPECT_THROW(
        BuildBatchDigest(
            fixture.params,
            std::span<const CiphertextEnvelope>(ciphertexts.data(), ciphertexts.size()),
            0,
            fixture.crs),
        std::runtime_error);
}

TEST(Trx, FastLagrangeMatchesNaive)
{
    const std::array<std::uint32_t, 6> indices = { 1, 3, 4, 7, 9, 12 };
    const std::vector<blst::Fr> fast = bte::building_blocks::lagrange_at_zero(indices);
    const std::vector<blst::Fr> naive = lagrange_at_zero_naive(indices);

    ASSERT_EQ(fast.size(), naive.size());
    for (std::size_t i = 0; i < fast.size(); ++i) {
        EXPECT_TRUE(blst::Fr::equal(fast[i], naive[i]));
    }
}

TEST(Trx, FkEvalProofsMatchNaiveCommitments)
{
    TrxFixture fixture;

    std::vector<EncryptedCiphertext> encrypted;
    for (std::size_t i = 0; i < 4; ++i) {
        encrypted.push_back(Encrypt(
            fixture.params,
            fixture.encryption_key,
            make_associated_data(static_cast<std::uint8_t>('f' + i)),
            [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); }));
    }

    const std::vector<CiphertextEnvelope> ciphertexts = to_envelopes(encrypted);
    const BatchDigest digest = BuildBatchDigest(
        fixture.params,
        std::span<const CiphertextEnvelope>(ciphertexts.data(), ciphertexts.size()),
        0,
        fixture.crs);

    const std::vector<blst::G1Point> fk_proofs = ComputeEvalProofs(digest, fixture.crs);
    const std::vector<std::vector<blst::Fr>> quotients =
        bte::building_blocks::exclusion_quotients(
            std::span<const blst::Fr>(digest.tags.data(), digest.tags.size()));

    ASSERT_EQ(fk_proofs.size(), quotients.size());
    for (std::size_t i = 0; i < fk_proofs.size(); ++i) {
        const blst::G1Point naive = bte::building_blocks::commit_polynomial(
            std::span<const blst::G1Point>(
                fixture.crs.contexts[digest.context_id].g1_powers.data(),
                fixture.crs.contexts[digest.context_id].g1_powers.size()),
            quotients[i]);
        EXPECT_TRUE(blst::G1Point::equal(fk_proofs[i], naive));
        EXPECT_TRUE(VerifyEvaluationProof(digest, i, fk_proofs[i], fixture.crs));
    }
}

TEST(Trx, VerifyCiphertextRejectsMalformedVerificationKey)
{
    TrxFixture fixture;
    EncryptedCiphertext encrypted = Encrypt(
        fixture.params,
        fixture.encryption_key,
        make_associated_data('v'),
        [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); });

    encrypted.ciphertext.vk_sig.fill(0xffU);
    EXPECT_FALSE(VerifyCiphertext(
        fixture.params,
        encrypted.ciphertext,
        std::span<const std::uint8_t>(encrypted.associated_data.data(), encrypted.associated_data.size())).has_value());
}

TEST(Trx, BatchSizeSweep)
{
    for (const std::size_t batch_size : { std::size_t(1), std::size_t(3), std::size_t(4), std::size_t(8) }) {
        TrxFixture fixture;

        std::vector<EncryptedCiphertext> encrypted;
        for (std::size_t i = 0; i < batch_size; ++i) {
            encrypted.push_back(Encrypt(
                fixture.params,
                fixture.encryption_key,
                make_associated_data(static_cast<std::uint8_t>('k' + i)),
                [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); }));
        }

        const std::vector<CiphertextEnvelope> ciphertexts = to_envelopes(encrypted);

        const BatchDigest digest = BuildBatchDigest(
            fixture.params,
            std::span<const CiphertextEnvelope>(ciphertexts.data(), ciphertexts.size()),
            0,
            fixture.crs);
        const std::vector<blst::G1Point> proofs = ComputeEvalProofs(digest, fixture.crs);
        const std::array<PartialDecryptionShare, 2> shares = {
            CreatePartialDecryption(fixture.params, fixture.secret_shares[0], digest, fixture.public_key_material),
            CreatePartialDecryption(fixture.params, fixture.secret_shares[1], digest, fixture.public_key_material)
        };
        const blst::G1Point sigma = CombineShares(
            std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
            digest,
            fixture.public_key_material);

        ASSERT_EQ(digest.valid_ciphertexts.size(), batch_size);
        ASSERT_EQ(proofs.size(), batch_size);
        for (std::size_t i = 0; i < batch_size; ++i) {
            EXPECT_TRUE(VerifyEvaluationProof(digest, i, proofs[i], fixture.crs));
            const auto recovered = DecapsulateOne(
                fixture.params,
                digest,
                i,
                proofs[i],
                sigma,
                fixture.crs,
                fixture.public_key_material);
            ASSERT_TRUE(recovered.has_value());
            EXPECT_EQ(*recovered, encrypted[i].session_key);
        }
    }
}

TEST(Trx, VerifyCiphertextRejectsMalformedCt2)
{
    TrxFixture fixture;
    EncryptedCiphertext encrypted = Encrypt(
        fixture.params,
        fixture.encryption_key,
        make_associated_data('2'),
        [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); });

    encrypted.ciphertext.ct2.fill(0xffU);
    EXPECT_FALSE(VerifyCiphertext(
        fixture.params,
        encrypted.ciphertext,
        std::span<const std::uint8_t>(encrypted.associated_data.data(), encrypted.associated_data.size())).has_value());
}

TEST(Trx, VerifyCiphertextRejectsWrongAssociatedData)
{
    TrxFixture fixture;
    EncryptedCiphertext encrypted = Encrypt(
        fixture.params,
        fixture.encryption_key,
        make_associated_data('q'),
        [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); });

    const Bytes wrong_ad = make_associated_data('r');
    EXPECT_FALSE(VerifyCiphertext(
        fixture.params,
        encrypted.ciphertext,
        std::span<const std::uint8_t>(wrong_ad.data(), wrong_ad.size())).has_value());
}

TEST(Trx, VerifyPartialDecryptionRejectsWrongValidatorIndex)
{
    TrxFixture fixture;
    const EncryptedCiphertext encrypted = Encrypt(
        fixture.params,
        fixture.encryption_key,
        make_associated_data('i'),
        [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); });

    const std::array<CiphertextEnvelope, 1> ciphertexts = { CiphertextEnvelope{ encrypted.ciphertext, encrypted.associated_data } };
    const BatchDigest digest = BuildBatchDigest(
        fixture.params,
        std::span<const CiphertextEnvelope>(ciphertexts.data(), ciphertexts.size()),
        0,
        fixture.crs);

    PartialDecryptionShare share = CreatePartialDecryption(
        fixture.params,
        fixture.secret_shares[0],
        digest,
        fixture.public_key_material);
    share.index = fixture.public_shares[1].index;

    EXPECT_FALSE(VerifyPartialDecryption(
        fixture.params,
        share,
        fixture.public_shares[0],
        digest,
        fixture.public_key_material));
}

TEST(Trx, VerifyPartialDecryptionRejectsWrongCommitment)
{
    TrxFixture fixture;
    const EncryptedCiphertext encrypted = Encrypt(
        fixture.params,
        fixture.encryption_key,
        make_associated_data('c'),
        [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); });

    const std::array<CiphertextEnvelope, 1> ciphertexts = { CiphertextEnvelope{ encrypted.ciphertext, encrypted.associated_data } };
    const BatchDigest digest = BuildBatchDigest(
        fixture.params,
        std::span<const CiphertextEnvelope>(ciphertexts.data(), ciphertexts.size()),
        0,
        fixture.crs);

    const PartialDecryptionShare share = CreatePartialDecryption(
        fixture.params,
        fixture.secret_shares[0],
        digest,
        fixture.public_key_material);
    BatchDigest wrong_digest = digest;
    wrong_digest.commitment = fixture.crs.contexts[1].g1_powers.front();

    EXPECT_FALSE(VerifyPartialDecryption(
        fixture.params,
        share,
        fixture.public_shares[0],
        wrong_digest,
        fixture.public_key_material));
}

TEST(Trx, DecapsulateRejectsShuffledProof)
{
    TrxFixture fixture;

    std::vector<EncryptedCiphertext> encrypted;
    for (std::size_t i = 0; i < 2; ++i) {
        encrypted.push_back(Encrypt(
            fixture.params,
            fixture.encryption_key,
            make_associated_data(static_cast<std::uint8_t>('s' + i)),
            [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); }));
    }

    const std::vector<CiphertextEnvelope> ciphertexts = to_envelopes(encrypted);
    const BatchDigest digest = BuildBatchDigest(
        fixture.params,
        std::span<const CiphertextEnvelope>(ciphertexts.data(), ciphertexts.size()),
        0,
        fixture.crs);
    const std::vector<blst::G1Point> proofs = ComputeEvalProofs(digest, fixture.crs);
    const std::array<PartialDecryptionShare, 2> shares = {
        CreatePartialDecryption(fixture.params, fixture.secret_shares[0], digest, fixture.public_key_material),
        CreatePartialDecryption(fixture.params, fixture.secret_shares[1], digest, fixture.public_key_material)
    };
    const blst::G1Point sigma = CombineShares(
        std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
        digest,
        fixture.public_key_material);

    const auto recovered = DecapsulateOne(
        fixture.params,
        digest,
        0,
        proofs[1],
        sigma,
        fixture.crs,
        fixture.public_key_material);
    EXPECT_FALSE(recovered.has_value());
}

TEST(Trx, WrongPkTauFailsRecovery)
{
    TrxFixture fixture;
    PublicKeyMaterial wrong_material = fixture.public_key_material;
    wrong_material.pk_tau = scale_g2(blst::G2Point::generator(), fr_from_u64(777));
    const EncryptionKey wrong_encryption_key = ToEncryptionKey(wrong_material);

    const EncryptedCiphertext encrypted = Encrypt(
        fixture.params,
        wrong_encryption_key,
        make_associated_data('p'),
        [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); });

    const std::array<CiphertextEnvelope, 1> ciphertexts = { CiphertextEnvelope{ encrypted.ciphertext, encrypted.associated_data } };
    const BatchDigest digest = BuildBatchDigest(
        fixture.params,
        std::span<const CiphertextEnvelope>(ciphertexts.data(), ciphertexts.size()),
        0,
        fixture.crs);
    const std::vector<blst::G1Point> proofs = ComputeEvalProofs(digest, fixture.crs);
    const std::array<PartialDecryptionShare, 2> shares = {
        CreatePartialDecryption(fixture.params, fixture.secret_shares[0], digest, fixture.public_key_material),
        CreatePartialDecryption(fixture.params, fixture.secret_shares[1], digest, fixture.public_key_material)
    };
    const blst::G1Point sigma = CombineShares(
        std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
        digest,
        fixture.public_key_material);

    const auto recovered = DecapsulateOne(
        fixture.params,
        digest,
        0,
        proofs[0],
        sigma,
        fixture.crs,
        fixture.public_key_material);
    ASSERT_TRUE(recovered.has_value());
    EXPECT_NE(*recovered, encrypted.session_key);
}

TEST(Trx, DeterministicVector)
{
    TrxFixture fixture;

    const EncryptedCiphertext encrypted = Encrypt(
        fixture.params,
        fixture.encryption_key,
        make_associated_data('0'),
        [&fixture](std::span<std::uint8_t> out) { fixture.fill_random(out); });

    const std::array<CiphertextEnvelope, 1> ciphertexts = {
        CiphertextEnvelope{ encrypted.ciphertext, encrypted.associated_data }
    };
    const BatchDigest digest = BuildBatchDigest(
        fixture.params,
        std::span<const CiphertextEnvelope>(ciphertexts.data(), ciphertexts.size()),
        0,
        fixture.crs);
    const std::vector<blst::G1Point> proofs = ComputeEvalProofs(digest, fixture.crs);
    const std::array<PartialDecryptionShare, 2> shares = {
        CreatePartialDecryption(fixture.params, fixture.secret_shares[0], digest, fixture.public_key_material),
        CreatePartialDecryption(fixture.params, fixture.secret_shares[1], digest, fixture.public_key_material)
    };
    const blst::G1Point sigma = CombineShares(
        std::span<const PartialDecryptionShare>(shares.data(), shares.size()),
        digest,
        fixture.public_key_material);
    const auto recovered = DecapsulateOne(
        fixture.params,
        digest,
        0,
        proofs[0],
        sigma,
        fixture.crs,
        fixture.public_key_material);
    ASSERT_TRUE(recovered.has_value());

    EXPECT_EQ(hex_g2(fixture.public_key_material.pk), "ad05ceb0be53d2624a796a7a033aec59d9463c18d672c451ec4f2e679daef882cab7d8dd88789065156a1340ca9d42650ef786ebdcda12e142a32f091307f2fedf52f6c36beb278b0007a03ad81bf9fee3710a04928e43e541d02c9be44722e8");
    EXPECT_EQ(hex_g2(fixture.public_key_material.pk_tau), "b98c4d04ca62a7a981706896735289f22a615ab9c57be627f320104b15399558d06a43c466a9f59f2c92a2ea08c269e2007d4236a298b10f50a33cdb5b70372186157eaaf817b249faea146ba3e283767e8c71a5ced75e709a5ac0a15e521f59");
    EXPECT_EQ(hex_fr(digest.tags[0]), "85f2281f5fc7ba96f2c1fe9f8d91da1a68eb1a84501dc01432e296c062c82b34");
    EXPECT_EQ(hex_fr(digest.polynomial[0]), "7c0dd7e09f3845690c9aff5f7512e3389dec8685b7ba791e169b0669f0dec13f");
    EXPECT_EQ(hex_fr(digest.polynomial[1]), "0100000000000000000000000000000000000000000000000000000000000000");
    EXPECT_EQ(hex_g1(digest.commitment), "8dc08fac8fe2132a73d92f523bd10f8b090d3400129cae66912b5eb0209d8f8133424b504c8c7dbcbbb6700f3834b6c2");
    EXPECT_EQ(hex_g1(proofs[0]), "aa0c82e3b27dcb7b7c173e8bee9d2c962e200f80036486c472f3c7b2f5258f0e6e9455eab71596e31b9d17b5d530b745");
    EXPECT_EQ(hex_g1(shares[0].value), "8ea3e2dd121e6f266ebf30c4ec55042d760fa23abc41e93917d19f54b84864e6840d026db302b6b177057644ae033760");
    EXPECT_EQ(hex_g1(sigma), "8eac71106293a2481447956b94ea0099802aa3606b4951bd895a89e518cdb82557e2dd2e320d4e83cb146bd1f08cd6f5");
    EXPECT_EQ(hex_bytes(encrypted.session_key), "9a6dc6905c2d722cd6bb46a706de1a7369ad74fa05e258620a05b34cd7ead539");
    EXPECT_EQ(hex_bytes(*recovered), "9a6dc6905c2d722cd6bb46a706de1a7369ad74fa05e258620a05b34cd7ead539");
}

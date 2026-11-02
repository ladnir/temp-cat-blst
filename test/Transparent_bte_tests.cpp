#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "BTE/building_blocks/blst_cpp.h"
#include "BTE/transparent_bte/core.h"

namespace
{
    blst::Fr fr(const std::uint64_t value)
    {
        return blst::Fr::from_uint64(value);
    }

    blst::GT pairing_generator()
    {
        blst::G1AffinePoint p{};
        blst::G2AffinePoint q{};
        blst::G1Point::to_affine(p, blst::G1Point::generator());
        blst::G2Point::to_affine(q, blst::G2Point::generator());
        blst::MillerResult miller{};
        blst::GT::miller_loop(miller, &q, &p);
        blst::GT out{};
        blst::GT::final_exp(out, miller);
        return out;
    }

    blst::GT gt_multiple(const blst::GT& point, const std::uint64_t scalar)
    {
        const blst::Fr scalar_fr = fr(scalar);
        std::array<std::uint8_t, 32> bytes{};
        blst::Fr::to_scaler(bytes.data(), scalar_fr);
        blst::GT out{};
        blst::GT::mult(out, point, bytes.data(), 255);
        return out;
    }

    bte::transparent_bte::SetupResult setup(const std::size_t servers = 6)
    {
        const std::vector<blst::Fr> sharing = { fr(19), fr(7), fr(3) };
        const std::array<std::uint8_t, 24> label = {
            '2','0','2','6','/','2','2','0','5','/','g','-','p','r','i','m','e','/','t','e','s','t','/','1'
        };
        const blst::G1Point g_prime = bte::building_blocks::hash_to_g1(
            label, "BTE_TRANSPARENT_2026_2205_G_PRIME_V1");
        return bte::transparent_bte::SetupFromSharingPolynomial(
            sharing, servers, g_prime);
    }

    struct BatchFixture
    {
        bte::transparent_bte::SetupResult keys;
        std::vector<blst::GT> messages;
        std::vector<bte::transparent_bte::Ciphertext> ciphertexts;
    };

    BatchFixture make_batch(const std::size_t count)
    {
        BatchFixture out{ setup(), {}, {} };
        const blst::GT base = pairing_generator();
        out.messages.reserve(count);
        out.ciphertexts.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            out.messages.push_back(gt_multiple(base, 1000 + i));
            out.ciphertexts.push_back(bte::transparent_bte::Encrypt(
                out.keys.public_parameters,
                fr(2000 + 17 * i),
                out.messages.back()));
        }
        return out;
    }

    void expect_gt_vectors_equal(
        const std::span<const blst::GT> lhs,
        const std::span<const blst::GT> rhs)
    {
        ASSERT_EQ(lhs.size(), rhs.size());
        for (std::size_t i = 0; i < lhs.size(); ++i) {
            EXPECT_NE(blst::GT::equal(lhs[i], rhs[i]), 0) << "slot " << i;
        }
    }
}

TEST(TransparentBte, CauchyTransformImplementationsMatchDirectOracle)
{
    for (const std::size_t count : { 1U, 2U, 3U, 8U, 17U, 33U }) {
        const BatchFixture batch = make_batch(count);
        const auto direct = bte::transparent_bte::PrecomputeBatchDecryptionNaive(
            batch.ciphertexts);
        const auto msm = bte::transparent_bte::PrecomputeBatchDecryptionMsm(
            batch.ciphertexts);
        const auto quadratic_pairing =
            bte::transparent_bte::PrecomputeBatchDecryptionQuadraticPairing(
                batch.ciphertexts);
        const auto fast = bte::transparent_bte::PrecomputeBatchDecryptionFast(
            batch.ciphertexts);
        expect_gt_vectors_equal(direct.cross_terms, msm.cross_terms);
        expect_gt_vectors_equal(direct.cross_terms, quadratic_pairing.cross_terms);
        expect_gt_vectors_equal(direct.cross_terms, fast.cross_terms);
    }
}

TEST(TransparentBte, ThresholdRoundTripAndShareVerification)
{
    for (const std::size_t count : { 1U, 4U, 19U }) {
        const BatchFixture batch = make_batch(count);
        std::vector<bte::transparent_bte::PartialDecryptionShare> shares;
        for (std::size_t i = 0; i < batch.keys.public_parameters.threshold; ++i) {
            shares.push_back(bte::transparent_bte::PartialDecrypt(
                batch.keys.secret_shares[i], batch.ciphertexts));
            EXPECT_TRUE(bte::transparent_bte::VerifyPartialDecryptionShare(
                batch.keys.public_parameters, shares.back(), batch.ciphertexts));
        }

        auto malformed = shares.front();
        blst::G1Point::dadd(
            malformed.value, malformed.value, blst::G1Point::generator());
        EXPECT_FALSE(bte::transparent_bte::VerifyPartialDecryptionShare(
            batch.keys.public_parameters, malformed, batch.ciphertexts));

        const blst::G1Point combined =
            bte::transparent_bte::CombinePartialDecryptionShares(
                batch.keys.public_parameters, shares);
        const auto precomputation =
            bte::transparent_bte::PrecomputeBatchDecryptionFast(batch.ciphertexts);
        const auto recovered = bte::transparent_bte::FinalizeBatchDecryption(
            precomputation, combined, batch.ciphertexts);
        expect_gt_vectors_equal(batch.messages, recovered);
    }
}

TEST(TransparentBte, RejectsDuplicateCiphertextPoints)
{
    BatchFixture batch = make_batch(3);
    batch.ciphertexts[2].c1 = batch.ciphertexts[0].c1;
    EXPECT_THROW(
        bte::transparent_bte::PrecomputeBatchDecryptionFast(batch.ciphertexts),
        std::invalid_argument);
}

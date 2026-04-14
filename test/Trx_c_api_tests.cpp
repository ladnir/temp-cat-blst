#include <array>
#include <vector>

#include <gtest/gtest.h>

#include "bindings/blst_trx.h"
#include "blst/fields/Fr.h"
#include "InsecureRng.h"

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

    blst_scalar fr_to_scalar(const blst::Fr& value)
    {
        blst_scalar out{};
        blst::Fr::to_scaler(out.b, value);
        return out;
    }

    std::vector<blst_trx_validator_secret_share> make_linear_secret_shares(
        const std::uint32_t count,
        const blst::Fr& secret,
        const blst::Fr& slope)
    {
        std::vector<blst_trx_validator_secret_share> shares;
        shares.reserve(count);
        for (std::uint32_t i = 1; i <= count; ++i) {
            const blst::Fr x = fr_from_u64(i);
            shares.push_back({
                i,
                fr_to_scalar(fr_add(secret, fr_mul(slope, x)))
            });
        }
        return shares;
    }

    struct TrxCApiFixture
    {
        InsecureRng rng{0x2468ace013579bdfULL};

        static bool fill_random(void* user_data, byte* out, const size_t out_len)
        {
            auto* self = static_cast<TrxCApiFixture*>(user_data);
            std::uint64_t word = 0;
            std::size_t remaining = 0;
            for (std::size_t i = 0; i < out_len; ++i) {
                if (remaining == 0) {
                    word = self->rng.next_u64();
                    remaining = sizeof(word);
                }
                out[i] = static_cast<byte>(word & 0xffU);
                word >>= 8;
                --remaining;
            }
            return true;
        }
    };
}

TEST(TrxCApi, EndToEndRoundTrip)
{
    TrxCApiFixture fixture;
    const blst_trx_params params{
        8,
        2,
        32,
        1,
        "BTE_TRX_H1_V1",
        "BTE_TRX_HF_V1",
        "BTE_TRX_KDF_SALT_V1",
        "BTE_TRX_KDF_INFO_V1",
        "BTE_TRX_OTS_BLS_SIG_V1",
    };

    blst_trx_crs* const crs = blst_trx_generate_crs(&params, &TrxCApiFixture::fill_random, &fixture);
    ASSERT_NE(crs, nullptr);

    const std::vector<blst_trx_validator_secret_share> secret_shares =
        make_linear_secret_shares(4, fr_from_u64(17), fr_from_u64(29));
    std::vector<blst_trx_validator_public_share> public_shares(secret_shares.size());
    for (std::size_t i = 0; i < secret_shares.size(); ++i) {
        ASSERT_TRUE(blst_trx_derive_validator_public_share(&secret_shares[i], crs, &public_shares[i]));
    }

    const std::array<std::uint32_t, 2> threshold_subset = { 1, 2 };
    blst_trx_public_key_material* const public_key_material =
        blst_trx_derive_public_key_material(
            public_shares.data(),
            public_shares.size(),
            threshold_subset.data(),
            threshold_subset.size());
    ASSERT_NE(public_key_material, nullptr);

    blst_trx_encryption_key* const encryption_key =
        blst_trx_to_encryption_key(public_key_material);
    ASSERT_NE(encryption_key, nullptr);

    const std::array<std::array<byte, 4>, 3> associated_data = {{
        { byte('a'), byte('d'), byte('-'), byte('0') },
        { byte('a'), byte('d'), byte('-'), byte('1') },
        { byte('a'), byte('d'), byte('-'), byte('2') },
    }};

    std::array<blst_trx_encrypted_ciphertext*, 3> encrypted = { nullptr, nullptr, nullptr };
    std::array<blst_trx_verified_ciphertext*, 3> verified = { nullptr, nullptr, nullptr };
    std::array<std::vector<byte>, 3> expected_session_keys;

    for (std::size_t i = 0; i < encrypted.size(); ++i) {
        encrypted[i] = blst_trx_encrypt(
            &params,
            encryption_key,
            associated_data[i].data(),
            associated_data[i].size(),
            &TrxCApiFixture::fill_random,
            &fixture);
        ASSERT_NE(encrypted[i], nullptr);

        expected_session_keys[i].resize(blst_trx_encrypted_ciphertext_session_key_size(encrypted[i]));
        ASSERT_TRUE(blst_trx_encrypted_ciphertext_copy_session_key(
            encrypted[i],
            expected_session_keys[i].data(),
            expected_session_keys[i].size()));

        verified[i] = blst_trx_verify_ciphertext(&params, encrypted[i]);
        ASSERT_NE(verified[i], nullptr);
    }

    std::array<const blst_trx_verified_ciphertext*, 3> verified_ptrs = {
        verified[0], verified[1], verified[2]
    };
    blst_trx_batch_digest* const digest =
        blst_trx_build_batch_digest(&params, verified_ptrs.data(), verified_ptrs.size(), 0, crs);
    ASSERT_NE(digest, nullptr);
    ASSERT_EQ(blst_trx_batch_digest_valid_ciphertext_count(digest), verified.size());

    std::vector<blst_p1> proofs(verified.size());
    ASSERT_TRUE(blst_trx_compute_eval_proofs(digest, crs, proofs.data(), proofs.size()));
    for (std::size_t i = 0; i < proofs.size(); ++i) {
        EXPECT_TRUE(blst_trx_verify_evaluation_proof(digest, i, &proofs[i], crs));
    }

    std::array<blst_trx_partial_decryption_share, 2> partial_shares = {};
    for (std::size_t i = 0; i < partial_shares.size(); ++i) {
        ASSERT_TRUE(blst_trx_create_partial_decryption(
            &params,
            &secret_shares[i],
            digest,
            public_key_material,
            &partial_shares[i]));
        EXPECT_TRUE(blst_trx_verify_partial_decryption(
            &params,
            &partial_shares[i],
            &public_shares[i],
            digest,
            public_key_material));
    }

    blst_p1 sigma = {};
    ASSERT_TRUE(blst_trx_combine_shares(
        partial_shares.data(),
        partial_shares.size(),
        digest,
        public_key_material,
        &sigma));

    for (std::size_t i = 0; i < proofs.size(); ++i) {
        std::vector<byte> recovered(expected_session_keys[i].size());
        ASSERT_TRUE(blst_trx_decapsulate_one(
            &params,
            digest,
            i,
            &proofs[i],
            &sigma,
            crs,
            public_key_material,
            recovered.data(),
            recovered.size()));
        EXPECT_EQ(recovered, expected_session_keys[i]);
    }

    blst_trx_batch_digest_free(digest);
    for (blst_trx_verified_ciphertext*& item : verified) {
        blst_trx_verified_ciphertext_free(item);
    }
    for (blst_trx_encrypted_ciphertext*& item : encrypted) {
        blst_trx_encrypted_ciphertext_free(item);
    }
    blst_trx_encryption_key_free(encryption_key);
    blst_trx_public_key_material_free(public_key_material);
    blst_trx_crs_free(crs);
}

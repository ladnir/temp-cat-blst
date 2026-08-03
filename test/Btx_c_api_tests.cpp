#include <vector>

#include <gtest/gtest.h>

#include "bindings/blst_btx.h"
#include "InsecureRng.h"

namespace
{
    struct BtxCApiFixture
    {
        InsecureRng rng{0x13579bdf2468ace0ULL};

        static bool fill_random(void* user_data, blst_byte* out, const size_t out_len)
        {
            auto* self = static_cast<BtxCApiFixture*>(user_data);
            std::uint64_t word = 0;
            std::size_t remaining = 0;
            for (std::size_t i = 0; i < out_len; ++i) {
                if (remaining == 0) {
                    word = self->rng.next_u64();
                    remaining = sizeof(word);
                }
                out[i] = static_cast<blst_byte>(word & 0xffU);
                word >>= 8;
                --remaining;
            }
            return true;
        }
    };
}

TEST(BtxCApi, PreparedFlowRecoversEncryptedKeys)
{
    BtxCApiFixture fixture;
    const blst_btx_params params{
        8,
        4,
        2,
        32,
        "BTE_BTX_KDF_SALT_V1",
        "BTE_BTX_KDF_INFO_V1",
    };

    blst_btx_setup_result* const setup =
        blst_btx_setup(&params, &BtxCApiFixture::fill_random, &fixture);
    ASSERT_NE(setup, nullptr);
    const blst_btx_public_parameters* const public_parameters =
        blst_btx_setup_result_public_parameters(setup);
    ASSERT_NE(public_parameters, nullptr);
    ASSERT_EQ(blst_btx_setup_result_secret_share_count(setup), 4U);

    blst_btx_prepared_public_parameters* const prepared =
        blst_btx_prepare_public_parameters(public_parameters);
    ASSERT_NE(prepared, nullptr);

    std::vector<blst_p1> ciphertexts(4);
    std::vector<blst_byte> expected_keys(4 * params.key_size);
    for (std::size_t i = 0; i < ciphertexts.size(); ++i) {
        ASSERT_TRUE(blst_btx_encrypt(
            &params,
            public_parameters,
            &BtxCApiFixture::fill_random,
            &fixture,
            &ciphertexts[i],
            expected_keys.data() + (i * params.key_size),
            params.key_size));
    }

    std::vector<blst_btx_partial_decryption_share> shares(2);
    for (std::size_t i = 0; i < shares.size(); ++i) {
        const blst_btx_server_secret_share* const secret_share =
            blst_btx_setup_result_secret_share(setup, i);
        ASSERT_NE(secret_share, nullptr);
        ASSERT_TRUE(blst_btx_create_partial_decryption_share(
            secret_share,
            ciphertexts.data(),
            ciphertexts.size(),
            &shares[i]));
    }

    blst_p1 aggregate_predecryption_key = {};
    ASSERT_TRUE(blst_btx_combine_partial_decryption_shares(
        shares.data(),
        shares.size(),
        public_parameters,
        &aggregate_predecryption_key));

    blst_btx_batch_precomputation* const precomputation =
        blst_btx_precompute_batch_decryption(
            prepared,
            ciphertexts.data(),
            ciphertexts.size(),
            BLST_BTX_CROSS_TERM_FFT);
    ASSERT_NE(precomputation, nullptr);
    ASSERT_EQ(blst_btx_batch_precomputation_batch_size(precomputation), ciphertexts.size());

    std::vector<blst_byte> finalized_keys(ciphertexts.size() * params.key_size);
    ASSERT_TRUE(blst_btx_finalize_batch_decryption(
        &params,
        prepared,
        precomputation,
        &aggregate_predecryption_key,
        finalized_keys.data(),
        finalized_keys.size()));
    EXPECT_EQ(finalized_keys, expected_keys);

    std::vector<blst_byte> recovered_keys(ciphertexts.size() * params.key_size);
    ASSERT_TRUE(blst_btx_recover_batch_keys(
        &params,
        prepared,
        ciphertexts.data(),
        ciphertexts.size(),
        &aggregate_predecryption_key,
        BLST_BTX_CROSS_TERM_FFT,
        recovered_keys.data(),
        recovered_keys.size()));
    EXPECT_EQ(recovered_keys, expected_keys);

    blst_btx_batch_precomputation_free(precomputation);
    blst_btx_prepared_public_parameters_free(prepared);
    blst_btx_setup_result_free(setup);
}

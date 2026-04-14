#include <gtest/gtest.h>

#include <array>
#include <cstring>

#include "bindings/blst.h"
#include "bindings/blst_btx.h"
#include "bindings/blst_trx.h"

TEST(BindingHeader, BlstUmbrellaExposesKeygenAndPkApis)
{
    std::array<byte, 32> ikm = {};
    for (size_t i = 0; i < ikm.size(); ++i) {
        ikm[i] = static_cast<byte>(i + 1);
    }

    blst_scalar sk = {};
    blst_p1 pk_g1 = {};
    blst_p2 pk_g2 = {};

    blst_keygen(&sk, ikm.data(), ikm.size(), nullptr, 0);
    blst_sk_to_pk_in_g1(&pk_g1, &sk);
    blst_sk_to_pk_in_g2(&pk_g2, &sk);

    EXPECT_FALSE(blst_p1_is_inf(&pk_g1));
    EXPECT_TRUE(blst_p1_in_g1(&pk_g1));
    EXPECT_FALSE(blst_p2_is_inf(&pk_g2));
    EXPECT_TRUE(blst_p2_in_g2(&pk_g2));
}

TEST(BindingHeader, KeygenIsDeterministicForSameInputs)
{
    std::array<byte, 64> ikm = {};
    for (size_t i = 0; i < ikm.size(); ++i) {
        ikm[i] = static_cast<byte>(0x7f ^ i);
    }

    blst_scalar sk1 = {};
    blst_scalar sk2 = {};
    blst_p1 pk1 = {};
    blst_p1 pk2 = {};

    blst_keygen(&sk1, ikm.data(), ikm.size(), nullptr, 0);
    blst_keygen(&sk2, ikm.data(), ikm.size(), nullptr, 0);

    EXPECT_EQ(std::memcmp(&sk1, &sk2, sizeof(sk1)), 0);

    blst_sk_to_pk_in_g1(&pk1, &sk1);
    blst_sk_to_pk_in_g1(&pk2, &sk2);
    EXPECT_TRUE(blst_p1_is_equal(&pk1, &pk2));
}

TEST(BindingHeader, DeriveMasterIsDeterministicForSameInputs)
{
    std::array<byte, 64> ikm = {};
    for (size_t i = 0; i < ikm.size(); ++i) {
        ikm[i] = static_cast<byte>(0xa5 + i);
    }

    blst_scalar sk1 = {};
    blst_scalar sk2 = {};

    blst_derive_master_eip2333(&sk1, ikm.data(), ikm.size());
    blst_derive_master_eip2333(&sk2, ikm.data(), ikm.size());

    EXPECT_EQ(std::memcmp(&sk1, &sk2, sizeof(sk1)), 0);
}

TEST(BindingHeader, ScalarHelperSurfaceIsExposedAndStable)
{
    EXPECT_EQ(static_cast<int>(BLST_BAD_SCALAR), 7);

    std::array<byte, 1> one_be = {1};
    blst_scalar one = {};
    blst_scalar two = {};
    blst_scalar zero = {};

    EXPECT_TRUE(blst_scalar_from_be_bytes(&one, one_be.data(), one_be.size()));
    EXPECT_TRUE(blst_sk_add_n_check(&two, &one, &one));
    EXPECT_TRUE(blst_sk_sub_n_check(&one, &two, &one));
    EXPECT_TRUE(blst_sk_mul_n_check(&two, &one, &two));
    EXPECT_FALSE(blst_scalar_from_be_bytes(&zero, nullptr, 0));
}

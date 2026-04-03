#include <gtest/gtest.h>

#include <array>

#include "bindings/blst.h"

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

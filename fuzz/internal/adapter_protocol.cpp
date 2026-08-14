// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter_common.h"

#include <vector>

namespace catblst::fuzz::adapter
{
    void run_protocol_case(std::span<const std::byte> input)
    {
        using detail::object_equal;
        using detail::require;
        input = input.first(std::min<std::size_t>(input.size(), 2048));
        const auto message = input.first(std::min<std::size_t>(input.size(), 1024));
        const auto dst = input.first(std::min<std::size_t>(input.size(), 512));
        const auto aug = input.last(std::min<std::size_t>(input.size(), 512));
        const auto info = input.first(std::min<std::size_t>(input.size(), 256));
        const auto salt = input.last(std::min<std::size_t>(input.size(), 256));

        std::array<blst_byte, 32> digest_a{};
        std::array<blst_byte, 32> digest_b{};
        ::blst_sha256(digest_a.data(), detail::bytes(message), message.size());
        ::blst_sha256(digest_b.data(), detail::bytes(message), message.size());
        require(digest_a == digest_b);

        constexpr std::array<std::size_t, 5> xmd_lengths{0, 1, 32, 8160, 8161};
        const std::size_t selector = input.empty()
            ? 0
            : std::to_integer<unsigned char>(input.front()) % xmd_lengths.size();
        const std::size_t xmd_length = xmd_lengths[selector];
        std::vector<blst_byte> xmd_a(xmd_length, 0xa5);
        std::vector<blst_byte> xmd_b(xmd_length, 0x5a);
        ::blst_expand_message_xmd(
            xmd_a.empty() ? nullptr : xmd_a.data(),
            xmd_a.size(),
            detail::bytes(message),
            message.size(),
            detail::bytes(dst),
            dst.size());
        ::blst_expand_message_xmd(
            xmd_b.empty() ? nullptr : xmd_b.data(),
            xmd_b.size(),
            detail::bytes(message),
            message.size(),
            detail::bytes(dst),
            dst.size());
        require(xmd_a == xmd_b);
        if (xmd_length > 8160)
            require(std::all_of(xmd_a.begin(), xmd_a.end(), [](blst_byte b) {
                return b == 0;
            }));

        std::array<blst_byte, 32> ikm{};
        for (std::size_t index = 0; index < ikm.size(); ++index)
            ikm[index] = input.empty()
                ? static_cast<blst_byte>(index + 1)
                : static_cast<blst_byte>(
                      std::to_integer<unsigned char>(input[index % input.size()]) ^
                      static_cast<unsigned char>(index + 1));

        blst_scalar secret{};
        ::blst_keygen(
            &secret,
            ikm.data(),
            ikm.size(),
            detail::bytes(info),
            info.size());
        require(::blst_sk_check(&secret));
        blst_scalar repeat{};
        ::blst_keygen(
            &repeat,
            ikm.data(),
            ikm.size(),
            detail::bytes(info),
            info.size());
        require(object_equal(secret, repeat));

        blst_scalar v3{};
        blst_scalar v45{};
        blst_scalar v5{};
        ::blst_keygen_v3(&v3, ikm.data(), ikm.size(), detail::bytes(info), info.size());
        ::blst_keygen_v4_5(
            &v45,
            ikm.data(),
            ikm.size(),
            detail::bytes(salt),
            salt.size(),
            detail::bytes(info),
            info.size());
        ::blst_keygen_v5(
            &v5,
            ikm.data(),
            ikm.size(),
            salt.empty() ? ikm.data() : detail::bytes(salt),
            salt.size(),
            detail::bytes(info),
            info.size());
        require(::blst_scalar_fr_check(&v3));
        require(::blst_sk_check(&v45));
        require(::blst_sk_check(&v5));

        blst_scalar rejected;
        std::memset(&rejected, 0xa5, sizeof(rejected));
        ::blst_keygen(&rejected, ikm.data(), 31, nullptr, 0);
        const blst_scalar zero_scalar{};
        require(object_equal(rejected, zero_scalar));
        std::memset(&rejected, 0xa5, sizeof(rejected));
        ::blst_keygen_v5(&rejected, ikm.data(), ikm.size(), nullptr, 0, nullptr, 0);
        require(object_equal(rejected, zero_scalar));

        blst_scalar master{};
        blst_scalar default_empty_info{};
        ::blst_derive_master_eip2333(&master, ikm.data(), ikm.size());
        ::blst_keygen(
            &default_empty_info, ikm.data(), ikm.size(), nullptr, 0);
        require(object_equal(master, default_empty_info));
        blst_scalar child{};
        const std::uint32_t child_index = input.size() >= 4
            ? (std::to_integer<std::uint32_t>(input[0]) |
               (std::to_integer<std::uint32_t>(input[1]) << 8) |
               (std::to_integer<std::uint32_t>(input[2]) << 16) |
               (std::to_integer<std::uint32_t>(input[3]) << 24))
            : 0;
        ::blst_derive_child_eip2333(&child, &master, child_index);
        require(::blst_sk_check(&child));

        blst_p1 encoded_g1{};
        blst_p1 hashed_g1{};
        blst_p2 encoded_g2{};
        blst_p2 hashed_g2{};
        ::blst_encode_to_g1(
            &encoded_g1,
            detail::bytes(message), message.size(),
            detail::bytes(dst), dst.size(),
            detail::bytes(aug), aug.size());
        ::blst_hash_to_g1(
            &hashed_g1,
            detail::bytes(message), message.size(),
            detail::bytes(dst), dst.size(),
            detail::bytes(aug), aug.size());
        ::blst_encode_to_g2(
            &encoded_g2,
            detail::bytes(message), message.size(),
            detail::bytes(dst), dst.size(),
            detail::bytes(aug), aug.size());
        ::blst_hash_to_g2(
            &hashed_g2,
            detail::bytes(message), message.size(),
            detail::bytes(dst), dst.size(),
            detail::bytes(aug), aug.size());
        require(::blst_p1_in_g1(&encoded_g1));
        require(::blst_p1_in_g1(&hashed_g1));
        require(::blst_p2_in_g2(&encoded_g2));
        require(::blst_p2_in_g2(&hashed_g2));

        const blst_fp u = detail::native<blst_fp>(blst_fp_from_le_bytes(message));
        const blst_fp v = detail::native<blst_fp>(blst_fp_from_le_bytes(aug));
        const blst_fp2 u2 = detail::native<blst_fp2>(reduce_fp2(message));
        const blst_fp2 v2 = detail::native<blst_fp2>(reduce_fp2(aug));
        blst_p1 mapped_g1{};
        blst_p2 mapped_g2{};
        ::blst_map_to_g1(&mapped_g1, &u, &v);
        ::blst_map_to_g2(&mapped_g2, &u2, &v2);
        require(::blst_p1_in_g1(&mapped_g1));
        require(::blst_p2_in_g2(&mapped_g2));

        blst_p1 pk_g1{};
        blst_p2 pk_g2{};
        blst_p2 signature_g2{};
        blst_p1 signature_g1{};
        ::blst_sk_to_pk_in_g1(&pk_g1, &secret);
        ::blst_sk_to_pk_in_g2(&pk_g2, &secret);
        ::blst_sign_pk_in_g1(&signature_g2, &hashed_g2, &secret);
        ::blst_sign_pk_in_g2(&signature_g1, &hashed_g1, &secret);
        require(::blst_p1_in_g1(&pk_g1));
        require(::blst_p2_in_g2(&pk_g2));
        require(::blst_p2_in_g2(&signature_g2));
        require(::blst_p1_in_g1(&signature_g1));

        blst_p1_affine pk_g1_affine{};
        blst_p2_affine pk_g2_affine{};
        blst_p2_affine signature_g2_affine{};
        blst_p1_affine signature_g1_affine{};
        ::blst_p1_to_affine(&pk_g1_affine, &pk_g1);
        ::blst_p2_to_affine(&pk_g2_affine, &pk_g2);
        ::blst_p2_to_affine(&signature_g2_affine, &signature_g2);
        ::blst_p1_to_affine(&signature_g1_affine, &signature_g1);

        require(::blst_core_verify_pk_in_g1(
                    &pk_g1_affine,
                    &signature_g2_affine,
                    true,
                    detail::bytes(message), message.size(),
                    detail::bytes(dst), dst.size(),
                    detail::bytes(aug), aug.size()) == BLST_SUCCESS);
        require(::blst_core_verify_pk_in_g2(
                    &pk_g2_affine,
                    &signature_g1_affine,
                    true,
                    detail::bytes(message), message.size(),
                    detail::bytes(dst), dst.size(),
                    detail::bytes(aug), aug.size()) == BLST_SUCCESS);

        std::array<blst_byte, 96> pk1_decorated{};
        std::array<blst_byte, 192> pk2_decorated{};
        std::array<blst_byte, 192> sig2_decorated{};
        std::array<blst_byte, 96> sig1_decorated{};
        blst_p1_affine pk1_direct{};
        blst_p2_affine pk2_direct{};
        blst_p2_affine sig2_direct{};
        blst_p1_affine sig1_direct{};
        ::blst_sk_to_pk2_in_g1(pk1_decorated.data(), &pk1_direct, &secret);
        ::blst_sk_to_pk2_in_g2(pk2_decorated.data(), &pk2_direct, &secret);
        ::blst_sign_pk2_in_g1(
            sig2_decorated.data(), &sig2_direct, &hashed_g2, &secret);
        ::blst_sign_pk2_in_g2(
            sig1_decorated.data(), &sig1_direct, &hashed_g1, &secret);
        require(::blst_p1_affine_is_equal(&pk1_direct, &pk_g1_affine));
        require(::blst_p2_affine_is_equal(&pk2_direct, &pk_g2_affine));
        require(::blst_p2_affine_is_equal(&sig2_direct, &signature_g2_affine));
        require(::blst_p1_affine_is_equal(&sig1_direct, &signature_g1_affine));
    }
}

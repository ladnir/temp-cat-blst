// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter_common.h"

#include <memory>
#include <new>
#include <vector>

namespace catblst::fuzz::adapter
{
    namespace
    {
        struct Context final
        {
            std::unique_ptr<void, AlignedDelete> storage;
            std::vector<blst_byte> dst;

            explicit Context(std::span<const std::byte> dst_input)
                : storage(
                      ::operator new(
                          ::blst_pairing_sizeof(),
                          std::align_val_t(::blst_pairing_alignof())),
                      AlignedDelete{::blst_pairing_alignof()}),
                  dst(dst_input.size())
            {
                if (!dst_input.empty())
                    std::memcpy(dst.data(), dst_input.data(), dst.size());
                ::blst_pairing_init(
                    get(), true, dst.empty() ? nullptr : dst.data(), dst.size());
            }

            [[nodiscard]] blst_pairing* get() noexcept
            {
                return static_cast<blst_pairing*>(storage.get());
            }
        };

        [[nodiscard]] std::array<blst_byte, 32> ikm(
            std::span<const std::byte> input,
            unsigned char domain) noexcept
        {
            std::array<blst_byte, 32> result{};
            for (std::size_t index = 0; index < result.size(); ++index)
            {
                const auto source = input.empty()
                    ? static_cast<unsigned char>(index)
                    : std::to_integer<unsigned char>(
                          input[(index + domain) % input.size()]);
                result[index] = static_cast<blst_byte>(
                    source ^ domain ^ static_cast<unsigned char>(index * 29U + 1U));
            }
            return result;
        }

        [[nodiscard]] std::vector<std::byte> mutate(
            std::span<const std::byte> message,
            unsigned char selector)
        {
            std::vector<std::byte> result(message.begin(), message.end());
            if (result.empty())
            {
                result.push_back(std::byte{static_cast<unsigned char>(
                    selector == 0 ? 1U : selector)});
            }
            else
            {
                const std::size_t index = selector % result.size();
                result[index] ^= std::byte{static_cast<unsigned char>(
                    1U << (selector % 8U))};
            }
            return result;
        }

        void check_pk_in_g1(
            const std::array<blst_scalar, 2>& secrets,
            const std::array<std::span<const std::byte>, 2>& messages,
            std::span<const std::byte> dst,
            std::span<const std::byte> aug,
            unsigned char selector)
        {
            using detail::require;
            std::array<blst_p1, 2> public_keys{};
            std::array<blst_p2, 2> hashes{};
            std::array<blst_p2, 2> signatures{};
            std::array<blst_p1_affine, 2> public_keys_affine{};
            std::array<blst_p2_affine, 2> signatures_affine{};

            for (std::size_t index = 0; index < 2; ++index)
            {
                ::blst_sk_to_pk_in_g1(&public_keys[index], &secrets[index]);
                ::blst_hash_to_g2(
                    &hashes[index],
                    detail::bytes(messages[index]), messages[index].size(),
                    detail::bytes(dst), dst.size(),
                    detail::bytes(aug), aug.size());
                ::blst_sign_pk_in_g1(
                    &signatures[index], &hashes[index], &secrets[index]);
                ::blst_p1_to_affine(&public_keys_affine[index], &public_keys[index]);
                ::blst_p2_to_affine(&signatures_affine[index], &signatures[index]);
                require(::blst_core_verify_pk_in_g1(
                            &public_keys_affine[index], &signatures_affine[index],
                            true,
                            detail::bytes(messages[index]), messages[index].size(),
                            detail::bytes(dst), dst.size(),
                            detail::bytes(aug), aug.size()) == BLST_SUCCESS);

                std::array<blst_byte, 48> pk_wire{};
                std::array<blst_byte, 96> signature_wire{};
                ::blst_p1_compress(pk_wire.data(), &public_keys[index]);
                ::blst_p2_compress(signature_wire.data(), &signatures[index]);
                blst_p1_affine pk_roundtrip{};
                blst_p2_affine signature_roundtrip{};
                require(::blst_p1_uncompress(&pk_roundtrip, pk_wire.data()) ==
                        BLST_SUCCESS);
                require(::blst_p2_uncompress(
                            &signature_roundtrip, signature_wire.data()) ==
                        BLST_SUCCESS);
                require(::blst_p1_affine_is_equal(
                    &pk_roundtrip, &public_keys_affine[index]));
                require(::blst_p2_affine_is_equal(
                    &signature_roundtrip, &signatures_affine[index]));
            }

            Context aggregate(dst);
            for (std::size_t index = 0; index < 2; ++index)
            {
                require(::blst_pairing_aggregate_pk_in_g1(
                            aggregate.get(),
                            &public_keys_affine[index], &signatures_affine[index],
                            detail::bytes(messages[index]), messages[index].size(),
                            detail::bytes(aug), aug.size()) == BLST_SUCCESS);
            }
            ::blst_pairing_commit(aggregate.get());
            require(::blst_pairing_finalverify(
                aggregate.get(), nullptr)); // FUZZ_PATH: bls_api.aggregate_verify

            Context left(dst);
            Context right(dst);
            require(::blst_pairing_aggregate_pk_in_g1(
                        left.get(), &public_keys_affine[0], &signatures_affine[0],
                        detail::bytes(messages[0]), messages[0].size(),
                        detail::bytes(aug), aug.size()) == BLST_SUCCESS);
            require(::blst_pairing_aggregate_pk_in_g1(
                        right.get(), &public_keys_affine[1], &signatures_affine[1],
                        detail::bytes(messages[1]), messages[1].size(),
                        detail::bytes(aug), aug.size()) == BLST_SUCCESS);
            ::blst_pairing_commit(left.get());
            ::blst_pairing_commit(right.get());
            require(::blst_pairing_merge(left.get(), right.get()) == BLST_SUCCESS);
            require(::blst_pairing_finalverify(left.get(), nullptr));

            const auto altered = mutate(messages[1], selector);
            Context rejected(dst);
            require(::blst_pairing_aggregate_pk_in_g1(
                        rejected.get(), &public_keys_affine[0], &signatures_affine[0],
                        detail::bytes(messages[0]), messages[0].size(),
                        detail::bytes(aug), aug.size()) == BLST_SUCCESS);
            require(::blst_pairing_aggregate_pk_in_g1(
                        rejected.get(), &public_keys_affine[1], &signatures_affine[1],
                        detail::bytes(altered), altered.size(),
                        detail::bytes(aug), aug.size()) == BLST_SUCCESS);
            ::blst_pairing_commit(rejected.get());
            require(!::blst_pairing_finalverify(rejected.get(), nullptr));
            require(::blst_core_verify_pk_in_g1(
                        &public_keys_affine[1], &signatures_affine[1], true,
                        detail::bytes(altered), altered.size(),
                        detail::bytes(dst), dst.size(),
                        detail::bytes(aug), aug.size()) != BLST_SUCCESS);
        }

        void check_pk_in_g2(
            const std::array<blst_scalar, 2>& secrets,
            const std::array<std::span<const std::byte>, 2>& messages,
            std::span<const std::byte> dst,
            std::span<const std::byte> aug,
            unsigned char selector)
        {
            using detail::require;
            std::array<blst_p2, 2> public_keys{};
            std::array<blst_p1, 2> hashes{};
            std::array<blst_p1, 2> signatures{};
            std::array<blst_p2_affine, 2> public_keys_affine{};
            std::array<blst_p1_affine, 2> signatures_affine{};

            for (std::size_t index = 0; index < 2; ++index)
            {
                ::blst_sk_to_pk_in_g2(&public_keys[index], &secrets[index]);
                ::blst_hash_to_g1(
                    &hashes[index],
                    detail::bytes(messages[index]), messages[index].size(),
                    detail::bytes(dst), dst.size(),
                    detail::bytes(aug), aug.size());
                ::blst_sign_pk_in_g2(
                    &signatures[index], &hashes[index], &secrets[index]);
                ::blst_p2_to_affine(&public_keys_affine[index], &public_keys[index]);
                ::blst_p1_to_affine(&signatures_affine[index], &signatures[index]);
                require(::blst_core_verify_pk_in_g2(
                            &public_keys_affine[index], &signatures_affine[index],
                            true,
                            detail::bytes(messages[index]), messages[index].size(),
                            detail::bytes(dst), dst.size(),
                            detail::bytes(aug), aug.size()) == BLST_SUCCESS);

                std::array<blst_byte, 96> pk_wire{};
                std::array<blst_byte, 48> signature_wire{};
                ::blst_p2_compress(pk_wire.data(), &public_keys[index]);
                ::blst_p1_compress(signature_wire.data(), &signatures[index]);
                blst_p2_affine pk_roundtrip{};
                blst_p1_affine signature_roundtrip{};
                require(::blst_p2_uncompress(&pk_roundtrip, pk_wire.data()) ==
                        BLST_SUCCESS);
                require(::blst_p1_uncompress(
                            &signature_roundtrip, signature_wire.data()) ==
                        BLST_SUCCESS);
                require(::blst_p2_affine_is_equal(
                    &pk_roundtrip, &public_keys_affine[index]));
                require(::blst_p1_affine_is_equal(
                    &signature_roundtrip, &signatures_affine[index]));
            }

            Context aggregate(dst);
            for (std::size_t index = 0; index < 2; ++index)
            {
                require(::blst_pairing_aggregate_pk_in_g2(
                            aggregate.get(),
                            &public_keys_affine[index], &signatures_affine[index],
                            detail::bytes(messages[index]), messages[index].size(),
                            detail::bytes(aug), aug.size()) == BLST_SUCCESS);
            }
            ::blst_pairing_commit(aggregate.get());
            require(::blst_pairing_finalverify(aggregate.get(), nullptr));

            Context left(dst);
            Context right(dst);
            require(::blst_pairing_aggregate_pk_in_g2(
                        left.get(), &public_keys_affine[0], &signatures_affine[0],
                        detail::bytes(messages[0]), messages[0].size(),
                        detail::bytes(aug), aug.size()) == BLST_SUCCESS);
            require(::blst_pairing_aggregate_pk_in_g2(
                        right.get(), &public_keys_affine[1], &signatures_affine[1],
                        detail::bytes(messages[1]), messages[1].size(),
                        detail::bytes(aug), aug.size()) == BLST_SUCCESS);
            ::blst_pairing_commit(left.get());
            ::blst_pairing_commit(right.get());
            require(::blst_pairing_merge(left.get(), right.get()) == BLST_SUCCESS);
            require(::blst_pairing_finalverify(left.get(), nullptr));

            const auto altered = mutate(messages[0], selector);
            Context rejected(dst);
            require(::blst_pairing_aggregate_pk_in_g2(
                        rejected.get(), &public_keys_affine[0], &signatures_affine[0],
                        detail::bytes(altered), altered.size(),
                        detail::bytes(aug), aug.size()) == BLST_SUCCESS);
            require(::blst_pairing_aggregate_pk_in_g2(
                        rejected.get(), &public_keys_affine[1], &signatures_affine[1],
                        detail::bytes(messages[1]), messages[1].size(),
                        detail::bytes(aug), aug.size()) == BLST_SUCCESS);
            ::blst_pairing_commit(rejected.get());
            require(!::blst_pairing_finalverify(rejected.get(), nullptr));
            require(::blst_core_verify_pk_in_g2(
                        &public_keys_affine[0], &signatures_affine[0], true,
                        detail::bytes(altered), altered.size(),
                        detail::bytes(dst), dst.size(),
                        detail::bytes(aug), aug.size()) != BLST_SUCCESS);
        }
    }

    void run_bls_api_case(std::span<const std::byte> input)
    {
        input = input.first(std::min<std::size_t>(input.size(), 512));
        const unsigned char selector = input.empty()
            ? 0
            : std::to_integer<unsigned char>(input.front());
        const auto body = input.empty() ? input : input.subspan(1);
        const std::size_t split = body.size() / 2;
        const std::array<std::span<const std::byte>, 2> messages{
            body.first(split), body.subspan(split)};
        const auto dst = input.first(std::min<std::size_t>(input.size(), 64));
        const auto aug = input.last(std::min<std::size_t>(input.size(), 64));

        const auto ikm_a = ikm(input, 0x31);
        const auto ikm_b = ikm(input, 0xa7);
        std::array<blst_scalar, 2> secrets{};
        ::blst_keygen(&secrets[0], ikm_a.data(), ikm_a.size(), nullptr, 0);
        ::blst_keygen(&secrets[1], ikm_b.data(), ikm_b.size(), nullptr, 0);
        detail::require(::blst_sk_check(&secrets[0]));
        detail::require(::blst_sk_check(&secrets[1]));

        check_pk_in_g1(secrets, messages, dst, aug, selector);
        check_pk_in_g2(secrets, messages, dst, aug, selector);
    }
}

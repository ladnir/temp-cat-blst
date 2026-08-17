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
            std::size_t size = 0;

            Context(bool hash_or_encode, std::span<const std::byte> dst_input)
                : storage(
                      ::operator new(
                          ::blst_pairing_sizeof(),
                          std::align_val_t(::blst_pairing_alignof())),
                      AlignedDelete{::blst_pairing_alignof()}),
                  dst(dst_input.size()),
                  size(::blst_pairing_sizeof())
            {
                if (!dst_input.empty())
                    std::memcpy(dst.data(), dst_input.data(), dst.size());
                ::blst_pairing_init(
                    get(), hash_or_encode,
                    dst.empty() ? nullptr : dst.data(), dst.size());
            }

            [[nodiscard]] blst_pairing* get() noexcept
            {
                return static_cast<blst_pairing*>(storage.get());
            }

            [[nodiscard]] const blst_pairing* get() const noexcept
            {
                return static_cast<const blst_pairing*>(storage.get());
            }

            [[nodiscard]] std::vector<std::byte> snapshot() const
            {
                std::vector<std::byte> result(size);
                std::memcpy(result.data(), storage.get(), size);
                return result;
            }
        };
    }

    void run_pairing_state_case(std::span<const std::byte> input)
    {
        using detail::require;
        input = input.first(std::min<std::size_t>(input.size(), 1024));
        const auto dst = input.first(std::min<std::size_t>(input.size(), 512));
        const auto message = input;
        require(::blst_pairing_sizeof() != 0);
        require(::blst_pairing_alignof() != 0);

        const blst_p1_affine p = *::blst_p1_affine_generator();
        const blst_p2_affine q = *::blst_p2_affine_generator();
        blst_fp12 miller{};
        blst_fp12 via_lines{};
        ::blst_miller_loop(&miller, &q, &p);
        std::array<blst_fp6, 68> lines{};
        ::blst_precompute_lines(lines.data(), &q);
        ::blst_miller_loop_lines(&via_lines, lines.data(), &p);
        require(::blst_fp12_is_equal(&miller, &via_lines));
        require(::blst_fp12_finalverify(&miller, &via_lines));

        const std::array<blst_p2_affine, 2> q_values{q, q};
        const std::array<blst_p1_affine, 2> p_values{p, p};
        const std::array<const blst_p2_affine*, 2> qs{q_values.data(), nullptr};
        const std::array<const blst_p1_affine*, 2> ps{p_values.data(), nullptr};
        blst_fp12 miller_n{};
        blst_fp12 miller_squared{};
        ::blst_miller_loop_n(&miller_n, qs.data(), ps.data(), 2);
        ::blst_fp12_mul(&miller_squared, &miller, &miller);
        require(::blst_fp12_is_equal(&miller_n, &miller_squared));

        Context raw(false, dst);
        require(::blst_pairing_get_dst(raw.get()) ==
                (raw.dst.empty() ? nullptr : raw.dst.data()));
        ::blst_pairing_raw_aggregate(raw.get(), &q, &p);
        blst_fp12 raw_result{};
        ::blst_pairing_get_miller_result(&raw_result, raw.get());
        require(::blst_fp12_is_equal(&raw_result, &miller));

        std::array<blst_byte, 32> ikm{};
        for (std::size_t index = 0; index < ikm.size(); ++index)
            ikm[index] = input.empty()
                ? static_cast<blst_byte>(index + 1)
                : static_cast<blst_byte>(
                      std::to_integer<unsigned char>(input[index % input.size()]) +
                      static_cast<unsigned char>(index + 1));
        blst_scalar secret{};
        ::blst_keygen(&secret, ikm.data(), ikm.size(), nullptr, 0);
        require(::blst_sk_check(&secret));

        blst_p1 pk1_projective{};
        blst_p2 pk2_projective{};
        blst_p1 hash1{};
        blst_p2 hash2{};
        blst_p1 sig1_projective{};
        blst_p2 sig2_projective{};
        ::blst_sk_to_pk_in_g1(&pk1_projective, &secret);
        ::blst_sk_to_pk_in_g2(&pk2_projective, &secret);
        ::blst_hash_to_g1(
            &hash1, detail::bytes(message), message.size(),
            detail::bytes(dst), dst.size(), nullptr, 0);
        ::blst_hash_to_g2(
            &hash2, detail::bytes(message), message.size(),
            detail::bytes(dst), dst.size(), nullptr, 0);
        ::blst_sign_pk_in_g2(&sig1_projective, &hash1, &secret);
        ::blst_sign_pk_in_g1(&sig2_projective, &hash2, &secret);
        blst_p1_affine pk1{};
        blst_p2_affine pk2{};
        blst_p1_affine sig1{};
        blst_p2_affine sig2{};
        ::blst_p1_to_affine(&pk1, &pk1_projective);
        ::blst_p2_to_affine(&pk2, &pk2_projective);
        ::blst_p1_to_affine(&sig1, &sig1_projective);
        ::blst_p2_to_affine(&sig2, &sig2_projective);

        Context g1(true, dst);
        require(::blst_pairing_aggregate_pk_in_g1(
                    g1.get(), &pk1, &sig2,
                    detail::bytes(message), message.size(), nullptr, 0) ==
                BLST_SUCCESS);
        ::blst_pairing_commit(g1.get());
        require(::blst_pairing_finalverify(g1.get(), nullptr));

        Context g2(true, dst);
        require(::blst_pairing_aggregate_pk_in_g2(
                    g2.get(), &pk2, &sig1,
                    detail::bytes(message), message.size(), nullptr, 0) ==
                BLST_SUCCESS);
        ::blst_pairing_commit(g2.get());
        require(::blst_pairing_finalverify(g2.get(), nullptr));

        const auto scalar = detail::fixed_bytes<32>(input);
        const std::size_t nbits = input.empty()
            ? 0
            : std::to_integer<unsigned char>(input.front()) % 257U;
        const std::size_t weight_nbits = std::max<std::size_t>(nbits, 1);
        const blst_byte* scalar_pointer = scalar.data();

        Context variants_g1(true, dst);
        require(::blst_pairing_chk_n_aggr_pk_in_g1(
                    variants_g1.get(), &pk1, true, &sig2, true,
                    detail::bytes(message), message.size(), nullptr, 0) ==
                BLST_SUCCESS);
        require(::blst_pairing_mul_n_aggregate_pk_in_g1(
                    variants_g1.get(), &pk1, &sig2,
                    scalar_pointer, weight_nbits,
                    detail::bytes(message), message.size(), nullptr, 0) ==
                BLST_SUCCESS);
        require(::blst_pairing_chk_n_mul_n_aggr_pk_in_g1(
                    variants_g1.get(), &pk1, true, &sig2, true,
                    scalar_pointer, weight_nbits,
                    detail::bytes(message), message.size(), nullptr, 0) ==
                BLST_SUCCESS);
        ::blst_pairing_commit(variants_g1.get());

        Context variants_g2(true, dst);
        require(::blst_pairing_chk_n_aggr_pk_in_g2(
                    variants_g2.get(), &pk2, true, &sig1, true,
                    detail::bytes(message), message.size(), nullptr, 0) ==
                BLST_SUCCESS);
        require(::blst_pairing_mul_n_aggregate_pk_in_g2(
                    variants_g2.get(), &pk2, &sig1,
                    scalar_pointer, weight_nbits,
                    detail::bytes(message), message.size(), nullptr, 0) ==
                BLST_SUCCESS);
        require(::blst_pairing_chk_n_mul_n_aggr_pk_in_g2(
                    variants_g2.get(), &pk2, true, &sig1, true,
                    scalar_pointer, weight_nbits,
                    detail::bytes(message), message.size(), nullptr, 0) ==
                BLST_SUCCESS);
        ::blst_pairing_commit(variants_g2.get());

        Context external(true, dst);
        require(::blst_pairing_aggregate_pk_in_g1(
                    external.get(), &pk1, nullptr,
                    detail::bytes(message), message.size(), nullptr, 0) ==
                BLST_SUCCESS);
        ::blst_pairing_commit(external.get());
        blst_fp12 external_signature{};
        ::blst_aggregated_in_g2(&external_signature, &sig2);
        require(::blst_pairing_finalverify(external.get(), &external_signature));

        Context unset_a(true, dst);
        Context unset_b(true, dst);
        ::blst_pairing_commit(unset_a.get());
        ::blst_pairing_commit(unset_b.get());
        require(::blst_pairing_merge(unset_a.get(), unset_b.get()) == BLST_SUCCESS);

        Context mismatch(true, dst);
        require(::blst_pairing_aggregate_pk_in_g1(
                    mismatch.get(), nullptr, nullptr,
                    nullptr, 0, nullptr, 0) == BLST_SUCCESS);
        const auto before = mismatch.snapshot();
        require(::blst_pairing_aggregate_pk_in_g2(
                    mismatch.get(), nullptr, nullptr,
                    nullptr, 0, nullptr, 0) == BLST_AGGR_TYPE_MISMATCH);
        const bool unchanged =
            std::memcmp(before.data(), mismatch.storage.get(), mismatch.size) == 0;
        require(unchanged); // FUZZ_PATH: pairing.type_mismatch_atomic
    }
}

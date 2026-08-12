/*
 * Copyright Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include "blst/common/defines.h"
#include "blst/fields/Fr.h"

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace blst
{
    namespace fft_detail
    {
        inline constexpr std::size_t fr_two_adicity = 32;

        // Primitive 2^k-th roots in Montgomery form, indexed by k. The table
        // is derived from 7^((r - 1) / 2^32) and repeated squaring. Keep the
        // independent scripts/audit/validate_fft_roots.py check in sync.
        inline constexpr std::array<Fr, fr_two_adicity + 1> fr_power_of_two_roots = {
            Fr{{ 0x00000001fffffffeULL, 0x5884b7fa00034802ULL, 0x998c4fefecbc4ff5ULL, 0x1824b159acc5056fULL }}, // 2^0
            Fr{{ 0xfffffffd00000003ULL, 0xfb38ec08fffb13fcULL, 0x99ad88181ce5880fULL, 0x5bc8f5f97cd877d8ULL }}, // 2^1
            Fr{{ 0xf3b05674aa89cfb1ULL, 0x072f01406006b9feULL, 0xce9a0dbf25667a26ULL, 0x4d2ce4052d598374ULL }}, // 2^2
            Fr{{ 0x6ecf49952c42eff8ULL, 0x6ec1fa834d5dfa57ULL, 0x3ff5cf385656b15cULL, 0x54b80267229e031bULL }}, // 2^3
            Fr{{ 0x1acf64e488202b65ULL, 0x5b18b637819f004fULL, 0xee649cab3107fdf5ULL, 0x27a8543dce9012a0ULL }}, // 2^4
            Fr{{ 0x30c78850dd04d52fULL, 0xe27f8794adaddbcbULL, 0x0a59f321ae22c531ULL, 0x1958b127ba7f8ca3ULL }}, // 2^5
            Fr{{ 0xbcf7e31ae8d6f162ULL, 0x15d1c587662fc7a4ULL, 0xb7ed5ecd71278a57ULL, 0x329d0d9340ccecb3ULL }}, // 2^6
            Fr{{ 0xdaf6401c4aed2336ULL, 0xa14c0cba41de4f6aULL, 0xdfdc7d04fc7f2e5bULL, 0x6f61ff4e296cff18ULL }}, // 2^7
            Fr{{ 0xb4b87dd6a9d4956cULL, 0x6a51d0bcd1482f03ULL, 0x7e8146f7e0dbcb1eULL, 0x21f4b238cccb85aeULL }}, // 2^8
            Fr{{ 0x1ff70724b651888bULL, 0xcd8166cbe56a5af4ULL, 0xe4a992a398e1e256ULL, 0x5705ea3686fd3b03ULL }}, // 2^9
            Fr{{ 0xa4a6e37a36d7822aULL, 0xfc83a2ce8787e5f1ULL, 0xe89553dae7ea8f80ULL, 0x6f6d38bb2a92ff1aULL }}, // 2^10
            Fr{{ 0xf4a46fbb8dc619e5ULL, 0x82fa97dcb8e77487ULL, 0xaf909c5bf198f018ULL, 0x00edba7d83a038faULL }}, // 2^11
            Fr{{ 0xf2df262c09458a39ULL, 0x048cdf5b99dff177ULL, 0x16857bc5c7cce57bULL, 0x043b3dbca4a915aeULL }}, // 2^12
            Fr{{ 0xa33d279ff0ccffc9ULL, 0x41fac79f59e91972ULL, 0x065d227fead1139bULL, 0x71db41abda03e055ULL }}, // 2^13
            Fr{{ 0x9e468a30f2a9a5d6ULL, 0x1fe34726ca6bcf1bULL, 0x54a621ae87eb73acULL, 0x100a9d2ebfd61a47ULL }}, // 2^14
            Fr{{ 0x6ebe497bf4e11f70ULL, 0x1533f4823c55cafcULL, 0x66a7f2da7eb4d017ULL, 0x5b6be79ab07f3e81ULL }}, // 2^15
            Fr{{ 0x4526e48b11c947f9ULL, 0xccc36140c614d3aaULL, 0xf2e56f1a03cd5f24ULL, 0x57db5a4c159e0248ULL }}, // 2^16
            Fr{{ 0x860cb52e2d964320ULL, 0x4c5059aa9678fcc9ULL, 0xd92f98d7cfe87e16ULL, 0x40329fd2d283e2d3ULL }}, // 2^17
            Fr{{ 0x56e7145ccb39b4a6ULL, 0x79825deadcbe3defULL, 0x3e516f05f02abbe0ULL, 0x532595f8bef5bd2dULL }}, // 2^18
            Fr{{ 0x3b801e7d90f9b375ULL, 0xaada3338524be4d6ULL, 0x5f0b11a6bdbd3f25ULL, 0x04063ff85adb38c4ULL }}, // 2^19
            Fr{{ 0x146391f8784685ceULL, 0x6975fe2172dbcd2cULL, 0x7bfaf2b2edc848d6ULL, 0x13c5c6b3afb23c0cULL }}, // 2^20
            Fr{{ 0x245f9cbc1f88b5f3ULL, 0x5dfbb5c68ef4a333ULL, 0xa9ac7cdbc68f208fULL, 0x44448a9838255605ULL }}, // 2^21
            Fr{{ 0x0db3e6ccbf7c975bULL, 0x9d463b98af4d41ecULL, 0x555e1ff59ab7a3e7ULL, 0x6154be202a31f630ULL }}, // 2^22
            Fr{{ 0xbcd2a6e09538df40ULL, 0x290e2f0c95e6fafdULL, 0x37b754afe7fe2f1dULL, 0x500b5365d11e6548ULL }}, // 2^23
            Fr{{ 0x6a4fac5435b444faULL, 0xa4c9d25df8a367daULL, 0xa33932bd4c0632ffULL, 0x038a3b43cfb507b5ULL }}, // 2^24
            Fr{{ 0x25acfadbd37328caULL, 0xa2157fa336ffb229ULL, 0x432bc6f21ac4d49eULL, 0x60eb5b8db0894969ULL }}, // 2^25
            Fr{{ 0x8a0719bd9fc1ac27ULL, 0x744144bb29716a74ULL, 0xffea9e0b541fcb9aULL, 0x4459ada388bb635fULL }}, // 2^26
            Fr{{ 0xe7b5ec2bfbd35ce8ULL, 0xa99a0f8a35e37e7bULL, 0x312c2051a70cacc4ULL, 0x415444dfc9407c7dULL }}, // 2^27
            Fr{{ 0x2c75ff0b8d27bcb7ULL, 0x00fc392e4413b9d7ULL, 0xef31d1fd8752020fULL, 0x6d5074c8e7a75dc5ULL }}, // 2^28
            Fr{{ 0x250a9d0e4b092216ULL, 0x525a653781be1b5eULL, 0xe4381082889cdfb0ULL, 0x02536cde1d9f1a1cULL }}, // 2^29
            Fr{{ 0x8ea72d687bb12408ULL, 0xa9ed127cb02387bcULL, 0x9c59d748e7669430ULL, 0x4e720080e9ced390ULL }}, // 2^30
            Fr{{ 0x1f2fda6f3141e583ULL, 0x8d44592227fe5b46ULL, 0xe33cee1b81e59a5fULL, 0x22c666251ad8ae31ULL }}, // 2^31
            Fr{{ 0xb9b58d8c5f0e466aULL, 0x5b1b4c801819d7ecULL, 0x0af53ae352a31e64ULL, 0x5bf3adda19e9b27bULL }}, // 2^32
        };

        template<typename>
        inline constexpr bool dependent_false = false;
    }

    /**
     * Return whether @p root has exact power-of-two multiplicative order @p n.
     * Invalid orders, zero, and roots of a proper subgroup return false.
     */
    template<typename F>
    [[nodiscard]] inline bool isPrimRootOfUnity(const u64 n, const F& root)
    {
        if (!std::has_single_bit(n) || root == F(0)) {
            return false;
        }
        if (root.pow(n) != F(1)) {
            return false;
        }
        return n == 1 || root.pow(n / 2) != F(1);
    }

    template<>
    [[nodiscard]] inline bool isPrimRootOfUnity<Fr>(const u64 n, const Fr& root)
    {
        if (!std::has_single_bit(n) || n > (u64{1} << fft_detail::fr_two_adicity)) {
            return false;
        }
        if (n == 1) {
            return root == Fr::one();
        }

        Fr power = root;
        const auto log_order = std::countr_zero(n);
        for (int i = 1; i < log_order; ++i) {
            Fr::sqr(power, power);
        }
        if (power == Fr::one()) {
            return false;
        }
        Fr::sqr(power, power);
        return power == Fr::one();
    }

    /**
     * Derive a primitive power-of-two @p n-th root from a caller-supplied
     * generator. The field order and exponent are public. Throws
     * std::invalid_argument before returning a value with the wrong order.
     */
    template<typename F>
    [[nodiscard]] inline F primRootOfUnity(const u64 n, const F& generator)
    {
        if (!std::has_single_bit(n)) {
            throw std::invalid_argument("primRootOfUnity: n must be a power of two");
        }

        const auto field_order_minus_one = F::order() - 1;
        if (field_order_minus_one % n != 0) {
            throw std::invalid_argument("primRootOfUnity: n does not divide the field order minus one");
        }

        F root = generator.pow(field_order_minus_one / n);
        if (!isPrimRootOfUnity(n, root)) {
            throw std::invalid_argument("primRootOfUnity: generator does not produce an exact-order root");
        }
        return root;
    }

    // Only Fr has a built-in generator/root table. Other fields must use the
    // explicit-generator overload so root construction is deterministic.
    template<typename F>
    [[nodiscard]] inline F primRootOfUnity(const u64)
    {
        static_assert(fft_detail::dependent_false<F>,
            "primRootOfUnity<F>(n) is available only for Fr; supply a generator");
        return F{};
    }

    template<>
    [[nodiscard]] inline Fr primRootOfUnity<Fr>(const u64 n)
    {
        if (!std::has_single_bit(n)) {
            throw std::invalid_argument("primRootOfUnity: n must be a power of two");
        }
        if (n > (u64{1} << fft_detail::fr_two_adicity)) {
            throw std::invalid_argument("primRootOfUnity: n exceeds Fr two-adicity");
        }
        return fft_detail::fr_power_of_two_roots[std::countr_zero(n)];
    }

    /**
     * Fill @p roots with consecutive powers of the exact primitive root @p psi.
     *
     * @pre The span size is a supported power of two and @p psi has that exact
     * multiplicative order. Invalid inputs throw std::invalid_argument before
     * any output element is modified. The input and output may not overlap.
     */
    template<typename F>
    inline void precomputeRootsOfUnity(const F& psi, const std::span<F> roots)
    {
        const std::size_t n = roots.size();
        if (!std::has_single_bit(n)) {
            throw std::invalid_argument("precomputeRootsOfUnity: root count must be a nonzero power of two");
        }

        const auto field_order_minus_one = F::order() - 1;
        if (n > field_order_minus_one || field_order_minus_one % n != 0) {
            throw std::invalid_argument("precomputeRootsOfUnity: root count is unsupported by the field");
        }
        if (!isPrimRootOfUnity(static_cast<u64>(n), psi)) {
            throw std::invalid_argument("precomputeRootsOfUnity: psi does not have the required exact order");
        }

        F root = F(1);
        for (F& destination : roots) {
            destination = root;
            root *= psi;
        }
    }
}

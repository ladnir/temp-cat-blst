/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_SCALAR_HASH_TO_FIELD_H
#define BLST_SCALAR_HASH_TO_FIELD_H

#include "blst/sign/sha256.h"
#include "blst/scalar/fields/Fp2.h"

namespace blst::scalar
{
    inline constexpr size_t xmd_sha256_block_bytes = 32;
    inline constexpr size_t xmd_sha256_max_bytes = 255 * xmd_sha256_block_bytes;

    namespace hash_to_field_detail
    {
        BLST_HD BLST_INLINE void zero_bytes(unsigned char* bytes, size_t len) noexcept
        {
            vec_zero(bytes, len);
        }

        BLST_HD BLST_INLINE void sha256_init_Zpad(Sha256Context* ctx) noexcept
        {
            ctx->h[0] = 0xda5698beU;
            ctx->h[1] = 0x17b9b469U;
            ctx->h[2] = 0x62335799U;
            ctx->h[3] = 0x779fbecaU;
            ctx->h[4] = 0x8ce5d491U;
            ctx->h[5] = 0xc0d26243U;
            ctx->h[6] = 0xbafef9eaU;
            ctx->h[7] = 0x1837a9d8U;
            ctx->total_bytes = 64;
            sha256_zero_buffer(ctx->buffer);
            ctx->buffered_bytes = 0;
        }

        BLST_HD BLST_FORCEINLINE void xor_32(unsigned char out[32],
            const unsigned char a[32], const unsigned char b[32]) noexcept
        {
            for (size_t i = 0; i < 32; ++i)
                out[i] = a[i] ^ b[i];
        }

        template<typename OutputField, size_t count>
        BLST_HD BLST_FORCEINLINE Fp& component(
            OutputField (&elems)[count], size_t index) noexcept
        {
            static_assert(std::is_same_v<OutputField, Fp> ||
                std::is_same_v<OutputField, Fp2>);

            if constexpr (std::is_same_v<OutputField, Fp>)
                return elems[index];
            else
                return elems[index / 2][index % 2];
        }
    }

#ifdef expand_message_xmd
    /*
     * Optional upstream-compatible override. The implementation must emit
     * exactly len_in_bytes bytes, accept lengths through
     * xmd_sha256_max_bytes, allocate no storage, and not throw.
     */
    BLST_HD void expand_message_xmd(unsigned char* bytes, size_t len_in_bytes,
        const unsigned char* aug, size_t aug_len,
        const unsigned char* msg, size_t msg_len,
        const unsigned char* DST, size_t DST_len) noexcept;
#else
    BLST_HD BLST_INLINE void expand_message_xmd(unsigned char* bytes, size_t len_in_bytes,
        const unsigned char* aug, size_t aug_len,
        const unsigned char* msg, size_t msg_len,
        const unsigned char* DST, size_t DST_len) noexcept
    {
        if (len_in_bytes == 0)
            return;

        BLST_ASSERT(bytes != nullptr);
        BLST_ASSERT(aug != nullptr || aug_len == 0);
        BLST_ASSERT(msg != nullptr || msg_len == 0);
        BLST_ASSERT(DST != nullptr || DST_len == 0);
        BLST_ASSERT(len_in_bytes <= xmd_sha256_max_bytes);
        if (bytes == nullptr || len_in_bytes > xmd_sha256_max_bytes ||
            (aug == nullptr && aug_len != 0) ||
            (msg == nullptr && msg_len != 0) ||
            (DST == nullptr && DST_len != 0)) {
            if (bytes != nullptr && len_in_bytes <= xmd_sha256_max_bytes)
                hash_to_field_detail::zero_bytes(bytes, len_in_bytes);
            return;
        }

        alignas(limb_t) unsigned char b_0[32];
        alignas(limb_t) unsigned char b_i[33 + 256 + 31];
        alignas(limb_t) unsigned char final_block[32];
        unsigned char* p;
        size_t i, b_i_bits, b_i_blocks;
        Sha256Context ctx;

        /*
         * compose template for 'strxor(b_0, b_(i-1)) || I2OSP(i, 1) || DST_prime'
         */
        if (DST_len > 255) {
            sha256_init(&ctx);
            sha256_update(&ctx, "H2C-OVERSIZE-DST-", 17);
            sha256_update(&ctx, DST, DST_len);
            sha256_final(b_0, &ctx);
            DST = b_0, DST_len = 32;
        }
        b_i_blocks = ((33 + DST_len + 1 + 9) + 63) & (size_t)-64;
        vec_zero(b_i + b_i_blocks - 64, 64);

        p = b_i + 33;
        for (i = 0; i < DST_len; i++)
            p[i] = DST[i];
        p[i++] = (unsigned char)DST_len;
        p[i++] = 0x80;
        p[i + 6] = p[i + 5] = p[i + 4] = p[i + 3] = p[i + 2] = p[i + 1] = p[i + 0] = 0;
        b_i_bits = (33 + DST_len + 1) * 8;
        p = b_i + b_i_blocks;
        p[-2] = (unsigned char)(b_i_bits >> 8);
        p[-1] = (unsigned char)(b_i_bits);

        hash_to_field_detail::sha256_init_Zpad(&ctx);   /* Z_pad | */
        sha256_update(&ctx, aug, aug_len);              /* | aug | */
        sha256_update(&ctx, msg, msg_len);              /* | msg | */
        /* | I2OSP(len_in_bytes, 2) || I2OSP(0, 1) || DST_prime    */
        b_i[30] = (unsigned char)(len_in_bytes >> 8);
        b_i[31] = (unsigned char)(len_in_bytes);
        b_i[32] = 0;
        sha256_update(&ctx, b_i + 30, 3 + DST_len + 1);
        sha256_final(b_0, &ctx);

        sha256_init_h(ctx.h);
        vec_copy(b_i, b_0, 32);
        ++b_i[32];
        blst_sha256_block_data_order_v(ctx.h, b_i, b_i_blocks / 64);

        const size_t ell = (len_in_bytes + 31) / 32;
        size_t remaining = len_in_bytes;
        unsigned char* output = bytes;
        const unsigned char* previous;
        if (remaining >= 32) {
            blst_sha256_emit_v(output, ctx.h);
            previous = output;
            output += 32;
            remaining -= 32;
        }
        else {
            blst_sha256_emit_v(final_block, ctx.h);
            blst_sha256_bcopy_v(output, final_block, remaining);
            return;
        }

        for (size_t block = 1; block < ell; ++block) {
            sha256_init_h(ctx.h);
            hash_to_field_detail::xor_32(b_i, b_0, previous);
            ++b_i[32];
            blst_sha256_block_data_order_v(ctx.h, b_i, b_i_blocks / 64);

            if (remaining >= 32) {
                blst_sha256_emit_v(output, ctx.h);
                previous = output;
                output += 32;
                remaining -= 32;
            }
            else {
                blst_sha256_emit_v(final_block, ctx.h);
                blst_sha256_bcopy_v(output, final_block, remaining);
                remaining = 0;
            }
        }
        BLST_ASSERT(remaining == 0);
    }
#endif

    /**
     * |nelems| is 'count * m' from spec. OutputField is Fp for m=1 and Fp2
     * for m=2; each write names its actual field subobject rather than
     * flattening adjacent Fp2 objects through pointer arithmetic. msg, DST,
     * and aug may be null exactly when their corresponding length is zero.
     * The result consists of reduced Montgomery-form field elements.
     */
	template<size_t nelems, typename OutputField, size_t count>
    BLST_HD BLST_INLINE void hash_to_field(OutputField (&elems)[count],
        const unsigned char* aug, size_t aug_len,
        const unsigned char* msg, size_t msg_len,
        const unsigned char* DST, size_t DST_len) noexcept
    {
        constexpr size_t extension_degree =
            std::is_same_v<OutputField, Fp2> ? 2 : 1;
        static_assert(nelems == count * extension_degree);

        constexpr size_t L = sizeof(Fp) + 128 / 8;  /* ceil((ceil(log2(p)) + k) / 8) */
        constexpr size_t len_in_bytes = L * nelems;   /* divisible by 64, hurray!      */
        static_assert(len_in_bytes % sizeof(limb_t) == 0);
        static_assert(len_in_bytes <= xmd_sha256_max_bytes,
            "hash_to_field exceeds the RFC 9380 XMD:SHA-256 output domain");

        limb_t pseudo_random[len_in_bytes / sizeof(limb_t)];
        unsigned char* bytes;
        UFp elem;

        aug_len = aug != NULL ? aug_len : 0;
        DST_len = DST != NULL ? DST_len : 0;

        expand_message_xmd((unsigned char*)pseudo_random, len_in_bytes,
            aug, aug_len, msg, msg_len, DST, DST_len);

        vec_zero(elem, sizeof(elem));
        bytes = (unsigned char*)pseudo_random;
        size_t output_index = 0;
		size_t nn = nelems;
        while (nn--) {
            limbs_from_be_bytes(elem, bytes, L);
            bytes += L;
            /*
             * L-bytes block % P, output is in Montgomery domain...
             */
            Fp& output = hash_to_field_detail::component(elems, output_index++);
            Fp::redc(output, elem);
            Fp::mul(output, output, get_BLST_RRRR());
        }
    }

    /**
     * Expands msg with RFC 9380 XMD:SHA-256 and writes exactly len_in_bytes.
     * The output is allocation-free. Invalid nonzero pointer/length pairs or
     * lengths above xmd_sha256_max_bytes clear a valid output buffer. A zero
     * output length is a no-op and permits bytes to be null.
     */
    BLST_HD BLST_INLINE void blst_expand_message_xmd(unsigned char* bytes, size_t len_in_bytes,
        const unsigned char* msg, size_t msg_len,
        const unsigned char* DST, size_t DST_len) noexcept
    {
        if (len_in_bytes == 0)
            return;

        BLST_ASSERT(bytes != nullptr);
        BLST_ASSERT(msg != nullptr || msg_len == 0);
        BLST_ASSERT(DST != nullptr || DST_len == 0);
        if (bytes == nullptr)
            return;

        if (len_in_bytes > xmd_sha256_max_bytes ||
            (msg == nullptr && msg_len != 0) ||
            (DST == nullptr && DST_len != 0)) {
            hash_to_field_detail::zero_bytes(bytes, len_in_bytes);
            return;
        }

        expand_message_xmd(bytes, len_in_bytes, nullptr, 0, msg, msg_len,
            DST, DST_len);
    }
}

#endif

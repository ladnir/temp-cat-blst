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
    BLST_HD BLST_INLINE void zero_bytes(unsigned char* bytes, size_t len)
    {
        vec_zero(bytes, len);
    }

#ifdef expand_message_xmd
    BLST_HD void expand_message_xmd(unsigned char* bytes, size_t len_in_bytes,
        const unsigned char* aug, size_t aug_len,
        const unsigned char* msg, size_t msg_len,
        const unsigned char* DST, size_t DST_len);
#else
    BLST_HD BLST_INLINE void sha256_init_Zpad(Sha256Context* ctx)
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

    BLST_HD BLST_INLINE void vec_xor(void* BLST_restrict ret, const void* BLST_restrict a,
        const void* BLST_restrict b, size_t num)
    {
        limb_t* rp = (limb_t*)ret;
        const limb_t* ap = (const limb_t*)a;
        const limb_t* bp = (const limb_t*)b;
        size_t i;

        num /= sizeof(limb_t);

        for (i = 0; i < num; i++)
            rp[i] = ap[i] ^ bp[i];
    }

    BLST_HD BLST_INLINE void expand_message_xmd(unsigned char* bytes, size_t len_in_bytes,
        const unsigned char* aug, size_t aug_len,
        const unsigned char* msg, size_t msg_len,
        const unsigned char* DST, size_t DST_len)
    {
        union { limb_t align; unsigned char c[32]; } b_0;
        union { limb_t align; unsigned char c[33 + 256 + 31]; } b_i;
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
            sha256_final(b_0.c, &ctx);
            DST = b_0.c, DST_len = 32;
        }
        b_i_blocks = ((33 + DST_len + 1 + 9) + 63) & (size_t)-64;
        vec_zero(b_i.c + b_i_blocks - 64, 64);

        p = b_i.c + 33;
        for (i = 0; i < DST_len; i++)
            p[i] = DST[i];
        p[i++] = (unsigned char)DST_len;
        p[i++] = 0x80;
        p[i + 6] = p[i + 5] = p[i + 4] = p[i + 3] = p[i + 2] = p[i + 1] = p[i + 0] = 0;
        b_i_bits = (33 + DST_len + 1) * 8;
        p = b_i.c + b_i_blocks;
        p[-2] = (unsigned char)(b_i_bits >> 8);
        p[-1] = (unsigned char)(b_i_bits);

        sha256_init_Zpad(&ctx);                         /* Z_pad | */
        sha256_update(&ctx, aug, aug_len);              /* | aug | */
        sha256_update(&ctx, msg, msg_len);              /* | msg | */
        /* | I2OSP(len_in_bytes, 2) || I2OSP(0, 1) || DST_prime    */
        b_i.c[30] = (unsigned char)(len_in_bytes >> 8);
        b_i.c[31] = (unsigned char)(len_in_bytes);
        b_i.c[32] = 0;
        sha256_update(&ctx, b_i.c + 30, 3 + DST_len + 1);
        sha256_final(b_0.c, &ctx);

        sha256_init_h(ctx.h);
        vec_copy(b_i.c, b_0.c, 32);
        ++b_i.c[32];
        blst_sha256_block_data_order_v(ctx.h, b_i.c, b_i_blocks / 64);
        blst_sha256_emit_v(bytes, ctx.h);

        len_in_bytes += 31; /* ell = ceil(len_in_bytes / b_in_bytes), with */
        len_in_bytes /= 32; /* caller being responsible for accordingly large
                             * buffer. hash_to_field passes one with length
                             * divisible by 64, remember? which works... */
        while (--len_in_bytes) {
            sha256_init_h(ctx.h);
            vec_xor(b_i.c, b_0.c, bytes, 32);
            bytes += 32;
            ++b_i.c[32];
            blst_sha256_block_data_order_v(ctx.h, b_i.c, b_i_blocks / 64);
            blst_sha256_emit_v(bytes, ctx.h);
        }
    }
#endif

    template<typename OutputField, size_t count>
    BLST_HD BLST_FORCEINLINE Fp& hash_to_field_component(
        OutputField (&elems)[count], size_t index)
    {
        static_assert(std::is_same_v<OutputField, Fp> ||
            std::is_same_v<OutputField, Fp2>);

        if constexpr (std::is_same_v<OutputField, Fp>)
            return elems[index];
        else
            return elems[index / 2][index % 2];
    }

    /*
     * |nelems| is 'count * m' from spec. OutputField is Fp for m=1 and Fp2
     * for m=2; each write names its actual field subobject rather than
     * flattening adjacent Fp2 objects through pointer arithmetic.
     */
	template<size_t nelems, typename OutputField, size_t count>
    BLST_HD BLST_INLINE void hash_to_field(OutputField (&elems)[count],
        const unsigned char* aug, size_t aug_len,
        const unsigned char* msg, size_t msg_len,
        const unsigned char* DST, size_t DST_len)
    {
        constexpr size_t extension_degree =
            std::is_same_v<OutputField, Fp2> ? 2 : 1;
        static_assert(nelems == count * extension_degree);

        constexpr size_t L = sizeof(Fp) + 128 / 8;  /* ceil((ceil(log2(p)) + k) / 8) */
        constexpr size_t len_in_bytes = L * nelems;   /* divisible by 64, hurray!      */
        static_assert(len_in_bytes % sizeof(limb_t) == 0);

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
            Fp& output = hash_to_field_component(elems, output_index++);
            Fp::redc(output, elem);
            Fp::mul(output, output, get_BLST_RRRR());
        }
    }

    BLST_HD BLST_INLINE void blst_expand_message_xmd(unsigned char* bytes, size_t len_in_bytes,
        const unsigned char* msg, size_t msg_len,
        const unsigned char* DST, size_t DST_len)
    {
        if (len_in_bytes == 0 || bytes == nullptr) {
            return;
        }

        if (len_in_bytes > static_cast<size_t>(-1) - 31) {
            zero_bytes(bytes, len_in_bytes);
            return;
        }

        size_t buf_len = (len_in_bytes + 31) & ((size_t)0 - 32);
        unsigned char* buf_ptr = bytes;

        if (buf_len > 255 * 32) {
            zero_bytes(bytes, len_in_bytes);
            return;
        }

        if (buf_len != len_in_bytes)
        {
			buf_ptr = (unsigned char*)malloc(buf_len);
            if (buf_ptr == nullptr) {
                zero_bytes(bytes, len_in_bytes);
                return;
            }
        }

        expand_message_xmd(buf_ptr, len_in_bytes, NULL, 0, msg, msg_len,
            DST, DST_len);
        if (buf_ptr != bytes) {
            unsigned char* ptr = buf_ptr;
            while (len_in_bytes--)
                *bytes++ = *ptr++;
            vec_zero(buf_ptr, buf_len);
			free(buf_ptr);
        }
    }
}

#endif

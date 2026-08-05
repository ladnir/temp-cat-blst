/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/common/bytes.h"
#include "blst/sign/sha256.h"
#include "blst/backend/backend.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fr.h"

//BLST_EXTERN_C_BEGIN

namespace blst
{
    namespace
    {
    struct HmacSha256Context
    {
        Sha256Context ctx;
        unsigned int h_ipad[8];
        unsigned int h_opad[8];
        /* Active limbs; unsigned-byte access names their object representation. */
        limb_t tail[64 / sizeof(limb_t)];
    };

    static_assert(64 % sizeof(limb_t) == 0);

    BLST_FORCEINLINE unsigned char* hmac_tail_bytes(
        HmacSha256Context* ctx) noexcept
    {
        return reinterpret_cast<unsigned char*>(ctx->tail);
    }

    BLST_FORCEINLINE void hmac_zero_tail(HmacSha256Context* ctx) noexcept
    {
        /* Functional initialization; the owning keygen scratch is erased later. */
        for (size_t i = 0; i < sizeof(ctx->tail) / sizeof(ctx->tail[0]); ++i)
            ctx->tail[i] = 0;
    }

    void hmac_init(HmacSha256Context* ctx, const void* K,
        size_t K_len) noexcept
    {
        BLST_ASSERT(ctx != nullptr);
        BLST_ASSERT(K != nullptr || K_len == 0);

        if (K == NULL) {            /* reuse h_ipad and h_opad */
            blst_sha256_hcopy_v(ctx->ctx.h, ctx->h_ipad);
            ctx->ctx.total_bytes = 64;
            sha256_zero_buffer(ctx->ctx.buffer);
            ctx->ctx.buffered_bytes = 0;

            return;
        }

        unsigned char* tail = hmac_tail_bytes(ctx);
        hmac_zero_tail(ctx);
        if (K_len > 64) {
            sha256_init(&ctx->ctx);
            sha256_update(&ctx->ctx, K, K_len);
            sha256_final(tail, &ctx->ctx);
        }
        else if (K_len != 0) {
            blst_sha256_bcopy_v(tail, K, K_len);
        }

        for (size_t i = 0; i < sizeof(ctx->tail) / sizeof(ctx->tail[0]); ++i)
            ctx->tail[i] ^= (limb_t)0x3636363636363636;

        sha256_init(&ctx->ctx);
        sha256_update(&ctx->ctx, tail, 64);
        blst_sha256_hcopy_v(ctx->h_ipad, ctx->ctx.h);

        for (size_t i = 0; i < sizeof(ctx->tail) / sizeof(ctx->tail[0]); ++i)
            ctx->tail[i] ^=
                (limb_t)(0x3636363636363636 ^ 0x5c5c5c5c5c5c5c5c);

        sha256_init_h(ctx->h_opad);
        blst_sha256_block_data_order_v(ctx->h_opad, tail, 1);

        hmac_zero_tail(ctx);
        tail[32] = 0x80;
        tail[62] = 3;        /* (64+32)*8 in big endian */
        tail[63] = 0;
    }

    void hmac_update(HmacSha256Context* ctx, const unsigned char* inp,
        size_t len) noexcept
    {
        BLST_ASSERT(ctx != nullptr);
        BLST_ASSERT(inp != nullptr || len == 0);
        sha256_update(&ctx->ctx, inp, len);
    }

    void hmac_final(unsigned char md[32], HmacSha256Context* ctx) noexcept
    {
        BLST_ASSERT(md != nullptr);
        BLST_ASSERT(ctx != nullptr);

        unsigned char* tail = hmac_tail_bytes(ctx);
        sha256_final(tail, &ctx->ctx);
        blst_sha256_hcopy_v(ctx->ctx.h, ctx->h_opad);
        blst_sha256_block_data_order_v(ctx->ctx.h, tail, 1);
        blst_sha256_emit_v(md, ctx->ctx.h);
    }

    void hkdf_extract(unsigned char PRK[32],
        const void* salt, size_t salt_len,
        const void* IKM, size_t IKM_len,
#ifndef BLST_HKDF_TESTMODE
        int IKM_fixup,
#endif
        HmacSha256Context* ctx) noexcept
    {
        BLST_ASSERT(PRK != nullptr);
        BLST_ASSERT(ctx != nullptr);
        BLST_ASSERT(salt != nullptr || salt_len == 0);
        BLST_ASSERT(IKM != nullptr || IKM_len == 0);

        unsigned char zero[1] = { 0 };

        hmac_init(ctx, salt != NULL ? salt : zero, salt_len);
        hmac_update(ctx, (const unsigned char*)IKM, IKM_len);
#ifndef BLST_HKDF_TESTMODE
        if (IKM_fixup) {
            /* Section 2.3 KeyGen in BLS-signature draft */
            hmac_update(ctx, zero, 1);
        }
#endif
        hmac_final(PRK, ctx);
    }

    bool hkdf_expand(unsigned char* OKM, size_t L,
        const unsigned char PRK[32],
        const void* info, size_t info_len,
#ifndef BLST_HKDF_TESTMODE
        int info_fixup,
#endif
        HmacSha256Context* ctx) noexcept
    {
        static constexpr size_t digest_size = 32;
        static constexpr size_t max_output_size = 255 * digest_size;

        BLST_ASSERT(PRK != nullptr || L == 0);
        BLST_ASSERT(ctx != nullptr);

        if (L > max_output_size || (L != 0 && OKM == NULL) ||
            (info_len != 0 && info == NULL)) {
            return false;
        }
        if (L == 0)
            return true;

        unsigned char suffix[3];
        size_t suffix_len = 0;
#ifndef BLST_HKDF_TESTMODE
        if (info_fixup) {
            /* Section 2.3 KeyGen in BLS-signature draft */
            suffix[suffix_len++] = (unsigned char)(L >> 8);
            suffix[suffix_len++] = (unsigned char)(L);
        }
#endif
        const size_t counter_offset = suffix_len++;

        for (unsigned int counter = 1; L != 0; ++counter) {
            if (counter == 1) {
                hmac_init(ctx, PRK, digest_size);
            }
            else {
                hmac_init(ctx, NULL, 0);
                hmac_update(ctx, hmac_tail_bytes(ctx), digest_size);
            }

            if (info_len != 0)
                hmac_update(ctx, (const unsigned char*)info, info_len);
            suffix[counter_offset] = (unsigned char)counter;
            hmac_update(ctx, suffix, suffix_len);
            hmac_final(hmac_tail_bytes(ctx), ctx);

            const size_t chunk = L < digest_size ? L : digest_size;
            blst_sha256_bcopy_v(OKM, hmac_tail_bytes(ctx), chunk);
            OKM += chunk;
            L -= chunk;
        }

        return true;
    }

    } // namespace

#ifndef BLST_HKDF_TESTMODE
    namespace
    {
    void derive_key(pow256 SK, const void* IKM, size_t IKM_len,
        const void* salt, size_t salt_len,
        const void* info, size_t info_len,
        int version) noexcept
    {
        BLST_ASSERT(SK != nullptr);
        BLST_ASSERT(IKM != nullptr || IKM_len < 32);
        BLST_ASSERT(version == 3 || version == 4 || version == 5);

        struct {
            HmacSha256Context ctx;
            unsigned char PRK[32], OKM[48];
            UFr ufr;
            Fr key;
        } scratch;
        unsigned char salt_prime[32] = "BLS-SIG-KEYGEN-SALT-";

        if (IKM_len < 32 || (version > 4 && salt == NULL)) {
            vec_zero(SK, sizeof(pow256));
            return;
        }

        /* Vet |info| before using a caller-supplied length. */
        info_len = info == NULL ? 0 : info_len;

        if (salt == NULL) {
            salt = salt_prime;
            salt_len = 20;
        }

        if (version == 4) {
            /* salt = H(salt) */
            sha256_init(&scratch.ctx.ctx);
            sha256_update(&scratch.ctx.ctx, salt, salt_len);
            sha256_final(salt_prime, &scratch.ctx.ctx);
            salt = salt_prime;
            salt_len = sizeof(salt_prime);
        }

        while (1) {
            /* PRK = HKDF-Extract(salt, IKM || I2OSP(0, 1)) */
            hkdf_extract(scratch.PRK, salt, salt_len,
                IKM, IKM_len, 1, &scratch.ctx);

            /* OKM = HKDF-Expand(PRK, key_info || I2OSP(L, 2), L) */
            if (!hkdf_expand(scratch.OKM, sizeof(scratch.OKM), scratch.PRK,
                info, info_len, 1, &scratch.ctx)) {
                vec_zero(SK, sizeof(pow256));
                vec_zero(&scratch, sizeof(scratch));
                return;
            }

            /* SK = OS2IP(OKM) mod r */
            static_assert(sizeof(scratch.ufr.data) / sizeof(scratch.ufr.data[0]) == 8);
            scratch.ufr.data[6] = 0;
            scratch.ufr.data[7] = 0;
            limbs_from_be_bytes(scratch.ufr, scratch.OKM, sizeof(scratch.OKM));
            Fr::redc(scratch.key, scratch.ufr);
            /*
             * Given that mul_mont_sparse_256 has special boundary conditions
             * it's appropriate to mention that redc_mont_256 output is fully
             * reduced at this point. Because we started with 384-bit input,
             * one with most significant half smaller than the modulus.
             */
            Fr::mul(scratch.key, scratch.key, get_BLST_rRR());

            if (version < 4 || !vec_is_zero(scratch.key.l, sizeof(Fr)))
                break;

            /* salt = H(salt) */
            sha256_init(&scratch.ctx.ctx);
            sha256_update(&scratch.ctx.ctx, salt, salt_len);
            sha256_final(salt_prime, &scratch.ctx.ctx);
            salt = salt_prime;
            salt_len = sizeof(salt_prime);
        }

        le_bytes_from_limbs(SK, scratch.key.l, sizeof(pow256));

        /*
         * scrub the stack just in case next callee inadvertently flashes
         * a fragment across application boundary...
         */
        vec_zero(&scratch, sizeof(scratch));
    }

    bool parent_sk_to_lamport_pk(pow256 PK, const pow256 parent_SK,
        unsigned int index) noexcept;

    } // namespace

    void blst_keygen(pow256 SK, const void* IKM, size_t IKM_len,
        const void* info, size_t info_len) noexcept
    {
        derive_key(SK, IKM, IKM_len, NULL, 0, info, info_len, 4);
    }

    void blst_keygen_v3(pow256 SK, const void* IKM, size_t IKM_len,
        const void* info, size_t info_len) noexcept
    {
        derive_key(SK, IKM, IKM_len, NULL, 0, info, info_len, 3);
    }

    void blst_keygen_v4_5(pow256 SK, const void* IKM, size_t IKM_len,
        const void* salt, size_t salt_len,
        const void* info, size_t info_len) noexcept
    {
        derive_key(SK, IKM, IKM_len, salt, salt_len, info, info_len, 4);
    }

    void blst_keygen_v5(pow256 SK, const void* IKM, size_t IKM_len,
        const void* salt, size_t salt_len,
        const void* info, size_t info_len) noexcept
    {
        derive_key(SK, IKM, IKM_len, salt, salt_len, info, info_len, 5);
    }

    /*
     * https://eips.ethereum.org/EIPS/eip-2333
     */
    void blst_derive_master_eip2333(pow256 SK, const void* seed,
        size_t seed_len) noexcept
    {
        derive_key(SK, seed, seed_len, NULL, 0, NULL, 0, 4);
    }

    namespace
    {
    bool parent_sk_to_lamport_pk(pow256 PK, const pow256 parent_SK,
        unsigned int index) noexcept
    {
        BLST_ASSERT(PK != nullptr);
        BLST_ASSERT(parent_SK != nullptr);

        size_t i;
        struct {
            HmacSha256Context ctx;
            Sha256Context ret;
            unsigned char PRK[32], IKM[32];
            unsigned char lamport[255][32];
        } scratch;

        /* salt = I2OSP(index, 4) */
        unsigned char salt[4] = { (unsigned char)(index >> 24),
                                  (unsigned char)(index >> 16),
                                  (unsigned char)(index >> 8),
                                  (unsigned char)(index) };

        /* IKM = I2OSP(parent_SK, 32) */
        for (i = 0; i < 32; i++)
            scratch.IKM[i] = parent_SK[31 - i];

        /* lamport_0 = IKM_to_lamport_SK(IKM, salt) */
        hkdf_extract(scratch.PRK, salt, sizeof(salt), scratch.IKM, 32, 0,
            &scratch.ctx);
        if (!hkdf_expand(scratch.lamport[0], sizeof(scratch.lamport),
            scratch.PRK, NULL, 0, 0, &scratch.ctx)) {
            vec_zero(PK, sizeof(pow256));
            vec_zero(&scratch, sizeof(scratch));
            return false;
        }

        sha256_zero_buffer(scratch.ctx.ctx.buffer);
        scratch.ctx.ctx.buffer[32] = 0x80;
        scratch.ctx.ctx.buffer[62] = 1;    /* 32*8 in big endian */
        scratch.ctx.ctx.buffer[63] = 0;
        for (i = 0; i < 255; i++) {
            /* lamport_PK = lamport_PK | SHA256(lamport_0[i]) */
            sha256_init_h(scratch.ctx.ctx.h);
            blst_sha256_bcopy_v(scratch.ctx.ctx.buffer, scratch.lamport[i], 32);
            blst_sha256_block_data_order_v(scratch.ctx.ctx.h, scratch.ctx.ctx.buffer, 1);
            blst_sha256_emit_v(scratch.lamport[i], scratch.ctx.ctx.h);
        }

        /* compressed_lamport_PK = SHA256(lamport_PK) */
        sha256_init(&scratch.ret);
        sha256_update(&scratch.ret, scratch.lamport, sizeof(scratch.lamport));

        /* not_IKM = flip_bits(IKM) */
        for (i = 0; i < 32; i++)
            scratch.IKM[i] = ~scratch.IKM[i];

        /* lamport_1 = IKM_to_lamport_SK(not_IKM, salt) */
        hkdf_extract(scratch.PRK, salt, sizeof(salt), scratch.IKM, 32, 0,
            &scratch.ctx);
        if (!hkdf_expand(scratch.lamport[0], sizeof(scratch.lamport),
            scratch.PRK, NULL, 0, 0, &scratch.ctx)) {
            vec_zero(PK, sizeof(pow256));
            vec_zero(&scratch, sizeof(scratch));
            return false;
        }

        sha256_zero_buffer(scratch.ctx.ctx.buffer);
        scratch.ctx.ctx.buffer[32] = 0x80;
        scratch.ctx.ctx.buffer[62] = 1;
        for (i = 0; i < 255; i++) {
            /* lamport_PK = lamport_PK | SHA256(lamport_1[i]) */
            sha256_init_h(scratch.ctx.ctx.h);
            blst_sha256_bcopy_v(scratch.ctx.ctx.buffer, scratch.lamport[i], 32);
            blst_sha256_block_data_order_v(scratch.ctx.ctx.h, scratch.ctx.ctx.buffer, 1);
            blst_sha256_emit_v(scratch.lamport[i], scratch.ctx.ctx.h);
        }

        /* compressed_lamport_PK = SHA256(lamport_PK) */
        sha256_update(&scratch.ret, scratch.lamport, sizeof(scratch.lamport));
        sha256_final(PK, &scratch.ret);

        /*
         * scrub the stack just in case next callee inadvertently flashes
         * a fragment across application boundary...
         */
        vec_zero(&scratch, sizeof(scratch));
        return true;
    }

    } // namespace

    void blst_derive_child_eip2333(pow256 SK, const pow256 parent_SK,
        unsigned int child_index) noexcept
    {
        if (!parent_sk_to_lamport_pk(SK, parent_SK, child_index)) {
            vec_zero(SK, sizeof(pow256));
            return;
        }
        derive_key(SK, SK, sizeof(pow256), NULL, 0, NULL, 0, 4);
    }

    //BLST_EXTERN_C_END
#endif

}

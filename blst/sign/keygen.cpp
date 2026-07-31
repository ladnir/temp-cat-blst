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

#include <limits>

//BLST_EXTERN_C_BEGIN

namespace blst
{

    typedef struct {
        SHA256_CTX ctx;
        unsigned int h_ipad[8];
        unsigned int h_opad[8];
        union { limb_t l[64 / sizeof(limb_t)]; unsigned char c[64]; } tail;
    } HMAC_SHA256_CTX;

    void HMAC_init(HMAC_SHA256_CTX* ctx, const void* K, size_t K_len)
    {
        size_t i;

        if (K == NULL) {            /* reuse h_ipad and h_opad */
            blst_sha256_hcopy_v(ctx->ctx.h, ctx->h_ipad);
            ctx->ctx.N = 64;
            vec_zero(ctx->ctx.buf, sizeof(ctx->ctx.buf));
            ctx->ctx.off = 0;

            return;
        }

        vec_zero(ctx->tail.c, sizeof(ctx->tail));
        if (K_len > 64) {
            sha256_init(&ctx->ctx);
            sha256_update(&ctx->ctx, K, K_len);
            sha256_final(ctx->tail.c, &ctx->ctx);
        }
        else if (K_len != 0) {
            blst_sha256_bcopy_v(ctx->tail.c, K, K_len);
        }

        for (i = 0; i < 64 / sizeof(limb_t); i++)
            ctx->tail.l[i] ^= (limb_t)0x3636363636363636;

        sha256_init(&ctx->ctx);
        sha256_update(&ctx->ctx, ctx->tail.c, 64);
        blst_sha256_hcopy_v(ctx->h_ipad, ctx->ctx.h);

        for (i = 0; i < 64 / sizeof(limb_t); i++)
            ctx->tail.l[i] ^= (limb_t)(0x3636363636363636 ^ 0x5c5c5c5c5c5c5c5c);

        sha256_init_h(ctx->h_opad);
        blst_sha256_block_data_order_v(ctx->h_opad, ctx->tail.c, 1);

        vec_zero(ctx->tail.c, sizeof(ctx->tail));
        ctx->tail.c[32] = 0x80;
        ctx->tail.c[62] = 3;        /* (64+32)*8 in big endian */
        ctx->tail.c[63] = 0;
    }

    void HMAC_update(HMAC_SHA256_CTX* ctx, const unsigned char* inp,
        size_t len)
    {
        sha256_update(&ctx->ctx, inp, len);
    }

    void HMAC_final(unsigned char md[32], HMAC_SHA256_CTX* ctx)
    {
        sha256_final(ctx->tail.c, &ctx->ctx);
        blst_sha256_hcopy_v(ctx->ctx.h, ctx->h_opad);
        blst_sha256_block_data_order_v(ctx->ctx.h, ctx->tail.c, 1);
        blst_sha256_emit_v(md, ctx->ctx.h);
    }

    void HKDF_Extract(unsigned char PRK[32],
        const void* salt, size_t salt_len,
        const void* IKM, size_t IKM_len,
#ifndef __BLST_HKDF_TESTMODE__
        int IKM_fixup,
#endif
        HMAC_SHA256_CTX* ctx)
    {
        unsigned char zero[1] = { 0 };

        HMAC_init(ctx, salt != NULL ? salt : zero, salt_len);
        HMAC_update(ctx, (const unsigned char*)IKM, IKM_len);
#ifndef __BLST_HKDF_TESTMODE__
        if (IKM_fixup) {
            /* Section 2.3 KeyGen in BLS-signature draft */
            HMAC_update(ctx, zero, 1);
        }
#endif
        HMAC_final(PRK, ctx);
    }

    bool HKDF_Expand(unsigned char* OKM, size_t L,
        const unsigned char PRK[32],
        const void* info, size_t info_len,
#ifndef __BLST_HKDF_TESTMODE__
        int info_fixup,
#endif
        HMAC_SHA256_CTX* ctx)
    {
        if (info_len > std::numeric_limits<size_t>::max() - 3) {
            return false;
        }

        unsigned char* info_prime = (unsigned char*)malloc(info_len + 3);
        if (info_prime == NULL) {
            return false;
        }

        HMAC_init(ctx, PRK, 32);

        if (info_len != 0)
            blst_sha256_bcopy_v(info_prime, info, info_len);
#ifndef __BLST_HKDF_TESTMODE__
        if (info_fixup) {
            /* Section 2.3 KeyGen in BLS-signature draft */
            info_prime[info_len + 0] = (unsigned char)(L >> 8);
            info_prime[info_len + 1] = (unsigned char)(L);
            info_len += 2;
        }
#endif
        info_prime[info_len] = 1;   /* counter */
        HMAC_update(ctx, info_prime, info_len + 1);
        HMAC_final(ctx->tail.c, ctx);
        while (L > 32) {
            blst_sha256_hcopy_v((unsigned int*)OKM, (const unsigned int*)ctx->tail.c);
            OKM += 32; L -= 32;
            ++info_prime[info_len]; /* counter */
            HMAC_init(ctx, NULL, 0);
            HMAC_update(ctx, ctx->tail.c, 32);
            HMAC_update(ctx, info_prime, info_len + 1);
            HMAC_final(ctx->tail.c, ctx);
        }
        blst_sha256_bcopy_v(OKM, ctx->tail.c, L);

        free(info_prime);
        return true;
    }

#ifndef __BLST_HKDF_TESTMODE__
    void keygen(pow256 SK, const void* IKM, size_t IKM_len,
        const void* salt, size_t salt_len,
        const void* info, size_t info_len,
        int version)
    {
        struct {
            HMAC_SHA256_CTX ctx;
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
            HKDF_Extract(scratch.PRK, salt, salt_len,
                IKM, IKM_len, 1, &scratch.ctx);

            /* OKM = HKDF-Expand(PRK, key_info || I2OSP(L, 2), L) */
            if (!HKDF_Expand(scratch.OKM, sizeof(scratch.OKM), scratch.PRK,
                info, info_len, 1, &scratch.ctx)) {
                vec_zero(SK, sizeof(pow256));
                vec_zero(&scratch, sizeof(scratch));
                return;
            }

            /* SK = OS2IP(OKM) mod r */
            vec_zero(scratch.ufr, sizeof(scratch.ufr));
            limbs_from_be_bytes(scratch.ufr, scratch.OKM, sizeof(scratch.OKM));
            Fr::redc(scratch.key, scratch.ufr);
            /*
             * Given that mul_mont_sparse_256 has special boundary conditions
             * it's appropriate to mention that redc_mont_256 output is fully
             * reduced at this point. Because we started with 384-bit input,
             * one with most significant half smaller than the modulus.
             */
            Fr::mul(scratch.key, scratch.key, get_BLST_rRR());

            if (version < 4 || !vec_is_zero(scratch.key.data, sizeof(Fr)))
                break;

            /* salt = H(salt) */
            sha256_init(&scratch.ctx.ctx);
            sha256_update(&scratch.ctx.ctx, salt, salt_len);
            sha256_final(salt_prime, &scratch.ctx.ctx);
            salt = salt_prime;
            salt_len = sizeof(salt_prime);
        }

        le_bytes_from_limbs(SK, scratch.key.data, sizeof(pow256));

        /*
         * scrub the stack just in case next callee inadvertently flashes
         * a fragment across application boundary...
         */
        vec_zero(&scratch, sizeof(scratch));
    }

    void blst_keygen(pow256 SK, const void* IKM, size_t IKM_len,
        const void* info, size_t info_len)
    {
        keygen(SK, IKM, IKM_len, NULL, 0, info, info_len, 4);
    }

    void blst_keygen_v3(pow256 SK, const void* IKM, size_t IKM_len,
        const void* info, size_t info_len)
    {
        keygen(SK, IKM, IKM_len, NULL, 0, info, info_len, 3);
    }

    void blst_keygen_v4_5(pow256 SK, const void* IKM, size_t IKM_len,
        const void* salt, size_t salt_len,
        const void* info, size_t info_len)
    {
        keygen(SK, IKM, IKM_len, salt, salt_len, info, info_len, 4);
    }

    void blst_keygen_v5(pow256 SK, const void* IKM, size_t IKM_len,
        const void* salt, size_t salt_len,
        const void* info, size_t info_len)
    {
        keygen(SK, IKM, IKM_len, salt, salt_len, info, info_len, 5);
    }

    /*
     * https://eips.ethereum.org/EIPS/eip-2333
     */
    void blst_derive_master_eip2333(pow256 SK, const void* seed, size_t seed_len)
    {
        keygen(SK, seed, seed_len, NULL, 0, NULL, 0, 4);
    }

    void parent_SK_to_lamport_PK(pow256 PK, const pow256 parent_SK,
        unsigned int index)
    {
        size_t i;
        struct {
            HMAC_SHA256_CTX ctx;
            SHA256_CTX ret;
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
        HKDF_Extract(scratch.PRK, salt, sizeof(salt), scratch.IKM, 32, 0,
            &scratch.ctx);
        if (!HKDF_Expand(scratch.lamport[0], sizeof(scratch.lamport),
            scratch.PRK, NULL, 0, 0, &scratch.ctx)) {
            vec_zero(PK, sizeof(pow256));
            vec_zero(&scratch, sizeof(scratch));
            return;
        }

        vec_zero(scratch.ctx.ctx.buf, sizeof(scratch.ctx.ctx.buf));
        scratch.ctx.ctx.buf[32] = 0x80;
        scratch.ctx.ctx.buf[62] = 1;    /* 32*8 in big endian */
        scratch.ctx.ctx.buf[63] = 0;
        for (i = 0; i < 255; i++) {
            /* lamport_PK = lamport_PK | SHA256(lamport_0[i]) */
            sha256_init_h(scratch.ctx.ctx.h);
            blst_sha256_bcopy_v(scratch.ctx.ctx.buf, scratch.lamport[i], 32);
            blst_sha256_block_data_order_v(scratch.ctx.ctx.h, scratch.ctx.ctx.buf, 1);
            blst_sha256_emit_v(scratch.lamport[i], scratch.ctx.ctx.h);
        }

        /* compressed_lamport_PK = SHA256(lamport_PK) */
        sha256_init(&scratch.ret);
        sha256_update(&scratch.ret, scratch.lamport, sizeof(scratch.lamport));

        /* not_IKM = flip_bits(IKM) */
        for (i = 0; i < 32; i++)
            scratch.IKM[i] = ~scratch.IKM[i];

        /* lamport_1 = IKM_to_lamport_SK(not_IKM, salt) */
        HKDF_Extract(scratch.PRK, salt, sizeof(salt), scratch.IKM, 32, 0,
            &scratch.ctx);
        if (!HKDF_Expand(scratch.lamport[0], sizeof(scratch.lamport),
            scratch.PRK, NULL, 0, 0, &scratch.ctx)) {
            vec_zero(PK, sizeof(pow256));
            vec_zero(&scratch, sizeof(scratch));
            return;
        }

        vec_zero(scratch.ctx.ctx.buf, sizeof(scratch.ctx.ctx.buf));
        scratch.ctx.ctx.buf[32] = 0x80;
        scratch.ctx.ctx.buf[62] = 1;
        for (i = 0; i < 255; i++) {
            /* lamport_PK = lamport_PK | SHA256(lamport_1[i]) */
            sha256_init_h(scratch.ctx.ctx.h);
            blst_sha256_bcopy_v(scratch.ctx.ctx.buf, scratch.lamport[i], 32);
            blst_sha256_block_data_order_v(scratch.ctx.ctx.h, scratch.ctx.ctx.buf, 1);
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
    }

    void blst_derive_child_eip2333(pow256 SK, const pow256 parent_SK,
        unsigned int child_index)
    {
        parent_SK_to_lamport_PK(SK, parent_SK, child_index);
        keygen(SK, SK, sizeof(pow256), NULL, 0, NULL, 0, 4);
    }

    //BLST_EXTERN_C_END
}

#endif

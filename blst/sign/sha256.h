/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_SHA256_H__
#define __BLS12_381_ASM_SHA256_H__

#include "blst/common/vect.h"
#include "blst/config.h"
#include "blst/backend/backend.h"


namespace blst
{

    /*
     * If SHA256_CTX conflicts with something, just redefine it to alternative
     * custom name prior including this header.
     */
    typedef struct {
        unsigned int h[8];
        unsigned long long N;
        unsigned char buf[64];
        size_t off;
    } SHA256_CTX;


    BLST_HD BLST_INLINE void sha256_init_h(unsigned int h[8])
    {
        h[0] = 0x6a09e667U;
        h[1] = 0xbb67ae85U;
        h[2] = 0x3c6ef372U;
        h[3] = 0xa54ff53aU;
        h[4] = 0x510e527fU;
        h[5] = 0x9b05688cU;
        h[6] = 0x1f83d9abU;
        h[7] = 0x5be0cd19U;
    }

    BLST_HD BLST_INLINE void sha256_init(SHA256_CTX* ctx)
    {
        sha256_init_h(ctx->h);
        ctx->N = 0;
        vec_zero(ctx->buf, sizeof(ctx->buf));
        ctx->off = 0;
    }

    BLST_HD BLST_INLINE void sha256_update(SHA256_CTX* ctx, const void* _inp, size_t len)
    {
        size_t n;
        const unsigned char* inp = (const unsigned char*)_inp;

        ctx->N += len;

        if ((len != 0) & ((n = ctx->off) != 0)) {
            size_t rem = sizeof(ctx->buf) - n;

            if (rem > len) {
                blst_sha256_bcopy_v(ctx->buf + n, inp, len);
                ctx->off += len;
                return;
            }
            else {
                blst_sha256_bcopy_v(ctx->buf + n, inp, rem);
                inp += rem;
                len -= rem;
                blst_sha256_block_data_order_v(ctx->h, ctx->buf, 1);
                vec_zero(ctx->buf, sizeof(ctx->buf));
                ctx->off = 0;
            }
        }

        n = len / sizeof(ctx->buf);
        if (n > 0) {
            blst_sha256_block_data_order_v(ctx->h, inp, n);
            n *= sizeof(ctx->buf);
            inp += n;
            len -= n;
        }

        if (len)
            blst_sha256_bcopy_v(ctx->buf, inp, ctx->off = len);
    }

#define __TOBE32(ptr, val) ((ptr)[0] = (unsigned char)((val)>>24), \
                            (ptr)[1] = (unsigned char)((val)>>16), \
                            (ptr)[2] = (unsigned char)((val)>>8),  \
                            (ptr)[3] = (unsigned char)(val))

    BLST_HD BLST_INLINE void sha256_final(unsigned char md[32], SHA256_CTX* ctx)
    {
        unsigned long long bits = ctx->N * 8;
        size_t n = ctx->off;
        unsigned char* tail;

        ctx->buf[n++] = 0x80;

        if (n > (sizeof(ctx->buf) - 8)) {
            blst_sha256_block_data_order_v(ctx->h, ctx->buf, 1);
            vec_zero(ctx->buf, sizeof(ctx->buf));
        }

        tail = ctx->buf + sizeof(ctx->buf) - 8;
        __TOBE32(tail, (unsigned int)(bits >> 32));
        __TOBE32(tail + 4, (unsigned int)bits);
        blst_sha256_block_data_order_v(ctx->h, ctx->buf, 1);
        blst_sha256_emit_v(md, ctx->h);
    }

#undef __TOBE32
}

#endif

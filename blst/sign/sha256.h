/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_SIGN_SHA256_H
#define BLST_SIGN_SHA256_H

#include <cstdint>
#include <limits>

#include "blst/common/vect.h"
#include "blst/config.h"
#include "blst/backend/backend.h"

namespace blst
{
    static_assert(std::numeric_limits<unsigned int>::digits == 32,
        "SHA-256 requires a 32-bit unsigned int");
    static_assert(std::numeric_limits<std::uint64_t>::digits == 64,
        "SHA-256 requires a 64-bit byte counter");
    static_assert(sizeof(size_t) == sizeof(std::uint64_t),
        "CatBLST supports only 64-bit size_t");

    /**
     * Native streaming SHA-256 state.
     *
     * Initialize before update or final. Final consumes the state; reinitialize
     * before reuse. Owners of secret-bearing state call `sha256_clear` when its
     * lifetime ends. Callers must not modify fields directly.
     */
    struct Sha256Context {
        unsigned int h[8];
        std::uint64_t total_bytes;
        unsigned char buffer[64];
        size_t buffered_bytes;
    };

    /* Frozen BTx/BTE still consumes this internal legacy spelling. */
    using SHA256_CTX = Sha256Context;

    constexpr std::uint64_t sha256_max_message_bytes =
        (std::numeric_limits<std::uint64_t>::max)() / 8;

    BLST_HD BLST_FORCEINLINE void sha256_zero_buffer(
        unsigned char buffer[64]) noexcept
    {
        /* Functional initialization: permit vector stores and dead-store folding. */
        for (size_t i = 0; i < 64; ++i)
            buffer[i] = 0;
    }

    BLST_HD BLST_FORCEINLINE void sha256_store_be32(unsigned char out[4],
        unsigned int value) noexcept
    {
        out[0] = (unsigned char)(value >> 24);
        out[1] = (unsigned char)(value >> 16);
        out[2] = (unsigned char)(value >> 8);
        out[3] = (unsigned char)value;
    }

    BLST_HD BLST_INLINE void sha256_init_h(unsigned int h[8]) noexcept
    {
        BLST_ASSERT(h != nullptr);

        h[0] = 0x6a09e667U;
        h[1] = 0xbb67ae85U;
        h[2] = 0x3c6ef372U;
        h[3] = 0xa54ff53aU;
        h[4] = 0x510e527fU;
        h[5] = 0x9b05688cU;
        h[6] = 0x1f83d9abU;
        h[7] = 0x5be0cd19U;
    }

    /** Initializes or reinitializes a context. No allocation or exception. */
    BLST_HD BLST_INLINE void sha256_init(Sha256Context* ctx) noexcept
    {
        BLST_ASSERT(ctx != nullptr);

        sha256_init_h(ctx->h);
        ctx->total_bytes = 0;
        sha256_zero_buffer(ctx->buffer);
        ctx->buffered_bytes = 0;
    }

    /**
     * Appends exactly len bytes. inp may be null only when len is zero.
     * A context must be initialized and not finalized. The SHA-256 message
     * domain is fewer than 2^64 bits; Debug builds enforce that bound.
     */
    BLST_HD BLST_INLINE void sha256_update(Sha256Context* ctx,
        const void* input, size_t len) noexcept
    {
        BLST_ASSERT(ctx != nullptr);
        BLST_ASSERT(input != nullptr || len == 0);
        BLST_ASSERT(ctx->buffered_bytes < sizeof(ctx->buffer));
        BLST_ASSERT(ctx->buffered_bytes ==
            (size_t)(ctx->total_bytes & (sizeof(ctx->buffer) - 1)));
        BLST_ASSERT(ctx->total_bytes <= sha256_max_message_bytes);
        BLST_ASSERT(len <= sha256_max_message_bytes - ctx->total_bytes);

        size_t n;
        const unsigned char* inp = (const unsigned char*)input;

        ctx->total_bytes += len;

        if ((len != 0) & ((n = ctx->buffered_bytes) != 0)) {
            size_t rem = sizeof(ctx->buffer) - n;

            if (rem > len) {
                blst_sha256_bcopy_v(ctx->buffer + n, inp, len);
                ctx->buffered_bytes += len;
                return;
            }
            else {
                blst_sha256_bcopy_v(ctx->buffer + n, inp, rem);
                inp += rem;
                len -= rem;
                blst_sha256_block_data_order_v(ctx->h, ctx->buffer, 1);
                sha256_zero_buffer(ctx->buffer);
                ctx->buffered_bytes = 0;
            }
        }

        n = len / sizeof(ctx->buffer);
        if (n > 0) {
            blst_sha256_block_data_order_v(ctx->h, inp, n);
            n *= sizeof(ctx->buffer);
            inp += n;
            len -= n;
        }

        if (len)
            blst_sha256_bcopy_v(ctx->buffer, inp,
                ctx->buffered_bytes = len);
    }

    /**
     * Publishes the 32-byte digest and consumes the current context state.
     * md must be writable. Reinitialize before any subsequent update/final.
     * This internal operation does not erase caller-owned context storage.
     */
    BLST_HD BLST_INLINE void sha256_final(unsigned char md[32],
        Sha256Context* ctx) noexcept
    {
        BLST_ASSERT(md != nullptr);
        BLST_ASSERT(ctx != nullptr);
        BLST_ASSERT(ctx->buffered_bytes < sizeof(ctx->buffer));
        BLST_ASSERT(ctx->buffered_bytes ==
            (size_t)(ctx->total_bytes & (sizeof(ctx->buffer) - 1)));
        BLST_ASSERT(ctx->total_bytes <= sha256_max_message_bytes);

        std::uint64_t bits = ctx->total_bytes * 8;
        size_t n = ctx->buffered_bytes;
        unsigned char* tail;

        ctx->buffer[n++] = 0x80;

        if (n > (sizeof(ctx->buffer) - 8)) {
            blst_sha256_block_data_order_v(ctx->h, ctx->buffer, 1);
            sha256_zero_buffer(ctx->buffer);
        }

        tail = ctx->buffer + sizeof(ctx->buffer) - 8;
        sha256_store_be32(tail, (unsigned int)(bits >> 32));
        sha256_store_be32(tail + 4, (unsigned int)bits);
        blst_sha256_block_data_order_v(ctx->h, ctx->buffer, 1);
        blst_sha256_emit_v(md, ctx->h);
    }

    /** Explicitly erases caller-owned context storage. */
    BLST_HD BLST_INLINE void sha256_clear(Sha256Context* ctx) noexcept
    {
        BLST_ASSERT(ctx != nullptr);
        vec_zero(ctx, sizeof(*ctx));
    }
}

#endif

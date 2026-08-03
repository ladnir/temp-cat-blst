/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "bindings/blst_scalar.h"

/* The public header must leave generic consumer typedef names available. */
typedef unsigned char byte;
typedef unsigned long long limb_t;
#include <stddef.h>

#ifdef __cplusplus
# define AUDIT_STATIC_ASSERT(condition) static_assert(condition)
# define AUDIT_ALIGNOF(type) alignof(type)
#else
# define AUDIT_STATIC_ASSERT(condition) _Static_assert(condition, #condition)
# define AUDIT_ALIGNOF(type) _Alignof(type)
#endif

AUDIT_STATIC_ASSERT(sizeof(blst_fr) == 32);
AUDIT_STATIC_ASSERT(sizeof(blst_fp) == 48);
AUDIT_STATIC_ASSERT(sizeof(blst_fp2) == 96);
AUDIT_STATIC_ASSERT(sizeof(blst_fp6) == 288);
AUDIT_STATIC_ASSERT(sizeof(blst_fp12) == 576);
AUDIT_STATIC_ASSERT(sizeof(blst_p1) == 144);
AUDIT_STATIC_ASSERT(sizeof(blst_p1_affine) == 96);
AUDIT_STATIC_ASSERT(sizeof(blst_p2) == 288);
AUDIT_STATIC_ASSERT(sizeof(blst_p2_affine) == 192);
AUDIT_STATIC_ASSERT(AUDIT_ALIGNOF(blst_fr) == AUDIT_ALIGNOF(blst_limb_t));
AUDIT_STATIC_ASSERT(AUDIT_ALIGNOF(blst_fp12) == AUDIT_ALIGNOF(blst_limb_t));
AUDIT_STATIC_ASSERT(AUDIT_ALIGNOF(blst_p2) == AUDIT_ALIGNOF(blst_limb_t));
AUDIT_STATIC_ASSERT(offsetof(blst_fr, l) == 0);
AUDIT_STATIC_ASSERT(offsetof(blst_fp, l) == 0);
AUDIT_STATIC_ASSERT(offsetof(blst_fp2, fp) == 0);
AUDIT_STATIC_ASSERT(offsetof(blst_fp6, fp2) == 0);
AUDIT_STATIC_ASSERT(offsetof(blst_fp12, fp6) == 0);
AUDIT_STATIC_ASSERT(offsetof(blst_p1, x) == 0);
AUDIT_STATIC_ASSERT(offsetof(blst_p1_affine, x) == 0);
AUDIT_STATIC_ASSERT(offsetof(blst_p2, x) == 0);
AUDIT_STATIC_ASSERT(offsetof(blst_p2_affine, x) == 0);

static void (*const scalar_from_u32)(blst_scalar*, const uint32_t*) =
    blst_scalar_from_uint32;
static void (*const u32_from_scalar)(uint32_t*, const blst_scalar*) =
    blst_uint32_from_scalar;
static void (*const scalar_from_u64)(blst_scalar*, const uint64_t*) =
    blst_scalar_from_uint64;
static void (*const u64_from_scalar)(uint64_t*, const blst_scalar*) =
    blst_uint64_from_scalar;
static void (*const scalar_from_be)(blst_scalar*, const blst_byte*) =
    blst_scalar_from_bendian;
static void (*const be_from_scalar)(blst_byte*, const blst_scalar*) =
    blst_bendian_from_scalar;
static void (*const scalar_from_le)(blst_scalar*, const blst_byte*) =
    blst_scalar_from_lendian;
static void (*const le_from_scalar)(blst_byte*, const blst_scalar*) =
    blst_lendian_from_scalar;

int main(void)
{
    return scalar_from_u32 == 0 || u32_from_scalar == 0 ||
        scalar_from_u64 == 0 || u64_from_scalar == 0 ||
        scalar_from_be == 0 || be_from_scalar == 0 ||
        scalar_from_le == 0 || le_from_scalar == 0;
}

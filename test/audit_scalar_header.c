/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "bindings/blst_scalar.h"

static void (*const scalar_from_u32)(blst_scalar*, const uint32_t*) =
    blst_scalar_from_uint32;
static void (*const u32_from_scalar)(uint32_t*, const blst_scalar*) =
    blst_uint32_from_scalar;
static void (*const scalar_from_u64)(blst_scalar*, const uint64_t*) =
    blst_scalar_from_uint64;
static void (*const u64_from_scalar)(uint64_t*, const blst_scalar*) =
    blst_uint64_from_scalar;
static void (*const scalar_from_be)(blst_scalar*, const byte*) =
    blst_scalar_from_bendian;
static void (*const be_from_scalar)(byte*, const blst_scalar*) =
    blst_bendian_from_scalar;
static void (*const scalar_from_le)(blst_scalar*, const byte*) =
    blst_scalar_from_lendian;
static void (*const le_from_scalar)(byte*, const blst_scalar*) =
    blst_lendian_from_scalar;

int main(void)
{
    return scalar_from_u32 == 0 || u32_from_scalar == 0 ||
        scalar_from_u64 == 0 || u64_from_scalar == 0 ||
        scalar_from_be == 0 || be_from_scalar == 0 ||
        scalar_from_le == 0 || le_from_scalar == 0;
}

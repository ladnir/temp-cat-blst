/*
 * Copyright Category Labs and Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "bindings/blst_avx.h"

static bool fill_random(void* context, blst_byte* out, size_t out_len)
{
    size_t i;
    const blst_byte seed = context != NULL
        ? *(const blst_byte*)context
        : 0;
    for (i = 0; i < out_len; ++i)
        out[i] = (blst_byte)(seed + i);
    return true;
}

bool blst_avx_batch_membership_c_header_probe(
    const blst_avx_p1* const* p1_points,
    const blst_avx_p2* const* p2_points,
    size_t count,
    blst_byte* seed,
    void* p2_workspace,
    size_t p2_workspace_bytes)
{
    return blst_avx_p1s_in_g1(p1_points, count, fill_random, seed) &&
        blst_avx_p2s_in_g2(
            p2_points, count, fill_random, seed,
            p2_workspace, p2_workspace_bytes);
}

bool blst_avx_secret_g1_msm_c_header_probe(
    blst_avx_p1* out,
    const blst_avx_p1* points,
    const blst_fr* scalars,
    size_t count,
    void* scratch,
    size_t scratch_bytes)
{
    return blst_avx_p1s_mult_ct_scratch_sizeof() != 0 &&
        blst_avx_p1s_mult_ct(
            out, points, scalars, count, scratch, scratch_bytes);
}

bool blst_avx_public_g1_msm_c_header_probe(
    blst_avx_p1* out,
    const blst_avx_p1* points,
    const blst_fr* scalars,
    size_t count,
    void* scratch,
    size_t scratch_bytes)
{
    return blst_avx_p1s_mult_fr_vartime_scratch_sizeof(count) <=
            scratch_bytes &&
        blst_avx_p1s_mult_fr_vartime(
            out, points, scalars, count, scratch, scratch_bytes);
}

bool blst_avx_public_g2_mult_c_header_probe(
    blst_avx_p2* out,
    const blst_avx_p2* point,
    const blst_fr* scalar)
{
    return blst_avx_p2_mult_fr_vartime(out, point, scalar);
}

bool blst_avx_public_g2_msm_c_header_probe(
    blst_avx_p2* out,
    const blst_avx_p2* points,
    const blst_fr* scalars,
    size_t count,
    void* scratch,
    size_t scratch_bytes)
{
    return blst_avx_p2s_mult_fr_vartime_scratch_sizeof(count) <=
            scratch_bytes &&
        blst_avx_p2s_mult_fr_vartime(
            out, points, scalars, count, scratch, scratch_bytes);
}

bool blst_avx_miller_loop_n_c_header_probe(
    blst_avx_fp12* out,
    const blst_avx_p1_affine* g1,
    const blst_avx_p2_affine* g2,
    size_t count)
{
    return blst_avx_miller_loop_n(out, g1, g2, count);
}

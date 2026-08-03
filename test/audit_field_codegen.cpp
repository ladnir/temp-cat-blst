/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/fields/Fp.h"
#include "blst/fields/Fr.h"

#if defined(_MSC_VER)
# define AUDIT_NOINLINE __declspec(noinline)
#else
# define AUDIT_NOINLINE __attribute__((noinline))
#endif

extern "C" AUDIT_NOINLINE void audit_fp_sqr_n_mul_3(
    blst::Fp* out,
    const blst::Fp* in)
{
    blst::Fp tmp;
    blst::Fp::sqr_n_mul(tmp, *in, 3, blst::Fp::one());
    *out = tmp;
}

extern "C" AUDIT_NOINLINE void audit_fr_exp_zero(
    blst::Fr* out,
    const blst::Fr* in)
{
    blst::Fr::exp(*out, *in, nullptr, 0);
}

extern "C" AUDIT_NOINLINE void audit_fp_exp_zero(
    blst::Fp* out,
    const blst::Fp* in)
{
    blst::Fp::exp(*out, *in, nullptr, 0);
}

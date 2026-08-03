/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/common/defines.h"

#if defined(_MSC_VER)
# define AUDIT_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define AUDIT_NOINLINE __attribute__((noinline))
#else
# define AUDIT_NOINLINE
#endif

extern "C" {

AUDIT_NOINLINE blst::limb_t audit_launder_select(
    const blst::limb_t when_true,
    const blst::limb_t when_false,
    blst::bool_t select_true)
{
    BLST_launder(select_true);
    const blst::limb_t mask = 0 - select_true;
    return when_false ^ ((when_true ^ when_false) & mask);
}

AUDIT_NOINLINE void audit_launder_swap(
    blst::limb_t* const lhs,
    blst::limb_t* const rhs,
    blst::bool_t swap)
{
    BLST_launder(swap);
    const blst::limb_t mask = 0 - swap;
    const blst::limb_t delta = (*lhs ^ *rhs) & mask;
    *lhs ^= delta;
    *rhs ^= delta;
}

AUDIT_NOINLINE blst::limb_t audit_launder_conditional_zero(
    const blst::limb_t value,
    blst::bool_t zero)
{
    BLST_launder(zero);
    const blst::limb_t keep_mask = 0 - (zero ^ 1);
    return value & keep_mask;
}

AUDIT_NOINLINE void audit_launder_secure_zero(
    void* const output,
    const std::size_t size)
{
    volatile unsigned char* current =
        static_cast<volatile unsigned char*>(output);
    for (std::size_t i = 0; i < size; ++i) {
        current[i] = 0;
    }
    BLST_launder(current);
}

AUDIT_NOINLINE void audit_launder_portable_mask(
    blst::limb_t output[4],
    const blst::limb_t reduced[4],
    const blst::limb_t candidate[4],
    const blst::limb_t carry,
    const blst::limb_t borrow)
{
    blst::limb_t mask = carry - borrow;
    BLST_launder(mask);
    for (std::size_t i = 0; i < 4; ++i) {
        output[i] = (reduced[i] & ~mask) | (candidate[i] & mask);
    }
}

}

#undef AUDIT_NOINLINE

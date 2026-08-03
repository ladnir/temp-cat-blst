// Copyright 2026 Peter Rindal.
// SPDX-License-Identifier: Apache-2.0

#include "blst/backend/backend.h"

extern "C" void audit_lshift_v256_by_3(blst::limb_t* ret,
    const blst::limb_t* a, const blst::limb_t* p)
{
    blst::lshift_mod_v256(ret, a, 3, p);
}

extern "C" void audit_rshift_v256_by_3(blst::limb_t* ret,
    const blst::limb_t* a, const blst::limb_t* p)
{
    blst::rshift_mod_v256(ret, a, 3, p);
}

extern "C" void audit_lshift_v384_by_3(blst::limb_t* ret,
    const blst::limb_t* a, const blst::limb_t* p)
{
    blst::lshift_mod_v384(ret, a, 3, p);
}

extern "C" void audit_rshift_v384_by_3(blst::limb_t* ret,
    const blst::limb_t* a, const blst::limb_t* p)
{
    blst::rshift_mod_v384(ret, a, 3, p);
}

extern "C" void audit_sha256_selected(unsigned int* state, const void* input,
    size_t blocks)
{
    blst::blst_sha256_block_data_order_v(state, input, blocks);
}

extern "C" void audit_sha256_one_block(unsigned int* state, const void* input)
{
    blst::blst_sha256_block_data_order_v(state, input, 1);
}

extern "C" void audit_sha256_copy_32(void* dst, const void* src)
{
    blst::blst_sha256_bcopy_v(dst, src, 32);
}

extern "C" blst::bool_t audit_vec_is_zero_32(const void* input)
{
    return blst::vec_is_zero(input, 32);
}

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

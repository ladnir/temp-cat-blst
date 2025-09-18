

#pragma once
#include "blst/curve/point.h"
#include "blst/fields/Fp12.h"

//BLST_EXTERN_C_BEGIN
void blst_miller_loop(vec384fp12 ret, const POINTonE2_affine* Q,
    const POINTonE1_affine* P);
void blst_final_exp(vec384fp12 ret, const vec384fp12 f);

void blst_miller_loop_n(vec384fp12 out, const POINTonE2_affine* const Qs[],
    const POINTonE1_affine* const Ps[],
    size_t n);

int blst_fp12_in_group(const vec384fp12 f);

bool_t is_cyclotomic(const vec384fp12 f);
void blst_precompute_lines(vec384fp6 Qlines[68], const POINTonE2_affine* Q);



void blst_miller_loop_lines(vec384fp12 ret, const vec384fp6 Qlines[68],
    const POINTonE1_affine* P);
//BLST_EXTERN_C_END

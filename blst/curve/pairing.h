

#pragma once
#include "blst/curve/point.h"
#include "blst/fields/Fp12.h"

namespace blst
{

    //BLST_EXTERN_C_BEGIN
    BLST_HD void blst_miller_loop(Fp12& ret, const POINTonE2_affine* Q,
        const POINTonE1_affine* P);
    BLST_HD  void blst_final_exp(Fp12& ret, const Fp12& f);

    BLST_HD void blst_miller_loop_n(
        Fp12& out, 
        const POINTonE2_affine* const Qs[],
        const POINTonE1_affine* const Ps[],
        size_t n);

    BLST_HD int blst_fp12_in_group(const Fp12& f);

    BLST_HD bool_t is_cyclotomic(const Fp12& f);
    BLST_HD void blst_precompute_lines(Fp6 Qlines[68], const POINTonE2_affine& Q);



    BLST_HD void blst_miller_loop_lines(Fp12& ret, const Fp6 Qlines[68],
        const POINTonE1_affine& P);
    //BLST_EXTERN_C_END

}


#pragma once
#include "blst/curve/point.h"
#include "blst/fields/Fp12.h"

namespace blst
{

    //BLST_EXTERN_C_BEGIN
    BLST_HD void blst_miller_loop(Fp12& ret, const G2AffinePoint* Q,
        const G1AffinePoint* P);
    //BLST_HD  void blst_final_exp(Fp12& ret, const Fp12& f);

    BLST_HD void blst_miller_loop_n(
        Fp12& out, 
        const G2AffinePoint* const Qs[],
        const G1AffinePoint* const Ps[],
        size_t n);

    BLST_HD void blst_precompute_lines(Fp6 Qlines[68], const G2AffinePoint& Q);



    BLST_HD void blst_miller_loop_lines(Fp12& ret, const Fp6 Qlines[68],
        const G1AffinePoint& P);
    //BLST_EXTERN_C_END

}
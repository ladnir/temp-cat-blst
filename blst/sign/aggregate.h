#pragma once

#include "blst/curve/point.h"
#include "blst/curve/GT.h"
//BLST_EXTERN_C_BEGIN

#ifndef N_MAX
# define N_MAX 8
#endif
namespace blst
{


    typedef struct {
        G1Point e1;
        G2Point e2;
    } AggregatedSignature;
    typedef struct {
        unsigned int ctrl;
        unsigned int nelems;
        const void* DST;
        size_t DST_len;
        GT gt;
        AggregatedSignature AggrSign;
        G2AffinePoint Q[N_MAX];
        G1AffinePoint P[N_MAX];
    } PAIRING;

    BLST_HD int blst_fp12_finalverify(const GT& GT1, const GT& GT2);

}
//BLST_EXTERN_C_END

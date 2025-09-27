#include "blst/curve/point.h"
//BLST_EXTERN_C_BEGIN

#ifndef N_MAX
# define N_MAX 8
#endif

namespace blst
{


    typedef union { POINTonE1 e1; POINTonE2 e2; } AggregatedSignature;
    typedef struct {
        unsigned int ctrl;
        unsigned int nelems;
        const void* DST;
        size_t DST_len;
        Fp12 GT;
        AggregatedSignature AggrSign;
        POINTonE2_affine Q[N_MAX];
        POINTonE1_affine P[N_MAX];
    } PAIRING;

    BLST_HD int blst_fp12_finalverify(const Fp12& GT1, const Fp12& GT2);

}
//BLST_EXTERN_C_END

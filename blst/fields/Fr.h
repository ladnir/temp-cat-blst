#ifndef BLST_FR_H
#define BLST_FR_H

#include "blst/scalar/fields/Fr.h"

namespace blst
{
    using scalar::Fr;

    using scalar::BLST_GT_LAMBDA2_HOST;
    using scalar::BLST_GT_LAMBDA_HOST;
    using scalar::BLST_r_HOST;
    using scalar::BLST_r_one_HOST;
    using scalar::BLST_rRR_HOST;
    using scalar::BLST_r_zero_HOST;

    using scalar::ct_inverse_mod_v256;
    using scalar::get_BLST_gt_lambda;
    using scalar::get_BLST_gt_lambda2;
    using scalar::get_BLST_r;
    using scalar::get_BLST_rRR;

    using scalar::sk_add_n_check;
    using scalar::sk_check;
    using scalar::sk_inverse;
    using scalar::sk_mul_n_check;
    using scalar::sk_sub_n_check;
}

#endif

#pragma once

#include "blst/scalar/curve/multi_scalar.h"

namespace blst
{
    template<typename T>
    using span = scalar::span<T>;

    using scalar::WindowExtract;
    using scalar::hash_point;
    using scalar::msm_add_or_sub;
    using scalar::msm_mult;
    using scalar::msm_mult_pippenger;
    using scalar::msm_mult_pippenger_scratch_sizeof_;
    using scalar::msm_mult_table;
    using scalar::msm_precompute_wbits;
    using scalar::msm_to_jacobian;
    using scalar::pippenger_window_size;
    using scalar::point_gather_booth_wbits;
    using scalar::point_precompute_row_wbits;
    using scalar::points_mult_wbits_precompute_sizeof;
    using scalar::points_mult_wbits_scratch_sizeof;
    using scalar::print_bits;
}

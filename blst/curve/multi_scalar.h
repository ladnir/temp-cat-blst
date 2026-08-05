#pragma once

#include "blst/scalar/curve/multi_scalar.h"

namespace blst
{
    /**
     * Selected native C++ MSM facade.
     *
     * Ordinary MSM operands are prevalidated G1/G2/GT subgroup elements and
     * scalars are public. Pippenger accepts at most 256 meaningful
     * little-endian bits and explicit windows in [1, 24]. Fixed-table windows
     * lie in [2, 14]. Caller-owned scratch must be sized with the matching
     * query; empty inputs produce the identity without dereferencing inputs.
     * Compact pointer descriptors use a non-null first entry and interpret a
     * later null entry as contiguous continuation. Exact n-entry descriptor
     * arrays are valid and need no sentinel.
     */
    template<typename T>
    using span = scalar::span<T>;

    using scalar::msm_mult_pippenger;
    using scalar::msm_mult_pippenger_scratch_sizeof_;
    using scalar::msm_mult_pippenger_window_scratch_sizeof_;
    using scalar::msm_mult_table;
    using scalar::msm_precompute_wbits;
    using scalar::pippenger_window_size;
    using scalar::points_mult_wbits_precompute_sizeof;
    using scalar::points_mult_wbits_scratch_sizeof;

#ifdef BLST_CUDA
    using scalar::WindowExtract;
    using scalar::hash_point;
    using scalar::msm_to_jacobian;
#endif
}

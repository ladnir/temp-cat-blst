#pragma once

#include "blst/common/defines.h"
#include "blst/fields/Fp12.h"
#include "blst/curve/point.h"
#if BLST_DEFAULT_AVX
#include "blst/avx/curve/GT.h"
#else
#include "blst/scalar/curve/GT.h"
#endif

namespace blst
{
#if BLST_DEFAULT_AVX
    /**
     * @brief Selected native C++ target-group type (`avx::GT`).
     *
     * Storage does not enforce membership. Arithmetic assumes operands already
     * lie in GT; scalar multiplication documents its scalar and timing contract
     * on the selected implementation type.
     */
    using avx::GT;
#else
    /**
     * @brief Selected native C++ target-group type (`scalar::GT`).
     *
     * Storage does not enforce membership. Arithmetic assumes operands already
     * lie in GT; scalar multiplication documents its scalar and timing contract
     * on the selected implementation type.
     */
    using scalar::GT;
#endif
}

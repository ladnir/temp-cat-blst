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
    using avx::GT;
#else
    using scalar::GT;
#endif
}

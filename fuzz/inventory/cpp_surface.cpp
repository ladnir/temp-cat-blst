// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

// This translation unit has no runtime role. Clang parses it to enumerate the
// reviewed native C++ API owners listed in fuzz/api_inventory.json.

#include "blst/fields/Fr.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fp12.h"
#include "blst/curve/point.h"
#include "blst/curve/GT.h"
#include "blst/curve/map_to_G.h"
#include "blst/curve/multi_scalar.h"
#include "blst/fft/context.h"
#include "blst/fft/domain.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/fr_interpolation.h"
#include "blst/fft/middle_product.h"
#include "blst/sign/aggregate.h"
#include "blst/sign/sha256.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
# include "blst/avx/curve/point.h"
# include "blst/avx/curve/GT.h"
# include "blst/avx/curve/map_to_G.h"
# include "blst/avx/curve/multi_scalar.h"
#endif

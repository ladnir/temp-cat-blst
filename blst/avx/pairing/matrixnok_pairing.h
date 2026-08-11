// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include "miller.h"
#include "pairing.h"
#include "../fields/conversion_inversion.h"
#include "../fields/matrixnok_fp12.h"

namespace blst::avx::matrixnok
{
    using MillerEngine = Miller<FpRing, Fp12Ring>;
    using PairingInverter = FP12RNSVartimeInverter<FpRing, Fp12Ring>;
} // namespace blst::avx::matrixnok

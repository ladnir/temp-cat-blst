// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include "miller.h"
#include "pairing.h"
#include "../fields/matrixnok_fp12.h"

namespace blst::avx::matrixnok
{
    // Isolated production-shaped pairing types for the fixed MatrixNoK field
    // tower. They deliberately do not participate in public dispatch while
    // correctness, code generation, and performance are being evaluated.
    using MillerCandidate = Miller<FpRingCandidate, Fp12RingCandidate>;
    using PairingInverterCandidate =
        FP12RNS_BLSTInverter<FpRingCandidate, Fp12RingCandidate>;
} // namespace blst::avx::matrixnok

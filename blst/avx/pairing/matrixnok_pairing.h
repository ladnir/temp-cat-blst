// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include "miller.h"
#include "pairing.h"
#include "../fields/matrixnok_fp12.h"

namespace blst::avx::matrixnok
{
    // Selected proof candidate for pairing over the fixed MatrixNoK field
    // tower. It remains deliberately absent from public dispatch until the
    // correspondence proof and release gates are complete.
    using MillerCandidate = Miller<FpRingCandidate, Fp12RingCandidate>;
    using PairingInverterCandidate =
        FP12RNS_BLSTInverter<FpRingCandidate, Fp12RingCandidate>;
} // namespace blst::avx::matrixnok

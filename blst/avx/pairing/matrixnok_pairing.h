// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include "miller.h"
#include "pairing.h"
#include "../fields/conversion_inversion.h"
#include "../fields/matrixnok_fp12.h"

namespace blst::avx::matrixnok
{
    // Exact proof-owner spellings for the selected production pairing tower.
    using MillerCandidate = Miller<FpRingCandidate, Fp12RingCandidate>;
    using PairingInverterCandidate =
        FP12RNSVartimeInverter<FpRingCandidate, Fp12RingCandidate>;

    using MillerEngine = MillerCandidate;
    using PairingInverter = PairingInverterCandidate;
} // namespace blst::avx::matrixnok

/*
 * Copyright 2026 Simon Langowski
 * Copyright 2026 Category Labs
 * SPDX-License-Identifier: Apache-2.0 AND MIT
 *
 * Selected aliases for the VROOM-derived Miller construction at commit
 * d6ecd534 and CatBLST's reviewed inversion policy.
 */
#pragma once

#include "miller.h"
#include "../fields/conversion_inversion.h"
#include "../fields/matrixnok_fp12.h"

namespace blst::avx::matrixnok
{
    /**
     * Finite-point production Miller core. Inputs must be validated subgroup
     * affine points and must not be infinity. The engine performs no runtime
     * validation and returns a MillerResult-domain Fp12 value, not GT.
     */
    using MillerEngine = Miller<FpRing, Fp12Ring>;

    /**
     * Variable-time inverter for nonzero public/intermediate Miller values.
     * It must not process secrets whose values may affect inversion timing.
     */
    using PairingInverter = FP12RNSVartimeInverter<FpRing, Fp12Ring>;
} // namespace blst::avx::matrixnok

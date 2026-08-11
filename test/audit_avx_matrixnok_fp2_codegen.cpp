// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// Intentionally small compile-only probes for the selected Fp2 kernel. Keep
// this translation unit independent of the broader MatrixNoK probe so field
// code-generation checks do not instantiate the complete curve and pairing
// stack.

#include "blst/avx/fields/fp2.h"
#include "blst/avx/rns/matrixnok_fp.h"

#if defined(_MSC_VER)
#define AUDIT_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define AUDIT_NOINLINE __attribute__((noinline))
#else
#define AUDIT_NOINLINE
#endif

#if !defined(BLST_AVX512)
#error "MatrixNoK Fp2 code-generation probes require native AVX-512"
#endif

namespace
{
    using FpRing = blst::avx::matrixnok::FpRingCandidate;
    using Fp2Ring = blst::avx::FP2Ring<FpRing>;
    using Fp2 = Fp2Ring::StandardElement;
}

extern "C" AUDIT_NOINLINE void audit_avx_matrixnok_fp2_multiply(
    Fp2* out, const Fp2* left, const Fp2* right)
{
    const FpRing ring;
    const Fp2Ring fp2ring(ring);
    *out = fp2ring.mul(*left, *right);
}

extern "C" AUDIT_NOINLINE void audit_avx_matrixnok_fp2_square(
    Fp2* out, const Fp2* value)
{
    const FpRing ring;
    const Fp2Ring fp2ring(ring);
    *out = fp2ring.square(*value);
}

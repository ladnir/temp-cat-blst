#ifndef BLST_FIELDS_FP4_AVX_H
#define BLST_FIELDS_FP4_AVX_H

/*
 * WARNING: This file is intentionally excluded from the active build during the
 * VROOM backend reset. It is retained only as reference and must not be wired
 * back into active code until the post-reset AVX API is defined.
 */

#include "blst/fields/Fp2_avx.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
    /**
     * AVX512 Fp4 helpers ported from avxbls.
     *
     * Field model:
     *   Fp4 = Fp2[v] / (v^2 - (u + 1))
     *   where Fp2 = Fp[u] / (u^2 + 1).
     *
     * Representation model:
     *   x = x0 + v*x1, with x0, x1 in Fp2.
     *   All components are Montgomery-form Fp2 values.
     */

    /**
     * Two packed Fp4 values in a `4x2x1w`-backed schedule.
     *
     * Let:
     *   X = X.c0 + v*X.c1
     *   Y = Y.c0 + v*Y.c1
     *
     * Pack construction used by tests/bench:
     *   slots = [X.c0, X.c1, Y.c0, Y.c1]
     *
     * Fp2 slot order inside this pack (high -> low):
     *   [Y.c1 | Y.c0 | X.c1 | X.c0]
     *
     * Fp-lane view (high -> low) via inherited `Fp2Pack4x2x1`:
     *   [Y.c1.im | Y.c1.re | Y.c0.im | Y.c0.re | X.c1.im | X.c1.re | X.c0.im | X.c0.re]
     *
     * All squaring variants compute, per element:
     *   out.c0 = in.c0^2 + (u+1)*in.c1^2
     *   out.c1 = 2*in.c0*in.c1
     */
    struct Fp4Pack2x2x2x1 : public Fp2Pack4x2x1
    {
        /// Widened-intermediate schedule matching original avxbls v1 path.
        static void sqr_v1(Fp4Pack2x2x2x1& r, const Fp4Pack2x2x2x1& a);
        /// Single-length schedule matching original avxbls v2 path.
        static void sqr_v2(Fp4Pack2x2x2x1& r, const Fp4Pack2x2x2x1& a);
        /// Karatsuba-style schedule matching original avxbls v3 path.
        static void sqr_v3(Fp4Pack2x2x2x1& r, const Fp4Pack2x2x2x1& a);
    };

    /// Widened companion storage for `Fp4Pack2x2x2x1`.
    struct UFp4Pack2x2x2x1 : public UFp2Pack4x2x1 {};

    /**
     * One packed Fp4 value in a `2x2x2w`-backed schedule.
     *
     * Let:
     *   X = X.c0 + v*X.c1
     *
     * Pack construction used by tests/bench:
     *   slots = [X.c0, X.c1]
     *
     * Fp2 slot order inside this pack (high -> low):
     *   [X.c1 | X.c0]
     *
     * Fp-lane groups follow inherited `Fp2Pack2x2x2` layout.
     *
     * All squaring variants compute:
     *   out.c0 = in.c0^2 + (u+1)*in.c1^2
     *   out.c1 = 2*in.c0*in.c1
     */
    struct Fp4Pack1x2x2x2 : public Fp2Pack2x2x2
    {
        /// Widened-intermediate schedule matching original avxbls v1 path.
        static void sqr_v1(Fp4Pack1x2x2x2& r, const Fp4Pack1x2x2x2& a);
        /// Single-length schedule matching original avxbls v2 path.
        static void sqr_v2(Fp4Pack1x2x2x2& r, const Fp4Pack1x2x2x2& a);
    };

    /// Widened companion storage for `Fp4Pack1x2x2x2`.
    struct UFp4Pack1x2x2x2 : public UFp2Pack2x2x2 {};

    static_assert(std::is_trivially_default_constructible_v<Fp4Pack2x2x2x1>);
    static_assert(std::is_trivially_copyable_v<Fp4Pack2x2x2x1>);
    static_assert(std::is_trivially_default_constructible_v<UFp4Pack2x2x2x1>);
    static_assert(std::is_trivially_copyable_v<UFp4Pack2x2x2x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp4Pack1x2x2x2>);
    static_assert(std::is_trivially_copyable_v<Fp4Pack1x2x2x2>);
    static_assert(std::is_trivially_default_constructible_v<UFp4Pack1x2x2x2>);
    static_assert(std::is_trivially_copyable_v<UFp4Pack1x2x2x2>);

}

#endif

#endif // BLST_FIELDS_FP4_AVX_H

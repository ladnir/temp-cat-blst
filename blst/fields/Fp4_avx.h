#ifndef BLST_FIELDS_FP4_AVX_H
#define BLST_FIELDS_FP4_AVX_H

#include "blst/fields/Fp2_avx.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
    /**
     * AVX512 Fp4 helpers ported from avxbls.
     *
     * Fp4 is represented as `a0 + a1 * v` with `v^2 = (u + 1)` over Fp2.
     * All component values are Montgomery-form Fp2 values.
     */

    /// Two packed Fp4 values in `4x2x1w`-backed schedule.
    using fp4_2x2x2x1w = fp2_4x2x1w;
    /// Widened companion layout for `fp4_2x2x2x1w`.
    using fp4x2_2x2x2x1w = fp2x2_4x2x1w;

    /// One packed Fp4 value in `2x2x2w`-backed schedule.
    using fp4_1x2x2x2w = fp2_2x2x2w;
    /// Widened companion layout for `fp4_1x2x2x2w`.
    using fp4x2_1x2x2x2w = fp2x2_2x2x2w;

    /**
     * Square Fp4 values in `fp4_2x2x2x1w` using widened intermediates.
     *
     * Computes:
     * - `r0 = a0^2 + (u+1)*a1^2`
     * - `r1 = 2*a0*a1`
     *
     * @pre `a` holds reduced Montgomery-form Fp2 components in `2x2x2x1w` schedule.
     * @post `r` holds reduced Montgomery-form Fp2 components in the same schedule.
     */
    void sqr_fp4_2x2x2x1w_v1(fp4_2x2x2x1w r, const fp4_2x2x2x1w a);

    /**
     * Square Fp4 values in `fp4_2x2x2x1w` using single-length intermediates.
     *
     * Formula and layout contract match `sqr_fp4_2x2x2x1w_v1`.
     */
    void sqr_fp4_2x2x2x1w_v2(fp4_2x2x2x1w r, const fp4_2x2x2x1w a);

    /**
     * Square Fp4 values in `fp4_2x2x2x1w` using a Karatsuba-style identity.
     *
     * Computes the same output as `sqr_fp4_2x2x2x1w_v1`/`v2`.
     */
    void sqr_fp4_2x2x2x1w_v3(fp4_2x2x2x1w r, const fp4_2x2x2x1w a);

    /**
     * Square Fp4 values in `fp4_1x2x2x2w` using widened intermediates.
     *
     * Computes:
     * - `r0 = a0^2 + (u+1)*a1^2`
     * - `r1 = 2*a0*a1`
     *
     * @pre `a` holds reduced Montgomery-form Fp2 components in `1x2x2x2w` schedule.
     * @post `r` holds reduced Montgomery-form Fp2 components in the same schedule.
     */
    void sqr_fp4_1x2x2x2w_v1(fp4_1x2x2x2w r, const fp4_1x2x2x2w a);

    /**
     * Square Fp4 values in `fp4_1x2x2x2w` using single-length intermediates.
     *
     * Formula and layout contract match `sqr_fp4_1x2x2x2w_v1`.
     */
    void sqr_fp4_1x2x2x2w_v2(fp4_1x2x2x2w r, const fp4_1x2x2x2w a);
}

#endif

#endif // BLST_FIELDS_FP4_AVX_H

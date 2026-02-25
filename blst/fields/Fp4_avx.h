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
    using Fp4Pack2x2x2x1 = Fp2Pack4x2x1;
    using fp4_2x2x2x1w = Fp4Pack2x2x2x1;
    /// Widened companion layout for `fp4_2x2x2x1w`.
    using UFp4Pack2x2x2x1 = UFp2Pack4x2x1;
    using fp4x2_2x2x2x1w = UFp4Pack2x2x2x1;

    /// One packed Fp4 value in `2x2x2w`-backed schedule.
    using Fp4Pack1x2x2x2 = Fp2Pack2x2x2;
    using fp4_1x2x2x2w = Fp4Pack1x2x2x2;
    /// Widened companion layout for `fp4_1x2x2x2w`.
    using UFp4Pack1x2x2x2 = UFp2Pack2x2x2;
    using fp4x2_1x2x2x2w = UFp4Pack1x2x2x2;

    static_assert(std::is_trivially_default_constructible_v<Fp4Pack2x2x2x1>);
    static_assert(std::is_trivially_copyable_v<Fp4Pack2x2x2x1>);
    static_assert(std::is_trivially_default_constructible_v<UFp4Pack2x2x2x1>);
    static_assert(std::is_trivially_copyable_v<UFp4Pack2x2x2x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp4Pack1x2x2x2>);
    static_assert(std::is_trivially_copyable_v<Fp4Pack1x2x2x2>);
    static_assert(std::is_trivially_default_constructible_v<UFp4Pack1x2x2x2>);
    static_assert(std::is_trivially_copyable_v<UFp4Pack1x2x2x2>);

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
    void sqr_fp4_2x2x2x1w_v1(Fp4Pack2x2x2x1& r, const Fp4Pack2x2x2x1& a);

    /**
     * Square Fp4 values in `fp4_2x2x2x1w` using single-length intermediates.
     *
     * Formula and layout contract match `sqr_fp4_2x2x2x1w_v1`.
     */
    void sqr_fp4_2x2x2x1w_v2(Fp4Pack2x2x2x1& r, const Fp4Pack2x2x2x1& a);

    /**
     * Square Fp4 values in `fp4_2x2x2x1w` using a Karatsuba-style identity.
     *
     * Computes the same output as `sqr_fp4_2x2x2x1w_v1`/`v2`.
     */
    void sqr_fp4_2x2x2x1w_v3(Fp4Pack2x2x2x1& r, const Fp4Pack2x2x2x1& a);

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
    void sqr_fp4_1x2x2x2w_v1(Fp4Pack1x2x2x2& r, const Fp4Pack1x2x2x2& a);

    /**
     * Square Fp4 values in `fp4_1x2x2x2w` using single-length intermediates.
     *
     * Formula and layout contract match `sqr_fp4_1x2x2x2w_v1`.
     */
    void sqr_fp4_1x2x2x2w_v2(Fp4Pack1x2x2x2& r, const Fp4Pack1x2x2x2& a);
}

#endif

#endif // BLST_FIELDS_FP4_AVX_H

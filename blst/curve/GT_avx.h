#ifndef BLST_CURVE_GT_AVX_H
#define BLST_CURVE_GT_AVX_H

#include "blst/curve/GT.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include "blst/fields/Fp12_avx.h"

namespace blst
{
    /**
     * Native SoA GT batch container for exactly 8 independent tuples.
     *
     * Group model:
     * - GT add maps to Fp12 multiply.
     * - GT doubling maps to Fp12 cyclotomic square.
     */
    struct GTPack8x1
    {
        Fp12Pack8x1 c;

        /// Pack 8 scalar GT values (AoS) into SoA.
        static void from_aos(GTPack8x1& r, const GT* in8);

        /// Unpack one GT SoA block into 8 scalar GT values (AoS).
        static void to_aos(GT* out8, const GTPack8x1& a);

        /// Partial AoS pack for tail handling.
        static void from_aos_partial(GTPack8x1& r, const GT* in, std::size_t n_active, const GT& fill);

        /// Partial AoS unpack for tail handling.
        static void to_aos_partial(GT* out, std::size_t n_active, const GTPack8x1& a);

        /// Group add (maps to Fp12 multiply): r = a + b.
        static void add(GTPack8x1& r, const GTPack8x1& a, const GTPack8x1& b);

        /// Group doubling (maps to cyclotomic square): r = 2*a.
        static void doubling(GTPack8x1& r, const GTPack8x1& a);

        /// Lane-wise select between two SoA values.
        static void select(GTPack8x1& r, const GTPack8x1& a, const GTPack8x1& b, u8 lane_mask);
    };
}

#endif

#endif // BLST_CURVE_GT_AVX_H

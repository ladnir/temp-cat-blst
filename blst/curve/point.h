#ifndef __BLS12_381_ASM_POINT_H__
#define __BLS12_381_ASM_POINT_H__

#include "blst/common/defines.h"
#if BLST_DEFAULT_AVX
#include "blst/avx/curve/point.h"
#else
#include "blst/scalar/curve/point.h"
#endif

namespace blst
{
    template<typename Field>
    using AffinePoint = scalar::AffinePoint<Field>;

    template<typename Field>
    using Point = scalar::Point<Field>;

    template<typename Field>
    using Pointxyzz = scalar::Pointxyzz<Field>;

#if BLST_DEFAULT_AVX
    using avx::G1AffinePoint;
    using avx::G1Point;
    using scalar::G1PointZZ;
    using avx::G2AffinePoint;
    using avx::G2Point;
    using scalar::G2PointZZ;
#else
    using scalar::G1AffinePoint;
    using scalar::G1Point;
    using scalar::G1PointZZ;
    using scalar::G2AffinePoint;
    using scalar::G2Point;
    using scalar::G2PointZZ;
#endif

    using scalar::get_B;
#if BLST_DEFAULT_AVX
    using avx::get_BLST_Aprime_E1;
    using avx::get_BLST_Aprime_E2;
#else
    using scalar::get_BLST_Aprime_E1;
    using scalar::get_BLST_Aprime_E2;
#endif
    using scalar::get_BLST_Bprime_E1;
    using scalar::get_BLST_Bprime_E2;
    using scalar::get_BLST_B_E1;
    using scalar::get_BLST_B_E2;
    using scalar::get_BLST_beta;
#if BLST_DEFAULT_AVX
    using avx::get_BLST_G1;
    using avx::get_BLST_G2;
    using avx::get_BLST_NEG_G1;
    using avx::get_BLST_NEG_G2;
#else
    using scalar::get_BLST_G1;
    using scalar::get_BLST_G2;
    using scalar::get_BLST_NEG_G1;
    using scalar::get_BLST_NEG_G2;
#endif

#if BLST_DEFAULT_AVX
    using avx::print;
    using avx::psi;
    using avx::times_minus_z;
#else
    using scalar::print;
    using scalar::psi;
    using scalar::times_minus_z;
#endif
}

#endif

#ifndef __BLS12_381_ASM_POINT_H__
#define __BLS12_381_ASM_POINT_H__

#include "blst/scalar/curve/point.h"

namespace blst
{
    template<typename Field>
    using AffinePoint = scalar::AffinePoint<Field>;

    template<typename Field>
    using Point = scalar::Point<Field>;

    template<typename Field>
    using Pointxyzz = scalar::Pointxyzz<Field>;

    using scalar::G1AffinePoint;
    using scalar::G1Point;
    using scalar::G1PointZZ;
    using scalar::G2AffinePoint;
    using scalar::G2Point;
    using scalar::G2PointZZ;

    using scalar::get_B;
    using scalar::get_BLST_Aprime_E1;
    using scalar::get_BLST_Aprime_E2;
    using scalar::get_BLST_Bprime_E1;
    using scalar::get_BLST_Bprime_E2;
    using scalar::get_BLST_B_E1;
    using scalar::get_BLST_B_E2;
    using scalar::get_BLST_beta;
    using scalar::get_BLST_G1;
    using scalar::get_BLST_G2;
    using scalar::get_BLST_NEG_G1;
    using scalar::get_BLST_NEG_G2;

    using scalar::print;
    using scalar::psi;
    using scalar::times_minus_z;
}

#endif

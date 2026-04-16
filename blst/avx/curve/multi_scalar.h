#pragma once

#include "blst/avx/curve/point.h"
#include "blst/curve/multi_scalar.h"

namespace blst::avx
{
    BLST_HD BLST_INLINE void msm_to_jacobian(G1Point& out, const G1Point& in)
    {
        out = in;
    }

    BLST_HD BLST_INLINE void msm_add_or_sub(
        G1Point& out,
        const G1Point& p1,
        const G1AffinePoint& p2,
        bool_t subtract)
    {
        G1AffinePoint adjusted = p2;
        if (subtract) {
            adjusted.Y = avx_curve_detail::curve_ring().standard_negate(adjusted.Y);
        }
        if (G1Point::is_inf(p1)) {
            G1Point::from_affine(out, adjusted);
            return;
        }
        if (G1AffinePoint::is_inf(adjusted)) {
            out = p1;
            return;
        }
        G1Point::add(out, p1, adjusted);
    }

    template<typename Group, typename AffineGroup, typename GroupZZ>
    BLST_HD BLST_CPP_INLINE void msm_mult_pippenger(
        Group& ret,
        const AffineGroup* const points[],
        size_t npoints,
        const byte* const scalars[],
        size_t nbits,
        GroupZZ buckets[],
        size_t window_size,
        bool print = false)
    {
        blst::scalar::msm_mult_pippenger<Group, AffineGroup, GroupZZ>(
            ret,
            points,
            npoints,
            scalars,
            nbits,
            buckets,
            window_size,
            print);
    }

    template<typename Group>
    BLST_HD BLST_INLINE size_t msm_mult_pippenger_scratch_sizeof_(size_t npoints)
    {
        return sizeof(Group) << (blst::scalar::pippenger_window_size(npoints) - 1);
    }
}

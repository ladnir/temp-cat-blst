#include "blst/fields/Fp12.h"

#include <cooperative_groups.h>
using namespace cooperative_groups;
namespace blst
{
    // ---------------------------------------------------------------
    // Minimal accessors so we can write clean code against your types.
    // Adjust if your operator[] overloads already provide this API.
    // ---------------------------------------------------------------

    __device__ __forceinline__ Fp2& fp12_get_fp2
    (
        Fp12& x,
        int fp12_idx,
        int fp6_idx
    )
    {
        return x.data[fp12_idx].data[fp6_idx];
    }

    __device__ __forceinline__ const Fp2& fp12_get_fp2
    (
        const Fp12& x,
        int fp12_idx,
        int fp6_idx
    )
    {
        return x.data[fp12_idx].data[fp6_idx];
    }

    // ---------------------------------------------------------------
    // Tiny copy helpers for Fp and Fp2 so ptxas can see lifetimes.
    // ---------------------------------------------------------------

    __device__ __forceinline__ void fp_copy
    (
        Fp& dst,
        const Fp& src
    )
    {
#pragma unroll
        for (int i = 0; i < 6; i++)
        {
            dst.data[i] = src.data[i];
        }
    }

    __device__ __forceinline__ void fp2_copy
    (
        Fp2& dst,
        const Fp2& src
    )
    {
        fp_copy(dst.data[0], src.data[0]);
        fp_copy(dst.data[1], src.data[1]);
    }

    // ---------------------------------------------------------------
    // 3-lane shuffles for uint64_t and for a whole Fp (6 limbs) / Fp2.
    // We broadcast from src_role in {0,1,2} to all 3 lanes in the tile.
    // ---------------------------------------------------------------

    __device__ __forceinline__ uint64_t shfl3_u64
    (
        uint64_t v,
        int      src_role,
        unsigned mask
    )
    {
        uint32_t lo;
        uint32_t hi;
        lo = __shfl_sync(mask, static_cast<uint32_t>(v & 0xFFFFFFFFu), src_role);
        hi = __shfl_sync(mask, static_cast<uint32_t>(v >> 32), src_role);
        return (static_cast<uint64_t>(hi) << 32) | static_cast<uint64_t>(lo);
    }

    __device__ __forceinline__ void shfl3_fp
    (
        Fp& dst,
        const Fp& src_from_this_lane,
        int       src_role,
        unsigned  mask
    )
    {
#pragma unroll
        for (int i = 0; i < 6; i++)
        {
            dst.data[i] = shfl3_u64(src_from_this_lane.data[i], src_role, mask);
        }
    }

    __device__ __forceinline__ void shfl3_fp2
    (
        Fp2& dst,
        const Fp2& src_from_this_lane,
        int        src_role,
        unsigned   mask
    )
    {
        shfl3_fp(dst.data[0], src_from_this_lane.data[0], src_role, mask);
        shfl3_fp(dst.data[1], src_from_this_lane.data[1], src_role, mask);
    }

    // ---------------------------------------------------------------
    // Fp4 square using your Fp2 building blocks (branch-free).
    // Input:  a0, a1   Output: out0, out1
    // Matches your non-ASM sqr_fp4 logic:
    //   s0 = a0^2
    //   s1 = a1^2
    //   out0 = (u+1)*s1 + s0
    //   out1 = (a0+a1)^2 - s0 - s1
    // ---------------------------------------------------------------

    __device__ __forceinline__ void sqr_fp4
    (
        Fp2& out0,
        Fp2& out1,
        const Fp2& a0,
        const Fp2& a1
    )
    {
        Fp2 s0;
        Fp2 s1;
        Fp2 sum;
        Fp2 tmp;

        Fp2::sqr(s0, a0);
        Fp2::sqr(s1, a1);
        Fp2::add(sum, a0, a1);
        Fp2::sqr(sum, sum);
        Fp2::mul_by_u_plus_1(out0, s1);
        Fp2::add(out0, out0, s0);
        Fp2::sub(tmp, sum, s0);
        Fp2::sub(out1, tmp, s1);
    }

    // ---------------------------------------------------------------
    // 3-lane cooperative cyclotomic square for a single FP12.
    // Role mapping per 3-lane group (G = lane/3, R = lane%3):
    //   R=0 computes t0 = sqr_fp4(a00, a11)
    //   R=1 computes t1 = sqr_fp4(a10, a02)
    //   R=2 computes t2 = sqr_fp4(a01, a12)
    // Recombine (uniform, branch-free):
    //   ret00 = (t0_0 - a00)*2 + t0_0
    //   ret01 = (t1_0 - a01)*2 + t1_0
    //   ret02 = (t2_0 - a02)*2 + t2_0
    //   tmp   =  mul_by_u_plus_1(t2_1)
    //   ret10 = (tmp + a10)*2 + tmp
    //   ret11 = (t0_1 + a11)*2 + t0_1
    //   ret12 = (t1_1 + a12)*2 + t1_1
    // Stores by role:
    //   R=0 -> out[0][0]=ret00, out[1][1]=ret11
    //   R=1 -> out[0][1]=ret01, out[1][2]=ret12
    //   R=2 -> out[0][2]=ret02, out[1][0]=ret10
    // ---------------------------------------------------------------

    __device__ void cyclotomic_sqr_3lane
    (
        Fp12& out,
        const Fp12& in
    )
    {
        auto tb = this_thread_block();
        auto t3 = tiled_partition<3>(tb);

        const int role = t3.thread_rank();
        const unsigned mask = t3.ballot(true);

        // Gather inputs per role.
        // role 0 owns a00,a11
        // role 1 owns a10,a02
        // role 2 owns a01,a12

        Fp2 a00;
        Fp2 a11;
        Fp2 a10;
        Fp2 a02;
        Fp2 a01;
        Fp2 a12;

        if (role == 0)
        {
            fp2_copy(a00, fp12_get_fp2(in, 0, 0));
            fp2_copy(a11, fp12_get_fp2(in, 1, 1));
        }

        if (role == 1)
        {
            fp2_copy(a10, fp12_get_fp2(in, 1, 0));
            fp2_copy(a02, fp12_get_fp2(in, 0, 2));
        }

        if (role == 2)
        {
            fp2_copy(a01, fp12_get_fp2(in, 0, 1));
            fp2_copy(a12, fp12_get_fp2(in, 1, 2));
        }

        // Compute the three Fp4 squares (one per role).

        Fp2 tX0;
        Fp2 tX1;

        if (role == 0)
        {
            sqr_fp4(tX0, tX1, a00, a11);
        }

        if (role == 1)
        {
            sqr_fp4(tX0, tX1, a10, a02);
        }

        if (role == 2)
        {
            sqr_fp4(tX0, tX1, a01, a12);
        }

        // Broadcast the six temporaries and the six needed inputs to all lanes.

        Fp2 t0_0;
        Fp2 t0_1;
        Fp2 t1_0;
        Fp2 t1_1;
        Fp2 t2_0;
        Fp2 t2_1;

        shfl3_fp2(t0_0, tX0, 0, mask);
        shfl3_fp2(t0_1, tX1, 0, mask);
        shfl3_fp2(t1_0, tX0, 1, mask);
        shfl3_fp2(t1_1, tX1, 1, mask);
        shfl3_fp2(t2_0, tX0, 2, mask);
        shfl3_fp2(t2_1, tX1, 2, mask);

        Fp2 a00_all;
        Fp2 a11_all;
        Fp2 a10_all;
        Fp2 a02_all;
        Fp2 a01_all;
        Fp2 a12_all;

        shfl3_fp2(a00_all, a00, 0, mask);
        shfl3_fp2(a11_all, a11, 0, mask);
        shfl3_fp2(a10_all, a10, 1, mask);
        shfl3_fp2(a02_all, a02, 1, mask);
        shfl3_fp2(a01_all, a01, 2, mask);
        shfl3_fp2(a12_all, a12, 2, mask);

        // Recombine, branch-free. Every lane computes all ret* values locally.

        Fp2 ret00;
        Fp2 ret01;
        Fp2 ret02;
        Fp2 ret10;
        Fp2 ret11;
        Fp2 ret12;

        Fp2 tmp;
        Fp2 x;

        // ret00 = (t0_0 - a00)*2 + t0_0

        Fp2::sub(x, t0_0, a00_all);
        Fp2::add(x, x, x);
        Fp2::add(ret00, x, t0_0);

        // ret01 = (t1_0 - a01)*2 + t1_0

        Fp2::sub(x, t1_0, a01_all);
        Fp2::add(x, x, x);
        Fp2::add(ret01, x, t1_0);

        // ret02 = (t2_0 - a02)*2 + t2_0

        Fp2::sub(x, t2_0, a02_all);
        Fp2::add(x, x, x);
        Fp2::add(ret02, x, t2_0);

        // tmp = mul_by_u_plus_1(t2_1)

        Fp2::mul_by_u_plus_1(tmp, t2_1);

        // ret10 = (tmp + a10)*2 + tmp

        Fp2::add(x, tmp, a10_all);
        Fp2::add(x, x, x);
        Fp2::add(ret10, x, tmp);

        // ret11 = (t0_1 + a11)*2 + t0_1

        Fp2::add(x, t0_1, a11_all);
        Fp2::add(x, x, x);
        Fp2::add(ret11, x, t0_1);

        // ret12 = (t1_1 + a12)*2 + t1_1

        Fp2::add(x, t1_1, a12_all);
        Fp2::add(x, x, x);
        Fp2::add(ret12, x, t1_1);

        // Stores by role. Each lane writes a disjoint pair, so no conflicts.

        if (role == 0)
        {
            fp2_copy(fp12_get_fp2(out, 0, 0), ret00);
            fp2_copy(fp12_get_fp2(out, 1, 1), ret11);
        }

        if (role == 1)
        {
            fp2_copy(fp12_get_fp2(out, 0, 1), ret01);
            fp2_copy(fp12_get_fp2(out, 1, 2), ret12);
        }

        if (role == 2)
        {
            fp2_copy(fp12_get_fp2(out, 0, 2), ret02);
            fp2_copy(fp12_get_fp2(out, 1, 0), ret10);
        }
    }
}
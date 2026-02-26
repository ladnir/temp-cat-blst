#include "blst/curve/GT_avx.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include <type_traits>

namespace blst
{
namespace
{
    static_assert(std::is_standard_layout_v<GT>);
    static_assert(sizeof(GT) == sizeof(Fp12));
    static_assert(alignof(GT) == alignof(Fp12));

    BLST_FORCEINLINE const Fp12* gt_as_fp12_ptr(const GT* in)
    {
        return reinterpret_cast<const Fp12*>(in);
    }

    BLST_FORCEINLINE Fp12* gt_as_fp12_ptr(GT* in)
    {
        return reinterpret_cast<Fp12*>(in);
    }
} // namespace

void GTPack8x1::from_aos(GTPack8x1& r, const GT* in8)
{
    BLST_ASSERT(in8 != nullptr);
    Fp12Pack8x1::from_aos(r.c, gt_as_fp12_ptr(in8));
}

void GTPack8x1::to_aos(GT* out8, const GTPack8x1& a)
{
    BLST_ASSERT(out8 != nullptr);
    Fp12Pack8x1::to_aos(gt_as_fp12_ptr(out8), a.c);
}

void GTPack8x1::from_aos_partial(
    GTPack8x1& r,
    const GT* in,
    const std::size_t n_active,
    const GT& fill)
{
    BLST_ASSERT(n_active <= 8);
    BLST_ASSERT(in != nullptr || n_active == 0);
    Fp12Pack8x1::from_aos_partial(r.c, gt_as_fp12_ptr(in), n_active, fill.data);
}

void GTPack8x1::to_aos_partial(GT* out, const std::size_t n_active, const GTPack8x1& a)
{
    BLST_ASSERT(n_active <= 8);
    BLST_ASSERT(out != nullptr || n_active == 0);
    Fp12Pack8x1::to_aos_partial(gt_as_fp12_ptr(out), n_active, a.c);
}

void GTPack8x1::add(GTPack8x1& r, const GTPack8x1& a, const GTPack8x1& b)
{
    Fp12Pack8x1::mul(r.c, a.c, b.c);
}

void GTPack8x1::doubling(GTPack8x1& r, const GTPack8x1& a)
{
    Fp12Pack8x1::cyclotomic_sqr(r.c, a.c);
}

void GTPack8x1::select(GTPack8x1& r, const GTPack8x1& a, const GTPack8x1& b, const u8 lane_mask)
{
    Fp12Pack8x1::select(r.c, a.c, b.c, lane_mask);
}

} // namespace blst

#endif

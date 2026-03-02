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

    BLST_FORCEINLINE u8 scalar_lane_mask(
        const byte* scalars8,
        const std::size_t scalar_stride,
        const std::size_t bit_index)
    {
        const std::size_t byte_index = bit_index >> 3;
        const std::size_t bit_in_byte = bit_index & 7;
        u8 lane_mask = 0;
        for (std::size_t lane = 0; lane < 8; ++lane) {
            const byte scalar_byte = scalars8[lane * scalar_stride + byte_index];
            lane_mask |= static_cast<u8>(((scalar_byte >> bit_in_byte) & 1u) << lane);
        }
        return lane_mask;
    }

    BLST_FORCEINLINE const GTPack8x1& gt_identity_pack()
    {
        static const GTPack8x1 kIdentityPack = []() {
            GTPack8x1 id{};
            GT lanes[8];
            for (std::size_t i = 0; i < 8; ++i) {
                lanes[i] = GT::identity();
            }
            GTPack8x1::from_aos(id, lanes);
            return id;
        }();
        return kIdentityPack;
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

void GTPack8x1::mult(
    GTPack8x1& out,
    const GTPack8x1& base,
    const byte* scalars8,
    const std::size_t scalar_stride,
    const std::size_t nbits)
{
    const std::size_t required_scalar_bytes = (nbits + 7) >> 3;
    BLST_ASSERT(scalars8 != nullptr || nbits == 0);
    BLST_ASSERT((scalar_stride > 0) || (nbits == 0));
    BLST_ASSERT((nbits == 0) || (scalar_stride >= required_scalar_bytes));

    if ((nbits != 0) && ((scalars8 == nullptr) || (scalar_stride < required_scalar_bytes))) {
        return;
    }

    const GTPack8x1 base_copy = base;
    out = gt_identity_pack();

    if (nbits == 0) {
        return;
    }

    GTPack8x1 add_candidate{};

    for (std::size_t bit = nbits; bit-- > 0;) {
        GTPack8x1::doubling(out, out);
        GTPack8x1::add(add_candidate, out, base_copy);
        const u8 lane_mask = scalar_lane_mask(scalars8, scalar_stride, bit);
        GTPack8x1::select(out, out, add_candidate, lane_mask);
    }
}

void GTPack8x1::add(
    std::span<GTPack8x1> out,
    std::span<const GTPack8x1> a,
    std::span<const GTPack8x1> b)
{
    if (out.empty()) {
        return;
    }

    const std::size_t n_packs = out.size();
    BLST_ASSERT(a.size() == n_packs);
    BLST_ASSERT(b.size() == n_packs);

    if ((a.size() != n_packs) || (b.size() != n_packs)) {
        return;
    }

    for (std::size_t i = 0; i < n_packs; ++i) {
        GTPack8x1::add(out[i], a[i], b[i]);
    }
}

void GTPack8x1::doubling(
    std::span<GTPack8x1> out,
    std::span<const GTPack8x1> a)
{
    if (out.empty()) {
        return;
    }

    const std::size_t n_packs = out.size();
    BLST_ASSERT(a.size() == n_packs);

    if (a.size() != n_packs) {
        return;
    }

    for (std::size_t i = 0; i < n_packs; ++i) {
        GTPack8x1::doubling(out[i], a[i]);
    }
}

void GTPack8x1::mult(
    std::span<GTPack8x1> out,
    std::span<const GTPack8x1> base,
    std::span<const byte> scalars,
    const std::size_t scalar_stride,
    const std::size_t nbits)
{
    const std::size_t required_scalar_bytes = (nbits + 7) >> 3;
    BLST_ASSERT((scalar_stride > 0) || (nbits == 0));
    BLST_ASSERT((nbits == 0) || (scalar_stride >= required_scalar_bytes));

    if (out.empty()) {
        return;
    }

    if ((nbits != 0) && (scalar_stride < required_scalar_bytes)) {
        return;
    }

    const std::size_t n_packs = out.size();
    BLST_ASSERT(base.size() == n_packs);

    if (base.size() != n_packs) {
        return;
    }

    const std::size_t required_scalar_span = (nbits == 0) ? 0 : (n_packs * 8 * scalar_stride);
    BLST_ASSERT(scalars.size() >= required_scalar_span);
    if (scalars.size() < required_scalar_span) {
        return;
    }

    for (std::size_t i = 0; i < n_packs; ++i) {
        const byte* pack_scalars = (nbits == 0) ? nullptr : (scalars.data() + (i * 8 * scalar_stride));
        GTPack8x1::mult(out[i], base[i], pack_scalars, scalar_stride, nbits);
    }
}

void GTPack8x1::select(GTPack8x1& r, const GTPack8x1& a, const GTPack8x1& b, const u8 lane_mask)
{
    Fp12Pack8x1::select(r.c, a.c, b.c, lane_mask);
}

} // namespace blst

#endif

#include "blst/curve/GT_avx.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include <array>
#include <type_traits>

namespace blst
{
namespace
{
    constexpr std::size_t kGtPackLanes = 8;
    constexpr std::size_t kGtGlv4Bits = 256;
    constexpr std::size_t kGtGlv4Limbs = 4;
    constexpr std::size_t kGtGlv4MaxSteps = 63;

    struct LaneGlv4Plan
    {
        u8 init_idx = 0;
        u8 total_steps = 0;
        std::array<u8, kGtGlv4MaxSteps + 1> mul_idx_at_step{};
    };

    static_assert(std::is_standard_layout_v<GT>);
    static_assert(sizeof(GT) == sizeof(Fp12));
    static_assert(alignof(GT) == alignof(Fp12));
    static_assert(BLST_LIMB_SIZE == 64);

    BLST_FORCEINLINE const Fp12* gt_as_fp12_ptr(const GT* in)
    {
        return reinterpret_cast<const Fp12*>(in);
    }

    BLST_FORCEINLINE Fp12* gt_as_fp12_ptr(GT* in)
    {
        return reinterpret_cast<Fp12*>(in);
    }

    BLST_FORCEINLINE limb_t load_le_limb(const byte* in)
    {
        limb_t v = 0;
        for (std::size_t i = 0; i < sizeof(limb_t); ++i) {
            v |= static_cast<limb_t>(in[i]) << (8 * i);
        }
        return v;
    }

    BLST_FORCEINLINE void load_scalar_limbs_256(limb_t out[kGtGlv4Limbs], const byte* scalar_32)
    {
        for (std::size_t i = 0; i < kGtGlv4Limbs; ++i) {
            out[i] = load_le_limb(scalar_32 + i * sizeof(limb_t));
        }
    }

    BLST_FORCEINLINE int scalar_lane_msb(
        const byte* scalars8,
        const std::size_t scalar_stride,
        const std::size_t lane,
        const std::size_t nbits)
    {
        if ((scalars8 == nullptr) || (nbits == 0)) {
            return -1;
        }

        const byte* lane_scalar = scalars8 + lane * scalar_stride;
        for (std::size_t bit = nbits; bit-- > 0;) {
            const std::size_t byte_idx = bit >> 3;
            const std::size_t bit_idx = bit & 7;
            if (((lane_scalar[byte_idx] >> bit_idx) & 1u) != 0) {
                return static_cast<int>(bit);
            }
        }
        return -1;
    }

    BLST_FORCEINLINE bool scalar_ge_group_order_256(const limb_t k[kGtGlv4Limbs])
    {
        const Fr& order = get_BLST_r();
        for (int limb = static_cast<int>(kGtGlv4Limbs) - 1; limb >= 0; --limb) {
            if (k[limb] > order[limb]) {
                return true;
            }
            if (k[limb] < order[limb]) {
                return false;
            }
        }
        return true;
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

    BLST_FORCEINLINE LaneGlv4Plan build_lane_glv4_plan(const GLV4_Exp& props, bool& ok)
    {
        LaneGlv4Plan plan{};
        if (props.sched_len == 0) {
            return plan;
        }

        plan.init_idx = static_cast<u8>(props.schedule[0] & 0x0F);
        std::size_t step = 0;
        for (std::size_t i = 1; i < props.sched_len; ++i) {
            const u8 op = props.schedule[i];
            const std::size_t sqr_cnt = op >> 4;
            const u8 mul_idx = static_cast<u8>(op & 0x0F);
            step += sqr_cnt;
            if (step > kGtGlv4MaxSteps) {
                ok = false;
                return LaneGlv4Plan{};
            }
            if (mul_idx != 0) {
                if (plan.mul_idx_at_step[step] != 0) {
                    ok = false;
                    return LaneGlv4Plan{};
                }
                plan.mul_idx_at_step[step] = mul_idx;
            }
        }
        plan.total_steps = static_cast<u8>(step);
        return plan;
    }

    BLST_FORCEINLINE bool gt_compile_glv4_plan(
        GTGlv4Pack8x1& plan,
        const byte* scalars8,
        const std::size_t scalar_stride,
        const std::size_t nbits)
    {
        BLST_ASSERT((nbits > 0) && (nbits <= kGtGlv4Bits));
        BLST_ASSERT((scalars8 != nullptr) && (scalar_stride >= ((nbits + 7) >> 3)));

        const std::size_t byte_len = (nbits + 7) >> 3;
        const std::size_t rem_bits = nbits & 7u;
        std::array<std::array<byte, kGtGlv4Bits / 8>, kGtPackLanes> scalar_bytes{};
        std::array<std::array<limb_t, kGtGlv4Limbs>, kGtPackLanes> scalar_limbs{};

        bool ok = true;
        std::size_t max_steps = 0;
        for (std::size_t lane = 0; lane < kGtPackLanes; ++lane) {
            const byte* lane_scalar = scalars8 + lane * scalar_stride;
            for (std::size_t i = 0; i < byte_len; ++i) {
                scalar_bytes[lane][i] = lane_scalar[i];
            }
            if ((rem_bits != 0) && (byte_len > 0)) {
                const byte mask = static_cast<byte>((1u << rem_bits) - 1u);
                scalar_bytes[lane][byte_len - 1] &= mask;
            }

            load_scalar_limbs_256(scalar_limbs[lane].data(), scalar_bytes[lane].data());
#ifndef NDEBUG
            BLST_ASSERT(!scalar_ge_group_order_256(scalar_limbs[lane].data()));
#endif
            plan.props[lane] = Fp12::glv4_compile(scalar_limbs[lane].data());
            const LaneGlv4Plan lane_plan = build_lane_glv4_plan(plan.props[lane], ok);
            if (!ok) {
                return false;
            }

            plan.init_idx[lane] = lane_plan.init_idx;
            plan.total_steps[lane] = lane_plan.total_steps;
            plan.mul_idx_at_step[lane] = lane_plan.mul_idx_at_step;
            if (lane_plan.total_steps > max_steps) {
                max_steps = lane_plan.total_steps;
            }
        }

        plan.max_steps = static_cast<u8>(max_steps);
        return true;
    }

    BLST_FORCEINLINE void gt_mult_basic(
        GTPack8x1& out,
        const GTPack8x1& base,
        const byte* scalars8,
        const std::size_t scalar_stride,
        const std::size_t nbits)
    {
        const GTPack8x1 base_copy = base;
        out = gt_identity_pack();

        if (nbits == 0) {
            return;
        }

        std::array<int, kGtPackLanes> lane_msb{};
        int max_msb = -1;
        for (std::size_t lane = 0; lane < kGtPackLanes; ++lane) {
            lane_msb[lane] = scalar_lane_msb(scalars8, scalar_stride, lane, nbits);
            if (lane_msb[lane] > max_msb) {
                max_msb = lane_msb[lane];
            }
        }
        if (max_msb < 0) {
            return;
        }

        GTPack8x1 sqr_candidate{};
        GTPack8x1 add_candidate{};
        for (int bit = max_msb; bit >= 0; --bit) {
            u8 active_mask = 0;
            for (std::size_t lane = 0; lane < kGtPackLanes; ++lane) {
                active_mask |= static_cast<u8>((lane_msb[lane] >= bit ? 1u : 0u) << lane);
            }

            if (active_mask != 0) {
                GTPack8x1::doubling(sqr_candidate, out);
                GTPack8x1::select(out, out, sqr_candidate, active_mask);
            }

            const u8 lane_mask =
                static_cast<u8>(scalar_lane_mask(scalars8, scalar_stride, static_cast<std::size_t>(bit)) &
                                active_mask);
            if (lane_mask != 0) {
                GTPack8x1::add(add_candidate, out, base_copy);
                GTPack8x1::select(out, out, add_candidate, lane_mask);
            }
        }
    }

    BLST_FORCEINLINE void gt_mult_glv4_precompiled(
        GTPack8x1& out,
        const GTPack8x1& base,
        const GTGlv4Pack8x1& plan)
    {
        std::array<Fp12Pack8x1, 4> basis{};
        basis[0] = base.c;
        Fp12Pack8x1::frobenius_map(basis[1], base.c, 1);
        Fp12Pack8x1::frobenius_map(basis[2], base.c, 2);
        Fp12Pack8x1::frobenius_map(basis[3], base.c, 3);

        for (std::size_t basis_idx = 0; basis_idx < basis.size(); ++basis_idx) {
            u8 conj_mask = 0;
            for (std::size_t lane = 0; lane < kGtPackLanes; ++lane) {
                const u8 sign_bit = static_cast<u8>((plan.props[lane].sign_mask >> basis_idx) & 1u);
                conj_mask |= static_cast<u8>(sign_bit << lane);
            }
            if (conj_mask != 0) {
                Fp12Pack8x1 conj{};
                Fp12Pack8x1::conjugate(conj, basis[basis_idx]);
                Fp12Pack8x1::select(basis[basis_idx], basis[basis_idx], conj, conj_mask);
            }
        }

        std::array<Fp12Pack8x1, 16> table{};
        table[0] = gt_identity_pack().c;
        table[1] = basis[0];
        table[2] = basis[1];
        table[4] = basis[2];
        table[8] = basis[3];

        Fp12Pack8x1::mul(table[3], table[1], table[2]);
        Fp12Pack8x1::mul(table[5], table[1], table[4]);
        Fp12Pack8x1::mul(table[6], table[2], table[4]);
        Fp12Pack8x1::mul(table[7], table[3], table[4]);

        Fp12Pack8x1::mul(table[9], table[1], table[8]);
        Fp12Pack8x1::mul(table[10], table[2], table[8]);
        Fp12Pack8x1::mul(table[11], table[3], table[8]);
        Fp12Pack8x1::mul(table[12], table[4], table[8]);
        Fp12Pack8x1::mul(table[13], table[5], table[8]);
        Fp12Pack8x1::mul(table[14], table[6], table[8]);
        Fp12Pack8x1::mul(table[15], table[7], table[8]);

        out = gt_identity_pack();
        for (u8 idx = 1; idx < 16; ++idx) {
            u8 lane_mask = 0;
            for (std::size_t lane = 0; lane < kGtPackLanes; ++lane) {
                lane_mask |= static_cast<u8>((plan.init_idx[lane] == idx ? 1u : 0u) << lane);
            }
            if (lane_mask != 0) {
                Fp12Pack8x1::select(out.c, out.c, table[idx], lane_mask);
            }
        }

        Fp12Pack8x1 sqr_candidate{};
        Fp12Pack8x1 mul_candidate{};
        Fp12Pack8x1 mul_operand{};
        for (std::size_t step = 1; step <= plan.max_steps; ++step) {
            u8 sqr_mask = 0;
            std::array<u8, 16> mul_masks{};
            for (std::size_t lane = 0; lane < kGtPackLanes; ++lane) {
                const bool lane_active =
                    (plan.init_idx[lane] != 0) && (step <= static_cast<std::size_t>(plan.total_steps[lane]));
                sqr_mask |= static_cast<u8>((lane_active ? 1u : 0u) << lane);

                if (lane_active) {
                    const u8 idx = plan.mul_idx_at_step[lane][step];
                    mul_masks[idx] |= static_cast<u8>(1u << lane);
                }
            }
            if (sqr_mask != 0) {
                Fp12Pack8x1::cyclotomic_sqr(sqr_candidate, out.c);
                Fp12Pack8x1::select(out.c, out.c, sqr_candidate, sqr_mask);
            }

            u8 nonzero_mul_mask = 0;
            for (u8 idx = 1; idx < 16; ++idx) {
                nonzero_mul_mask |= mul_masks[idx];
            }
            if (nonzero_mul_mask != 0) {
                mul_operand = table[0];
                for (u8 idx = 1; idx < 16; ++idx) {
                    if (mul_masks[idx] != 0) {
                        Fp12Pack8x1::select(mul_operand, mul_operand, table[idx], mul_masks[idx]);
                    }
                }
                Fp12Pack8x1::mul(mul_candidate, out.c, mul_operand);
                out.c = mul_candidate;
            }
        }
    }

    BLST_FORCEINLINE bool gt_mult_glv4(
        GTPack8x1& out,
        const GTPack8x1& base,
        const byte* scalars8,
        const std::size_t scalar_stride,
        const std::size_t nbits)
    {
        GTGlv4Pack8x1 plan{};
        if (!gt_compile_glv4_plan(plan, scalars8, scalar_stride, nbits)) {
            return false;
        }
        gt_mult_glv4_precompiled(out, base, plan);
        return true;
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

    if ((nbits > 128) && (nbits <= kGtGlv4Bits)) {
        if (gt_mult_glv4(out, base, scalars8, scalar_stride, nbits)) {
            return;
        }
    }

    gt_mult_basic(out, base, scalars8, scalar_stride, nbits);
}

bool GTPack8x1::compile_glv4(
    GTGlv4Pack8x1& out_plan,
    const byte* scalars8,
    const std::size_t scalar_stride,
    const std::size_t nbits)
{
    const std::size_t required_scalar_bytes = (nbits + 7) >> 3;
    BLST_ASSERT((nbits > 0) && (nbits <= kGtGlv4Bits));
    BLST_ASSERT(scalars8 != nullptr);
    BLST_ASSERT(scalar_stride >= required_scalar_bytes);

    if ((nbits == 0) || (nbits > kGtGlv4Bits) || (scalars8 == nullptr) || (scalar_stride < required_scalar_bytes)) {
        return false;
    }

    return gt_compile_glv4_plan(out_plan, scalars8, scalar_stride, nbits);
}

void GTPack8x1::mult(
    GTPack8x1& out,
    const GTPack8x1& base,
    const GTGlv4Pack8x1& plan)
{
    gt_mult_glv4_precompiled(out, base, plan);
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

bool GTPack8x1::compile_glv4(
    std::span<GTGlv4Pack8x1> out_plan,
    std::span<const byte> scalars,
    const std::size_t scalar_stride,
    const std::size_t nbits)
{
    const std::size_t required_scalar_bytes = (nbits + 7) >> 3;
    BLST_ASSERT((nbits > 0) && (nbits <= kGtGlv4Bits));
    BLST_ASSERT(scalar_stride >= required_scalar_bytes);

    if ((nbits == 0) || (nbits > kGtGlv4Bits) || (scalar_stride < required_scalar_bytes)) {
        return false;
    }

    if (out_plan.empty()) {
        return true;
    }

    const std::size_t n_packs = out_plan.size();
    const std::size_t required_scalar_span = n_packs * 8 * scalar_stride;
    BLST_ASSERT(scalars.size() >= required_scalar_span);
    if (scalars.size() < required_scalar_span) {
        return false;
    }

    for (std::size_t i = 0; i < n_packs; ++i) {
        const byte* pack_scalars = scalars.data() + (i * 8 * scalar_stride);
        if (!GTPack8x1::compile_glv4(out_plan[i], pack_scalars, scalar_stride, nbits)) {
            return false;
        }
    }
    return true;
}

void GTPack8x1::mult(
    std::span<GTPack8x1> out,
    std::span<const GTPack8x1> base,
    std::span<const GTGlv4Pack8x1> plan)
{
    if (out.empty()) {
        return;
    }

    const std::size_t n_packs = out.size();
    BLST_ASSERT(base.size() == n_packs);
    BLST_ASSERT(plan.size() == n_packs);
    if ((base.size() != n_packs) || (plan.size() != n_packs)) {
        return;
    }

    for (std::size_t i = 0; i < n_packs; ++i) {
        GTPack8x1::mult(out[i], base[i], plan[i]);
    }
}

void GTPack8x1::select(GTPack8x1& r, const GTPack8x1& a, const GTPack8x1& b, const u8 lane_mask)
{
    Fp12Pack8x1::select(r.c, a.c, b.c, lane_mask);
}

} // namespace blst

#endif

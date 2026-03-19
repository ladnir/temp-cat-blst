#pragma once

#include "blst/common/defines.h"
#include "blst/curve/GT.h"
#include "blst/fields/Fr.h"
#include "blst/vroom/src/bounded_ring.hpp"
#include "blst/vroom/src/conversion_inversion.hpp"
#include "blst/vroom/src/fp12.hpp"
#include "blst/vroom/src/gt_glv4.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace blst {

inline constexpr char kVroomGt52BaseModulus[] =
    "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787";

struct VroomGt52
{
    using Ring = BoundedRing<>;
    using Fp12Ring = FP12Ring<Ring>;
    using Element = typename Fp12Ring::StandardElement;

    Element data = {};

    VroomGt52() = default;
    explicit VroomGt52(const Element& value) : data(value) {}

    static const Ring& ring()
    {
        static const Ring instance(BigInt(kVroomGt52BaseModulus, 10));
        return instance;
    }

    static const Fp12Ring& fp12ring()
    {
        static const Fp12Ring instance(ring());
        return instance;
    }

    static VroomGt52 identity()
    {
        return VroomGt52(fp12ring().one(ring()));
    }

    static VroomGt52 from_scalar(const GT& value)
    {
        return VroomGt52(vroom_blst_direct_bridge::scalar_fp12_to_ring(ring(), value.data));
    }

    static GT to_scalar(const VroomGt52& value)
    {
        GT out = {};
        out.data = vroom_blst_direct_bridge::ring_fp12_to_scalar(ring(), value.data);
        return out;
    }
};

inline std::array<uint64_t, 4> fr_to_limbs_le(const Fr& value)
{
    Fr canonical = {};
    Fr::from_mont(canonical, value);
    return {
        canonical.data[0],
        canonical.data[1],
        canonical.data[2],
        canonical.data[3],
    };
}

struct VroomGt52PreparedRoots
{
    std::span<const Fr> scalars;
    std::vector<std::array<uint64_t, 4>> limbs;
    std::vector<VroomGLV4Exp> compiled;
    Fr scale_scalar = {};
    std::array<uint64_t, 4> scale_limbs = {};
    VroomGLV4Exp scale_compiled = {};
    bool has_scale = false;

    explicit VroomGt52PreparedRoots(std::span<const Fr> input_scalars = {})
        : scalars(input_scalars)
    {
        limbs.reserve(scalars.size());
        compiled.reserve(scalars.size());
        for (const Fr& scalar : scalars) {
            limbs.push_back(fr_to_limbs_le(scalar));
            compiled.push_back(vroom_glv4_compile(limbs.back().data()));
        }
    }

    void set_scale_from_size(std::size_t n)
    {
        BLST_ASSERT(n != 0);
        Fr::from_uint64(scale_scalar, static_cast<uint64_t>(n));
        Fr::reciprocal(scale_scalar, scale_scalar);
        scale_limbs = fr_to_limbs_le(scale_scalar);
        scale_compiled = vroom_glv4_compile(scale_limbs.data());
        has_scale = true;
    }
};

struct VroomGt52CoeffCtx
{
    class ScopedRoots
    {
    public:
        explicit ScopedRoots(const VroomGt52PreparedRoots& roots)
            : previous_(prepared_roots_)
        {
            prepared_roots_ = &roots;
        }

        ~ScopedRoots()
        {
            prepared_roots_ = previous_;
        }

        ScopedRoots(const ScopedRoots&) = delete;
        ScopedRoots& operator=(const ScopedRoots&) = delete;

    private:
        const VroomGt52PreparedRoots* previous_;
    };

    static auto plus(const VroomGt52& a, const VroomGt52& b)
    {
        return VroomGt52(VroomGt52::fp12ring().fp12_flat_mul(a.data, b.data, VroomGt52::ring()));
    }

    static auto minus(const VroomGt52& a, const VroomGt52& b)
    {
        const auto b_inv = VroomGt52::fp12ring().fp12_flat_conjugate(b.data, VroomGt52::ring());
        return VroomGt52(VroomGt52::fp12ring().fp12_flat_mul(a.data, b_inv, VroomGt52::ring()));
    }

    static auto mult(const VroomGt52& a, const Fr& b)
    {
        const auto& ring = VroomGt52::ring();
        const auto& fp12ring = VroomGt52::fp12ring();

        if (prepared_roots_ != nullptr && !prepared_roots_->scalars.empty()) {
            const Fr* base = prepared_roots_->scalars.data();
            if (&b >= base && &b < base + prepared_roots_->scalars.size()) {
                const std::size_t idx = static_cast<std::size_t>(&b - base);
                return VroomGt52(vroom_cyclotomic_exp_glv4(
                    a.data,
                    prepared_roots_->compiled[idx],
                    fp12ring,
                    ring));
            }
        }

        if (prepared_roots_ != nullptr &&
            prepared_roots_->has_scale &&
            Fr::equal(b, prepared_roots_->scale_scalar)) {
            return VroomGt52(vroom_cyclotomic_exp_glv4(
                a.data,
                prepared_roots_->scale_compiled,
                fp12ring,
                ring));
        }

        const std::array<uint64_t, 4> scalar_limbs = fr_to_limbs_le(b);
        return VroomGt52(vroom_cyclotomic_exp_glv4(a.data, scalar_limbs.data(), fp12ring, ring));
    }

    static auto inverse(const Fr& a)
    {
        Fr r = {};
        Fr::reciprocal(r, a);
        return r;
    }

    template<typename F>
    static auto from_u64(uint64_t a)
    {
        Fr r = {};
        Fr::from_uint64(r, a);
        return r;
    }

private:
    inline static thread_local const VroomGt52PreparedRoots* prepared_roots_ = nullptr;
};

inline std::vector<VroomGt52> scalar_gt_vector_to_vroom52(std::span<const GT> values)
{
    std::vector<VroomGt52> out;
    out.reserve(values.size());
    for (const GT& value : values) {
        out.push_back(VroomGt52::from_scalar(value));
    }
    return out;
}

inline std::vector<GT> vroom52_vector_to_scalar_gt(std::span<const VroomGt52> values)
{
    std::vector<GT> out;
    out.reserve(values.size());
    for (const VroomGt52& value : values) {
        out.push_back(VroomGt52::to_scalar(value));
    }
    return out;
}

}  // namespace blst

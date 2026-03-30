#include "bindings/blst_middle_product.h"

#include "blst/fft/context.h"
#include "blst/fft/middle_product.h"

#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using BoundaryG1 = blst::scalar::G1Point;
using BoundaryG2 = blst::scalar::G2Point;
using BoundaryFp12 = blst::Fp12;
using PublicGt = blst::GT;

struct BlstP1Layout
{
    blst::Fp x;
    blst::Fp y;
    blst::Fp z;
};

struct BlstP2Layout
{
    blst::Fp2 x;
    blst::Fp2 y;
    blst::Fp2 z;
};

struct BlstFp12Layout
{
    blst::Fp6 fp6[2];
};

static_assert(sizeof(BlstP1Layout) == sizeof(BoundaryG1));
static_assert(alignof(BlstP1Layout) == alignof(BoundaryG1));
static_assert(sizeof(BlstP2Layout) == sizeof(BoundaryG2));
static_assert(alignof(BlstP2Layout) == alignof(BoundaryG2));
static_assert(sizeof(BlstFp12Layout) == sizeof(BoundaryFp12));
static_assert(alignof(BlstFp12Layout) == alignof(BoundaryFp12));
static_assert(std::is_standard_layout_v<BoundaryG1>);
static_assert(std::is_standard_layout_v<BoundaryG2>);
static_assert(std::is_standard_layout_v<BoundaryFp12>);
static_assert(std::is_trivially_copyable_v<BoundaryG1>);
static_assert(std::is_trivially_copyable_v<BoundaryG2>);
static_assert(std::is_trivially_copyable_v<BoundaryFp12>);

inline BoundaryFp12 gt_to_boundary(const PublicGt& value)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    {
        unsigned char serialized[12 * 48] = {};
        PublicGt::serialize(serialized, value);
        blst::scalar::GT scalar_gt = {};
        blst::scalar::GT::deserialize(scalar_gt, serialized);
        return scalar_gt.data;
    }
#else
    if constexpr (std::is_same_v<PublicGt, blst::scalar::GT>) {
        return value.data;
    } else {
        unsigned char serialized[12 * 48] = {};
        PublicGt::serialize(serialized, value);
        blst::scalar::GT scalar_gt = {};
        blst::scalar::GT::deserialize(scalar_gt, serialized);
        return scalar_gt.data;
    }
#endif
}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#define BLST_MIDDLE_PRODUCT_USE_AVX 1
using NativeG1 = blst::avx::G1Point;
using NativeG2 = blst::avx::G2Point;
using NativePrecomputedB = blst::MiddleProductPrecomputedBAvx;
using NativeWorkspace = blst::MiddleProductWorkspaceAvx;

inline std::vector<NativeG1> convert_g1_inputs(std::span<const BoundaryG1> input)
{
    std::vector<NativeG1> out;
    out.reserve(input.size());
    for (const BoundaryG1& point : input) {
        out.push_back(blst::avx::avx_map_to_G_detail::from_scalar(point));
    }
    return out;
}

inline std::vector<NativeG2> convert_g2_inputs(std::span<const BoundaryG2> input)
{
    std::vector<NativeG2> out;
    out.reserve(input.size());
    for (const BoundaryG2& point : input) {
        out.push_back(blst::avx::avx_map_to_G_detail::from_scalar(point));
    }
    return out;
}

inline NativePrecomputedB precompute_native(
    std::span<const BoundaryG2> b,
    std::size_t n)
{
    const std::vector<NativeG2> converted = convert_g2_inputs(b);
    return blst::detail::precompute_middle_product_b_impl<blst::MiddleProductAvxContext>(
        std::span<const NativeG2>(converted),
        n,
        nullptr,
        blst::MiddleProductPairingMode::Direct);
}

inline void run_native(
    std::span<const BoundaryG1> a,
    const NativePrecomputedB& precomputed_b,
    std::span<PublicGt> out)
{
    const std::vector<NativeG1> converted = convert_g1_inputs(a);
    NativeWorkspace workspace;
    blst::detail::middle_product_impl<blst::MiddleProductAvxContext>(
        std::span<const NativeG1>(converted),
        precomputed_b,
        out,
        workspace,
        nullptr);
}
#else
#define BLST_MIDDLE_PRODUCT_USE_AVX 0
using NativePrecomputedB = blst::MiddleProductPrecomputedBScalar;
using NativeWorkspace = blst::MiddleProductWorkspaceScalar;

inline NativePrecomputedB precompute_native(
    std::span<const BoundaryG2> b,
    std::size_t n)
{
    return blst::detail::precompute_middle_product_b_impl<blst::MiddleProductScalarContext>(
        std::span<const BoundaryG2>(b),
        n,
        nullptr,
        blst::MiddleProductPairingMode::Direct);
}

inline void run_native(
    std::span<const BoundaryG1> a,
    const NativePrecomputedB& precomputed_b,
    std::span<PublicGt> out)
{
    NativeWorkspace workspace;
    blst::detail::middle_product_impl<blst::MiddleProductScalarContext>(
        std::span<const BoundaryG1>(a),
        precomputed_b,
        out,
        workspace,
        nullptr);
}
#endif

} // namespace

struct blst_middle_product_precomputed_b
{
    std::size_t n;
    NativePrecomputedB value;

    explicit blst_middle_product_precomputed_b(
        std::size_t input_n,
        NativePrecomputedB&& input_value)
        : n(input_n), value(std::move(input_value))
    {
    }
};

extern "C" {

blst_middle_product_precomputed_b* blst_middle_product_precompute_b(
    const blst_p2* b,
    size_t n)
{
    if (b == nullptr || n == 0) {
        return nullptr;
    }

    try {
        auto precomputed = precompute_native(
            std::span<const BoundaryG2>(reinterpret_cast<const BoundaryG2*>(b), 2 * n),
            n);
        return new blst_middle_product_precomputed_b(n, std::move(precomputed));
    } catch (...) {
        return nullptr;
    }
}

bool blst_middle_product(
    const blst_p1* a,
    size_t n,
    const blst_middle_product_precomputed_b* precomputed_b,
    blst_fp12* out)
{
    if (a == nullptr || out == nullptr || precomputed_b == nullptr || n == 0) {
        return false;
    }
    if (precomputed_b->n != n) {
        return false;
    }

    try {
        std::vector<PublicGt> native_out(n, PublicGt::identity());
        run_native(
            std::span<const BoundaryG1>(reinterpret_cast<const BoundaryG1*>(a), n),
            precomputed_b->value,
            std::span<PublicGt>(native_out));

        for (std::size_t i = 0; i < n; ++i) {
            reinterpret_cast<BoundaryFp12*>(out)[i] = gt_to_boundary(native_out[i]);
        }
        return true;
    } catch (...) {
        return false;
    }
}

void blst_middle_product_precomputed_b_free(
    blst_middle_product_precomputed_b* precomputed_b)
{
    delete precomputed_b;
}

}

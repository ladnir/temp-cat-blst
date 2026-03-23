#pragma once

#include "blst/avx/fields/conversion_inversion.h"
#include "blst/avx/fields/final_exponentiation.h"
#include "blst/avx/fields/fp12.h"
#include "blst/avx/gt/glv4.h"
#include "blst/avx/rns/bounded_ring.h"
#include "blst/scalar/curve/GT.h"

#include <array>
#include <cstddef>

namespace blst::avx {

namespace avx_gt_detail {

using Ring = BoundedRing<>;
using FP12 = FP12Ring<Ring>;

inline const Ring& ring()
{
    static const Ring value(
        BigInt("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787",
               10));
    return value;
}

inline const FP12& fp12ring()
{
    static const FP12 value(ring());
    return value;
}

inline std::array<uint64_t, 4> scalar_to_u64_le(const byte* scalar, size_t nbits)
{
    BLST_ASSERT(nbits <= 256);
    std::array<uint64_t, 4> out = {};
    const size_t nbytes = (nbits + 7) / 8;
    for (size_t i = 0; i < nbytes; ++i) {
        out[i / 8] |= static_cast<uint64_t>(scalar[i]) << ((i % 8) * 8);
    }
    return out;
}

inline blst::Fp12 to_scalar(const typename FP12::StandardElement& value)
{
    return avx_blst_direct_bridge::ring_fp12_to_scalar(ring(), value);
}

inline typename FP12::StandardElement from_scalar(const blst::Fp12& value)
{
    return avx_blst_direct_bridge::scalar_fp12_to_ring(ring(), value);
}

template<class Element>
inline typename FP12::StandardElement to_standard(const Element& value)
{
    return fp12ring().fp12_flat_mul(fp12ring().one(ring()), value, ring());
}

}  // namespace avx_gt_detail

struct GT
{
    using Field_type = avx_gt_detail::FP12::StandardElement;

    Field_type data = {};

    constexpr GT() = default;
    constexpr GT(const Field_type& value) : data(value) {}

    auto operator[](size_t idx) -> decltype(data[idx]) { return data[idx]; }
    auto operator[](size_t idx) const -> decltype(data[idx]) { return data[idx]; }

    static void add(GT& out, const GT& a, const GT& b)
    {
        out.data = avx_gt_detail::fp12ring().fp12_flat_mul(a.data, b.data, avx_gt_detail::ring());
    }

    static void dadd(GT& out, const GT& a, const GT& b)
    {
        add(out, a, b);
    }

    static void sub(GT& out, const GT& a, const GT& b)
    {
        const auto b_inv = avx_gt_detail::fp12ring().fp12_flat_conjugate(b.data, avx_gt_detail::ring());
        out.data = avx_gt_detail::fp12ring().fp12_flat_mul(a.data, b_inv, avx_gt_detail::ring());
    }

    static void doubling(GT& out, const GT& a)
    {
        out.data = avx_gt_detail::to_standard(
            avx_gt_detail::fp12ring().fp12_flat_cyclotomic_sqr(a.data, avx_gt_detail::ring()));
    }

    static void add_n_dbl(GT& out, const GT& a, size_t n)
    {
        add(out, out, a);
        for (size_t i = 0; i < n; ++i) {
            doubling(out, out);
        }
    }

    static void mult(GT& out, const GT& a, const byte* scalar, size_t nbits)
    {
        const auto words = avx_gt_detail::scalar_to_u64_le(scalar, nbits);
        if (nbits == 256) {
            const auto props = avx_glv4_compile(words.data());
            out.data = avx_cyclotomic_exp_glv4(a.data, props, avx_gt_detail::fp12ring(), avx_gt_detail::ring());
            return;
        }
        out.data = avx_cyclotomic_exp_basic(a.data, words.data(), avx_gt_detail::fp12ring(), avx_gt_detail::ring());
    }

    static void cneg(GT& p, limb_t cbit)
    {
        if (cbit) {
            p.data = avx_gt_detail::fp12ring().fp12_flat_conjugate(p.data, avx_gt_detail::ring());
        }
    }

    static void neg(GT& out, const GT& p)
    {
        out.data = avx_gt_detail::fp12ring().fp12_flat_conjugate(p.data, avx_gt_detail::ring());
    }

    static void cswap(GT& a, GT& b, bool_t cbit)
    {
        vec_cswap(&a, &b, sizeof(GT), cbit);
    }

    static void ccopy(GT& a, const GT& b, bool_t cbit)
    {
        vec_select(&a, &b, &a, sizeof(GT), cbit);
    }

    static limb_t equal(const GT& a, const GT& b)
    {
        return blst::scalar::GT::equal(
            blst::scalar::GT(avx_gt_detail::to_scalar(a.data)),
            blst::scalar::GT(avx_gt_detail::to_scalar(b.data)));
    }

    static const GT& generator()
    {
        static const GT g(identity());
        return g;
    }

    static GT identity()
    {
        return GT(avx_gt_detail::fp12ring().one(avx_gt_detail::ring()));
    }

    static int is_identity(const GT& p)
    {
        return equal(p, identity());
    }

    static bool_t in_group(const GT& p)
    {
        return blst::scalar::GT::in_group(blst::scalar::GT(avx_gt_detail::to_scalar(p.data)));
    }

    static bool_t is_cyclotomic(const GT& p)
    {
        return blst::scalar::GT::is_cyclotomic(blst::scalar::GT(avx_gt_detail::to_scalar(p.data)));
    }

    static void serialize(unsigned char serialized[12 * 48], const GT& p)
    {
        blst::scalar::GT::serialize(serialized, blst::scalar::GT(avx_gt_detail::to_scalar(p.data)));
    }

    static void deserialize(GT& out, const unsigned char serialized[12 * 48])
    {
        blst::scalar::GT scalar = {};
        blst::scalar::GT::deserialize(scalar, serialized);
        out.data = avx_gt_detail::from_scalar(scalar.data);
    }

    static void final_mult(GT& out, const GT& f)
    {
        FP12RNS_BLSTInverter<avx_gt_detail::Ring> inverter;
        out.data = final_exp(f.data, inverter, avx_gt_detail::fp12ring(), avx_gt_detail::ring());
    }

    static void mult_by_z(GT& out, const GT& a)
    {
        out.data = avx_gt_detail::to_standard(
            raise_to_z(a.data, avx_gt_detail::fp12ring(), avx_gt_detail::ring()));
    }

    static void mult_by_z_div_by_2(GT& out, const GT& a)
    {
        out.data = avx_gt_detail::to_standard(
            raise_to_z_div_by_2(a.data, avx_gt_detail::fp12ring(), avx_gt_detail::ring()));
    }

    static void frobenius_map(GT& out, const GT& a, size_t n)
    {
        switch (n) {
        case 1:
            out.data = avx_gt_detail::fp12ring().fp12_flat_frobenius_map1(a.data, avx_gt_detail::ring());
            break;
        case 2:
            out.data = avx_gt_detail::fp12ring().fp12_flat_frobenius_map2(a.data, avx_gt_detail::ring());
            break;
        case 3:
            out.data = avx_gt_detail::fp12ring().fp12_flat_frobenius_map3(a.data, avx_gt_detail::ring());
            break;
        default:
            BLST_ASSERT(false && "Unsupported Frobenius power");
        }
    }

    static void conjugate(GT& a)
    {
        a.data = avx_gt_detail::fp12ring().fp12_flat_conjugate(a.data, avx_gt_detail::ring());
    }
};

}  // namespace blst::avx

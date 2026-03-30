#pragma once

#include "fp12.h"
#include "../cpu/vector/vector_impl.hpp"

#include "blst/fields/Fp.h"
#include "blst/fields/Fp12.h"

#include <array>
#include <cstddef>
#include <cstring>

namespace blst::avx {

namespace avx_blst_scalar_bridge {

inline blst::Fp bigint_to_scalar_fp(const BigInt& value)
{
    blst::Fp standard = {};
    size_t written = 0;
    mpz_export(standard.data, &written, -1, sizeof(standard.data[0]), 0, 0, value.get_mpz().get_mpz_t());

    blst::Fp mont = {};
    blst::Fp::to_mont(mont, standard);
    return mont;
}

inline BigInt scalar_fp_to_bigint(const blst::Fp& value)
{
    blst::Fp standard = {};
    blst::Fp::from_mont(standard, value);

    mpz_class out;
    mpz_import(out.get_mpz_t(), blst::Fp::num_limb(), -1, sizeof(standard.data[0]), 0, 0, standard.data);
    return BigInt(out);
}

template<class Ring>
inline blst::Fp ring_to_scalar_fp(const Ring& ring, const typename Ring::StandardElement& value)
{
    return bigint_to_scalar_fp(ring.to_bigint(value));
}

template<class Ring>
inline typename Ring::StandardElement scalar_fp_to_ring(const Ring& ring, const blst::Fp& value)
{
    return ring.from_bigint(scalar_fp_to_bigint(value));
}

template<class Ring>
inline blst::Fp12 ring_fp12_to_scalar(
    const Ring& ring,
    const std::array<typename Ring::StandardElement, 12>& value)
{
    blst::Fp12 out = {};
    out[0][0][0] = ring_to_scalar_fp(ring, value[0]);
    out[0][0][1] = ring_to_scalar_fp(ring, value[1]);
    out[0][1][0] = ring_to_scalar_fp(ring, value[2]);
    out[0][1][1] = ring_to_scalar_fp(ring, value[3]);
    out[0][2][0] = ring_to_scalar_fp(ring, value[4]);
    out[0][2][1] = ring_to_scalar_fp(ring, value[5]);
    out[1][0][0] = ring_to_scalar_fp(ring, value[6]);
    out[1][0][1] = ring_to_scalar_fp(ring, value[7]);
    out[1][1][0] = ring_to_scalar_fp(ring, value[8]);
    out[1][1][1] = ring_to_scalar_fp(ring, value[9]);
    out[1][2][0] = ring_to_scalar_fp(ring, value[10]);
    out[1][2][1] = ring_to_scalar_fp(ring, value[11]);
    return out;
}

template<class Ring>
inline std::array<typename Ring::StandardElement, 12> scalar_fp12_to_ring(
    const Ring& ring,
    const blst::Fp12& value)
{
    return {
        scalar_fp_to_ring(ring, value[0][0][0]),
        scalar_fp_to_ring(ring, value[0][0][1]),
        scalar_fp_to_ring(ring, value[0][1][0]),
        scalar_fp_to_ring(ring, value[0][1][1]),
        scalar_fp_to_ring(ring, value[0][2][0]),
        scalar_fp_to_ring(ring, value[0][2][1]),
        scalar_fp_to_ring(ring, value[1][0][0]),
        scalar_fp_to_ring(ring, value[1][0][1]),
        scalar_fp_to_ring(ring, value[1][1][0]),
        scalar_fp_to_ring(ring, value[1][1][1]),
        scalar_fp_to_ring(ring, value[1][2][0]),
        scalar_fp_to_ring(ring, value[1][2][1]),
    };
}

}  // namespace avx_blst_scalar_bridge

namespace avx_blst_direct_bridge {

template<int limbs>
AVXVector<limbs> to_vec(uint8_t* le_buf)
{
    static_assert(limbs == 8, "to_vec currently only supports 8 limbs");

    uint64_t* limbs_64 = reinterpret_cast<uint64_t*>(le_buf);
    std::array<uint64_t, 8> digits = {};

    digits[0] = limbs_64[0] & ((1ULL << 52) - 1);
    digits[1] = ((limbs_64[0] >> 52) | (limbs_64[1] << 12)) & ((1ULL << 52) - 1);
    digits[2] = ((limbs_64[1] >> 40) | (limbs_64[2] << 24)) & ((1ULL << 52) - 1);
    digits[3] = ((limbs_64[2] >> 28) | (limbs_64[3] << 36)) & ((1ULL << 52) - 1);
    digits[4] = ((limbs_64[3] >> 16) | (limbs_64[4] << 48)) & ((1ULL << 52) - 1);
    digits[5] = (limbs_64[4] >> 4) & ((1ULL << 52) - 1);
    digits[6] = ((limbs_64[4] >> 56) | (limbs_64[5] << 8)) & ((1ULL << 52) - 1);
    digits[7] = (limbs_64[5] >> 44) & ((1ULL << 52) - 1);

    return AVXVector<limbs>(digits);
}

template<bool montgomery = false>
void from_vec(AVXVector<8>& vec_lo, AVXVector<8>& vec_hi, blst::Fp& ret)
{
    std::array<uint64_t, 8> digits_lo = {};
    std::array<uint64_t, 8> digits_hi = {};
    vec_hi.store(digits_hi.data());
    vec_lo.store(digits_lo.data());

    std::array<uint64_t, 9> digits = {};
    digits[0] = digits_lo[0];
    for (int i = 1; i < 8; ++i) {
        digits[i] = digits_hi[i - 1] + digits_lo[i];
    }
    digits[8] = digits_hi[7];

    const blst::Fp& kRR = blst::get_BLST_RR();
    blst::UFp n_data_even_64 = {};
    blst::UFp n_data_odd_64 = {};
    uint8_t* n_data_even = reinterpret_cast<uint8_t*>(n_data_even_64.data);
    uint8_t* n_data_odd = reinterpret_cast<uint8_t*>(n_data_odd_64.data);

    for (int i = 0; i < 9; i += 2) {
        *reinterpret_cast<uint64_t*>(n_data_even + (i / 2) * 13) = digits[i];
    }
    for (int i = 1; i < 8; i += 2) {
        *reinterpret_cast<uint64_t*>(n_data_odd + (i / 2) * 13 + 6) = digits[i] << 4;
    }

    blst::UFp sum = {};
    uint64_t carry = 0;
    for (int i = 0; i < 9; ++i) {
        uint64_t a = n_data_even_64.data[i];
        uint64_t b = n_data_odd_64.data[i];
        uint64_t s = a + b + carry;
        sum.data[i] = s;
        carry = (s < a) || (carry && s == a) ? 1 : 0;
    }

    blst::Fp temp = {};
    blst::Fp::redc(temp, sum);
    blst::Fp::mul(ret, temp, kRR);

    if constexpr (montgomery) {
        blst::Fp mont = {};
        blst::Fp::to_mont(mont, ret);
        ret = mont;
    }
}

template<class Ring, bool montgomery = false>
void convert_to_int(const Ring& ring, typename Ring::StandardElement& a, blst::Fp& ret)
{
    AVXVector<8> out_hi = AVXVector<8>(0);
    AVXVector<8> out_lo = AVXVector<8>(0);
    ring.sys.convert_from.conversion_matrix.rns_reduce_raw(a.m2.data, out_hi, out_lo);
    from_vec<montgomery>(out_lo, out_hi, ret);
}

template<class Ring, bool montgomery = false>
auto convert_to_rns(const Ring& ring, const blst::Fp& value)
{
    blst::Fp canonical = {};
    if constexpr (montgomery) {
        blst::Fp::from_mont(canonical, value);
    } else {
        canonical = value;
    }

    AVXVector<8> l = to_vec<8>(reinterpret_cast<uint8_t*>(canonical.data));
    return ring.template convert_to_rns_batch<1>(std::array<AVXVector<8>, 1>{l})[0];
}

template<class Ring>
blst::Fp ring_to_scalar_fp(const Ring& ring, const typename Ring::StandardElement& value)
{
    typename Ring::StandardElement mutable_value = value;
    blst::Fp mont = {};
    convert_to_int<Ring, true>(ring, mutable_value, mont);
    return mont;
}

template<class Ring>
typename Ring::StandardElement scalar_fp_to_ring(const Ring& ring, const blst::Fp& value)
{
    return convert_to_rns<Ring, true>(ring, value);
}

template<class Ring>
blst::Fp12 ring_fp12_to_scalar(
    const Ring& ring,
    const std::array<typename Ring::StandardElement, 12>& value)
{
    blst::Fp12 out = {};
    out[0][0][0] = ring_to_scalar_fp(ring, value[0]);
    out[0][0][1] = ring_to_scalar_fp(ring, value[1]);
    out[0][1][0] = ring_to_scalar_fp(ring, value[2]);
    out[0][1][1] = ring_to_scalar_fp(ring, value[3]);
    out[0][2][0] = ring_to_scalar_fp(ring, value[4]);
    out[0][2][1] = ring_to_scalar_fp(ring, value[5]);
    out[1][0][0] = ring_to_scalar_fp(ring, value[6]);
    out[1][0][1] = ring_to_scalar_fp(ring, value[7]);
    out[1][1][0] = ring_to_scalar_fp(ring, value[8]);
    out[1][1][1] = ring_to_scalar_fp(ring, value[9]);
    out[1][2][0] = ring_to_scalar_fp(ring, value[10]);
    out[1][2][1] = ring_to_scalar_fp(ring, value[11]);
    return out;
}

template<class Ring>
std::array<typename Ring::StandardElement, 12> scalar_fp12_to_ring(
    const Ring& ring,
    const blst::Fp12& value)
{
    return {
        scalar_fp_to_ring(ring, value[0][0][0]),
        scalar_fp_to_ring(ring, value[0][0][1]),
        scalar_fp_to_ring(ring, value[0][1][0]),
        scalar_fp_to_ring(ring, value[0][1][1]),
        scalar_fp_to_ring(ring, value[0][2][0]),
        scalar_fp_to_ring(ring, value[0][2][1]),
        scalar_fp_to_ring(ring, value[1][0][0]),
        scalar_fp_to_ring(ring, value[1][0][1]),
        scalar_fp_to_ring(ring, value[1][1][0]),
        scalar_fp_to_ring(ring, value[1][1][1]),
        scalar_fp_to_ring(ring, value[1][2][0]),
        scalar_fp_to_ring(ring, value[1][2][1]),
    };
}

}  // namespace avx_blst_direct_bridge

template<class Ring>
auto invert_via_blst(const Ring& ring, typename Ring::StandardElement& a)
{
    blst::Fp scalar = avx_blst_direct_bridge::ring_to_scalar_fp(ring, a);
    blst::Fp inverse = {};
    blst::Fp::reciprocal(inverse, scalar);
    return avx_blst_direct_bridge::scalar_fp_to_ring(ring, inverse);
}

template<class Ring>
class BLSTInverter {
public:
    using StandardElement = typename Ring::StandardElement;

    StandardElement invert(const StandardElement& a, const Ring& ring) const
    {
        StandardElement mutable_a = a;
        return invert_via_blst(ring, mutable_a);
    }
};

template<class Ring>
std::array<typename Ring::StandardElement, 12> invert_fp12_via_blst(
    const Ring& ring,
    const std::array<typename Ring::StandardElement, 12>& fp12_ring)
{
    blst::Fp12 scalar = avx_blst_direct_bridge::ring_fp12_to_scalar(ring, fp12_ring);
    blst::Fp12 inverse = {};
    blst::Fp12::inverse(inverse, scalar);
    return avx_blst_direct_bridge::scalar_fp12_to_ring(ring, inverse);
}

template<class Ring>
std::array<typename Ring::StandardElement, 12> invert_fp12_via_blst_bigint(
    const Ring& ring,
    const std::array<typename Ring::StandardElement, 12>& fp12_ring)
{
    return invert_fp12_via_blst(ring, fp12_ring);
}

template<class Ring>
class FP12BLSTInverter {
public:
    using StandardElement = std::array<typename Ring::StandardElement, 12>;
    using FP12R = FP12Ring<Ring>;

    StandardElement invert(const StandardElement& f, const FP12R&, const Ring& ring) const
    {
        return invert_fp12_via_blst(ring, f);
    }
};

template<class Ring>
class FP12RNS_BLSTInverter {
public:
    using FP12R = FP12Ring<Ring>;
    using StandardElement = typename FP12R::StandardElement;
    using Inverter = BLSTInverter<Ring>;

    Inverter inverter;

    StandardElement invert(const StandardElement& f, const FP12R& fp12ring, const Ring& ring) const
    {
        return fp12ring.template fp12_flat_inverse<Inverter>(f, inverter, ring);
    }
};

}  // namespace blst::avx

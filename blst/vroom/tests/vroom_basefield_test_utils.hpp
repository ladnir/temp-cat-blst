#pragma once

#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fp12.h"
#include "blst/vroom/cpu/precompute/gmp_wrapper.hpp"
#include "blst/vroom/src/bounded_ring.hpp"
#include "blst/vroom/src/fp2.hpp"
#include "blst/vroom/src/fp12.hpp"

#include <array>
#include <cstdint>
#include <utility>

namespace blst_vroom_test {

using VroomRing = BoundedRing<>;
using VroomFp2Ring = FP2Ring<VroomRing>;
using VroomFp2 = VroomFp2Ring::StandardElement;
using VroomFp12Ring = FP12Ring<VroomRing>;
using VroomFp12 = VroomFp12Ring::StandardElement;

inline constexpr char kBls12_381BaseModulus[] =
    "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787";

inline BigInt modulus_bigint()
{
    return BigInt(kBls12_381BaseModulus, 10);
}

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

inline blst::Fp2 bigint_pair_to_scalar_fp2(const BigInt& x, const BigInt& y)
{
    blst::Fp2 out = {};
    out[0] = bigint_to_scalar_fp(x);
    out[1] = bigint_to_scalar_fp(y);
    return out;
}

inline std::pair<BigInt, BigInt> scalar_fp2_to_bigints(const blst::Fp2& value)
{
    blst::Fp2 standard = {};
    blst::Fp::from_mont(standard[0], value[0]);
    blst::Fp::from_mont(standard[1], value[1]);

    mpz_class x;
    mpz_import(x.get_mpz_t(), blst::Fp::num_limb(), -1, sizeof(standard[0].data[0]), 0, 0, standard[0].data);

    mpz_class y;
    mpz_import(y.get_mpz_t(), blst::Fp::num_limb(), -1, sizeof(standard[1].data[0]), 0, 0, standard[1].data);

    return {BigInt(x), BigInt(y)};
}

inline VroomFp2 make_vroom_fp2(const VroomRing& ring, const BigInt& x, const BigInt& y)
{
    return VroomFp2(ring.from_bigint(x), ring.from_bigint(y));
}

inline std::array<BigInt, 12> scalar_fp12_to_bigints(const blst::Fp12& value)
{
    return {
        scalar_fp_to_bigint(value[0][0][0]),
        scalar_fp_to_bigint(value[0][0][1]),
        scalar_fp_to_bigint(value[0][1][0]),
        scalar_fp_to_bigint(value[0][1][1]),
        scalar_fp_to_bigint(value[0][2][0]),
        scalar_fp_to_bigint(value[0][2][1]),
        scalar_fp_to_bigint(value[1][0][0]),
        scalar_fp_to_bigint(value[1][0][1]),
        scalar_fp_to_bigint(value[1][1][0]),
        scalar_fp_to_bigint(value[1][1][1]),
        scalar_fp_to_bigint(value[1][2][0]),
        scalar_fp_to_bigint(value[1][2][1]),
    };
}

inline blst::Fp12 bigint_array_to_scalar_fp12(const std::array<BigInt, 12>& coeffs)
{
    blst::Fp12 out = {};
    out[0][0][0] = bigint_to_scalar_fp(coeffs[0]);
    out[0][0][1] = bigint_to_scalar_fp(coeffs[1]);
    out[0][1][0] = bigint_to_scalar_fp(coeffs[2]);
    out[0][1][1] = bigint_to_scalar_fp(coeffs[3]);
    out[0][2][0] = bigint_to_scalar_fp(coeffs[4]);
    out[0][2][1] = bigint_to_scalar_fp(coeffs[5]);
    out[1][0][0] = bigint_to_scalar_fp(coeffs[6]);
    out[1][0][1] = bigint_to_scalar_fp(coeffs[7]);
    out[1][1][0] = bigint_to_scalar_fp(coeffs[8]);
    out[1][1][1] = bigint_to_scalar_fp(coeffs[9]);
    out[1][2][0] = bigint_to_scalar_fp(coeffs[10]);
    out[1][2][1] = bigint_to_scalar_fp(coeffs[11]);
    return out;
}

inline std::array<BigInt, 12> seeded_fp12_coeffs(uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);

    std::array<BigInt, 12> coeffs;
    for (BigInt& coeff : coeffs) {
        coeff = BigInt(rng.get_z_range(modulus.get_mpz()));
    }
    return coeffs;
}

inline VroomFp12 make_vroom_fp12(const VroomRing& ring, const std::array<BigInt, 12>& coeffs)
{
    VroomFp12 out = {};
    for (size_t i = 0; i < coeffs.size(); ++i) {
        out[i] = ring.from_bigint(coeffs[i]);
    }
    return out;
}

inline VroomFp12 make_vroom_fp12_from_seed(const VroomRing& ring, uint64_t seed)
{
    return make_vroom_fp12(ring, seeded_fp12_coeffs(seed));
}

inline blst::Fp12 make_scalar_fp12_from_seed(uint64_t seed)
{
    return bigint_array_to_scalar_fp12(seeded_fp12_coeffs(seed));
}

inline blst::Fp12 make_scalar_cyclotomic_fp12_from_seed(uint64_t seed)
{
    blst::Fp12 out = make_scalar_fp12_from_seed(seed);
    blst::Fp12::final_exp(out, out);
    return out;
}

inline VroomFp12 make_vroom_cyclotomic_fp12_from_seed(const VroomRing& ring, uint64_t seed)
{
    return make_vroom_fp12(ring, scalar_fp12_to_bigints(make_scalar_cyclotomic_fp12_from_seed(seed)));
}

inline std::array<BigInt, 6> seeded_sparse_fp6_coeffs(uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);

    std::array<BigInt, 6> coeffs;
    for (BigInt& coeff : coeffs) {
        coeff = BigInt(rng.get_z_range(modulus.get_mpz()));
    }
    return coeffs;
}

inline blst::Fp6 bigint_array_to_scalar_fp6(const std::array<BigInt, 6>& coeffs)
{
    blst::Fp6 out = {};
    out[0] = bigint_pair_to_scalar_fp2(coeffs[0], coeffs[1]);
    out[1] = bigint_pair_to_scalar_fp2(coeffs[2], coeffs[3]);
    out[2] = bigint_pair_to_scalar_fp2(coeffs[4], coeffs[5]);
    return out;
}

inline blst::Fp6 make_scalar_sparse_fp6_from_seed(uint64_t seed)
{
    return bigint_array_to_scalar_fp6(seeded_sparse_fp6_coeffs(seed));
}

inline std::array<VroomRing::StandardElement, 6> make_vroom_sparse_fp6_coeffs(const VroomRing& ring, const std::array<BigInt, 6>& coeffs)
{
    std::array<VroomRing::StandardElement, 6> out = {};
    for (size_t i = 0; i < out.size(); ++i) {
        out[i] = ring.from_bigint(coeffs[i]);
    }
    return out;
}

inline std::array<VroomRing::StandardElement, 6> make_vroom_sparse_fp6_from_seed(const VroomRing& ring, uint64_t seed)
{
    return make_vroom_sparse_fp6_coeffs(ring, seeded_sparse_fp6_coeffs(seed));
}

inline VroomFp2 make_vroom_fp2_from_sparse6(const std::array<VroomRing::StandardElement, 6>& coeffs, size_t pair_idx)
{
    return VroomFp2(coeffs[pair_idx * 2], coeffs[pair_idx * 2 + 1]);
}

template<class VroomFp12Like>
inline std::array<BigInt, 12> vroom_fp12_to_bigints(const VroomRing& ring, const VroomFp12Like& value)
{
    std::array<BigInt, 12> out;
    for (size_t i = 0; i < out.size(); ++i) {
        typename VroomRing::StandardElement canonical = {};
        canonical.m2 = value[i].m2;
        out[i] = ring.to_bigint(canonical);
    }
    return out;
}

}  // namespace blst_vroom_test

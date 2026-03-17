#pragma once

#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fp12.h"
#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/vroom/cpu/precompute/gmp_wrapper.hpp"
#include "blst/vroom/src/ec.hpp"
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
using VroomAffineFp = AffinePoint<VroomRing::StandardElement>;
using VroomAffineFp2 = AffinePoint<VroomFp2>;
using VroomProjectiveFp2 = ProjectivePoint<VroomFp2>;

inline constexpr char kBls12_381BaseModulus[] =
    "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787";
inline constexpr char kG1GeneratorX[] =
    "3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507";
inline constexpr char kG1GeneratorY[] =
    "1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569";
inline constexpr char kG2GeneratorXC0[] =
    "352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160";
inline constexpr char kG2GeneratorXC1[] =
    "3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758";
inline constexpr char kG2GeneratorYC0[] =
    "1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905";
inline constexpr char kG2GeneratorYC1[] =
    "927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582";

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

inline blst::G1AffinePoint make_scalar_g1_generator_affine()
{
    blst::G1AffinePoint affine = {};
    affine.X = bigint_to_scalar_fp(BigInt(kG1GeneratorX, 10));
    affine.Y = bigint_to_scalar_fp(BigInt(kG1GeneratorY, 10));
    return affine;
}

inline blst::G2AffinePoint make_scalar_g2_generator_affine()
{
    blst::G2AffinePoint affine = {};
    affine.X = bigint_pair_to_scalar_fp2(BigInt(kG2GeneratorXC0, 10), BigInt(kG2GeneratorXC1, 10));
    affine.Y = bigint_pair_to_scalar_fp2(BigInt(kG2GeneratorYC0, 10), BigInt(kG2GeneratorYC1, 10));
    return affine;
}

inline blst::G1AffinePoint make_scalar_g1_doubled_affine()
{
    const blst::G1Point generator = blst::G1Point::generator();
    blst::G1Point doubled = {};
    blst::G1Point::doubling(doubled, generator);

    blst::G1AffinePoint affine = {};
    blst::G1Point::to_affine(affine, doubled);
    return affine;
}

inline blst::G2AffinePoint make_scalar_g2_doubled_affine()
{
    const blst::G2Point generator = blst::G2Point::generator();
    blst::G2Point doubled = {};
    blst::G2Point::doubling(doubled, generator);

    blst::G2AffinePoint affine = {};
    blst::G2Point::to_affine(affine, doubled);
    return affine;
}

inline VroomAffineFp make_vroom_g1_affine(const VroomRing& ring, const blst::G1AffinePoint& point)
{
    return {
        ring.from_bigint(scalar_fp_to_bigint(point.X)),
        ring.from_bigint(scalar_fp_to_bigint(point.Y)),
    };
}

inline VroomAffineFp2 make_vroom_g2_affine(const VroomRing& ring, const blst::G2AffinePoint& point)
{
    const auto x = scalar_fp2_to_bigints(point.X);
    const auto y = scalar_fp2_to_bigints(point.Y);
    return {
        make_vroom_fp2(ring, x.first, x.second),
        make_vroom_fp2(ring, y.first, y.second),
    };
}

inline VroomAffineFp make_vroom_g1_generator_affine(const VroomRing& ring)
{
    return make_vroom_g1_affine(ring, make_scalar_g1_generator_affine());
}

inline VroomAffineFp2 make_vroom_g2_generator_affine(const VroomRing& ring)
{
    return make_vroom_g2_affine(ring, make_scalar_g2_generator_affine());
}

inline blst::G1AffinePoint make_scalar_px2(const blst::G1AffinePoint& point)
{
    blst::G1AffinePoint px2 = {};
    blst::Fp::add(px2.X, point.X, point.X);
    blst::Fp::neg(px2.X, px2.X);
    blst::Fp::add(px2.Y, point.Y, point.Y);
    return px2;
}

inline blst::Fp12 make_scalar_miller_first_step(const blst::G1AffinePoint& p, const blst::G2AffinePoint& q)
{
    blst::Fp6 qlines[68] = {};
    blst::GT::precompute_lines(qlines, q);

    const blst::G1AffinePoint px2 = make_scalar_px2(p);
    blst::Fp6 line = {};
    line[0] = qlines[0][0];
    line[1] = qlines[0][1];
    line[2] = qlines[0][2];

    blst::Fp::mul(line[1][0], line[1][0], px2.X);
    blst::Fp::mul(line[1][1], line[1][1], px2.X);
    blst::Fp::mul(line[2][0], line[2][0], px2.Y);
    blst::Fp::mul(line[2][1], line[2][1], px2.Y);

    blst::Fp12 out = {};
    out[0][0] = line[0];
    out[0][1] = line[1];
    out[1][1] = line[2];
    return out;
}

inline blst::GT make_scalar_miller_loop(const blst::G1AffinePoint& p, const blst::G2AffinePoint& q)
{
    blst::GT result = {};
    blst::GT::miller_loop(result, &q, &p);
    return result;
}

inline blst::GT make_scalar_final_pairing(const blst::G1AffinePoint& p, const blst::G2AffinePoint& q)
{
    blst::GT result = make_scalar_miller_loop(p, q);
    blst::GT::final_mult(result, result);
    return result;
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

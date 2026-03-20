#pragma once

#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fp12.h"
#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/avx/cpu/precompute/gmp_wrapper.hpp"
#include "blst/avx/curve/ec.h"
#include "blst/avx/rns/bounded_ring.h"
#include "blst/avx/fp2/fp2.h"
#include "blst/avx/fp12/fp12.h"

#include <array>
#include <cstdint>
#include <utility>

namespace blst_avx_test {

using namespace blst::avx;

using AvxRing = BoundedRing<>;
using AvxFp2Ring = FP2Ring<AvxRing>;
using AvxFp2 = AvxFp2Ring::StandardElement;
using AvxFp12Ring = FP12Ring<AvxRing>;
using AvxFp12 = AvxFp12Ring::StandardElement;
using Avx50Ring = BoundedRing<381, 8, 50>;
using Avx50Fp2Ring = FP2Ring<Avx50Ring>;
using Avx50Fp2 = Avx50Fp2Ring::StandardElement;
using Avx50Fp12Ring = FP12Ring<Avx50Ring>;
using Avx50Fp12 = Avx50Fp12Ring::StandardElement;
using AvxAffineFp = AffinePoint<AvxRing::StandardElement>;
using AvxAffineFp2 = AffinePoint<AvxFp2>;
using AvxProjectiveFp2 = ProjectivePoint<AvxFp2>;
using AvxProjectiveFp = ProjectivePoint<AvxRing::StandardElement>;
using AvxG1 = G1<AvxRing>;
using AvxG2 = G2<AvxRing>;
using AvxCurveRing = BoundedRing<381, 8, 52, -1932, 2377, 12>;
using AvxCurveFp2Ring = FP2Ring<AvxCurveRing>;
using AvxCurveFp2 = AvxCurveFp2Ring::StandardElement;
using AvxCurveAffineFp = AffinePoint<AvxCurveRing::StandardElement>;
using AvxCurveAffineFp2 = AffinePoint<AvxCurveFp2>;
using AvxCurveProjectiveFp = ProjectivePoint<AvxCurveRing::StandardElement>;
using AvxCurveProjectiveFp2 = ProjectivePoint<AvxCurveFp2>;
using AvxCurveG1 = G1<AvxCurveRing>;
using AvxCurveG2 = G2<AvxCurveRing>;
using AvxCurve50Ring = BoundedRing<381, 8, 50, -1932, 2377, 12>;
using AvxCurve50Fp2Ring = FP2Ring<AvxCurve50Ring>;
using AvxCurve50Fp2 = AvxCurve50Fp2Ring::StandardElement;
using AvxCurve50AffineFp = AffinePoint<AvxCurve50Ring::StandardElement>;
using AvxCurve50AffineFp2 = AffinePoint<AvxCurve50Fp2>;
using AvxCurve50ProjectiveFp = ProjectivePoint<AvxCurve50Ring::StandardElement>;
using AvxCurve50ProjectiveFp2 = ProjectivePoint<AvxCurve50Fp2>;
using AvxCurve50G1 = G1<AvxCurve50Ring>;
using AvxCurve50G2 = G2<AvxCurve50Ring>;

inline constexpr char kBls12_381BaseModulus[] =
    "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787";
inline constexpr char kBls12_381ScalarModulusHex[] =
    "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";
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

inline BigInt scalar_modulus_bigint()
{
    return BigInt(kBls12_381ScalarModulusHex, 16);
}

template<size_t N>
inline std::array<uint8_t, N> bigint_to_bytes_le_array(const BigInt& value)
{
    std::array<uint8_t, N> bytes = {};
    BigInt current = value;
    for (size_t i = 0; i < N; ++i) {
        bytes[i] = static_cast<uint8_t>(current.to_ulong() & 0xff);
        current = current >> 8;
    }
    return bytes;
}

inline BigInt seeded_scalar_bigint(uint64_t seed)
{
    const BigInt modulus = scalar_modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return BigInt(rng.get_z_range(modulus.get_mpz()));
}

inline std::array<uint8_t, 32> seeded_scalar_bytes_le(uint64_t seed)
{
    return bigint_to_bytes_le_array<32>(seeded_scalar_bigint(seed));
}

inline std::array<uint64_t, 4> scalar_bytes_to_limbs_le(const std::array<uint8_t, 32>& bytes)
{
    std::array<uint64_t, 4> limbs = {};
    for (size_t limb_idx = 0; limb_idx < limbs.size(); ++limb_idx) {
        uint64_t limb = 0;
        for (size_t byte_idx = 0; byte_idx < sizeof(uint64_t); ++byte_idx) {
            limb |= static_cast<uint64_t>(bytes[limb_idx * sizeof(uint64_t) + byte_idx]) << (byte_idx * 8);
        }
        limbs[limb_idx] = limb;
    }
    return limbs;
}

inline std::array<uint64_t, 4> seeded_scalar_limbs_le(uint64_t seed)
{
    return scalar_bytes_to_limbs_le(seeded_scalar_bytes_le(seed));
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

template<class Ring>
inline typename FP2Ring<Ring>::StandardElement make_avx_fp2(const Ring& ring, const BigInt& x, const BigInt& y)
{
    return typename FP2Ring<Ring>::StandardElement(ring.from_bigint(x), ring.from_bigint(y));
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

inline blst::G1Point make_scalar_g1_projective_from_seed(uint64_t seed)
{
    const std::array<uint8_t, 32> scalar_bytes = seeded_scalar_bytes_le(seed);
    blst::G1Point point = {};
    blst::G1Point::mult(point, blst::G1Point::generator(), scalar_bytes.data(), 256);
    return point;
}

inline blst::G2Point make_scalar_g2_projective_from_seed(uint64_t seed)
{
    const std::array<uint8_t, 32> scalar_bytes = seeded_scalar_bytes_le(seed);
    blst::G2Point point = {};
    blst::G2Point::mult(point, blst::G2Point::generator(), scalar_bytes.data(), 256);
    return point;
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

template<class Ring>
inline AffinePoint<typename Ring::StandardElement> make_avx_g1_affine(const Ring& ring, const blst::G1AffinePoint& point)
{
    return {
        ring.from_bigint(scalar_fp_to_bigint(point.X)),
        ring.from_bigint(scalar_fp_to_bigint(point.Y)),
    };
}

template<class Ring>
inline AffinePoint<typename FP2Ring<Ring>::StandardElement> make_avx_g2_affine(const Ring& ring, const blst::G2AffinePoint& point)
{
    const auto x = scalar_fp2_to_bigints(point.X);
    const auto y = scalar_fp2_to_bigints(point.Y);
    return {
        typename FP2Ring<Ring>::StandardElement(ring.from_bigint(x.first), ring.from_bigint(x.second)),
        typename FP2Ring<Ring>::StandardElement(ring.from_bigint(y.first), ring.from_bigint(y.second)),
    };
}

template<class Ring>
inline AffinePoint<typename Ring::StandardElement> make_avx_g1_generator_affine(const Ring& ring)
{
    return make_avx_g1_affine(ring, make_scalar_g1_generator_affine());
}

template<class Ring>
inline AffinePoint<typename FP2Ring<Ring>::StandardElement> make_avx_g2_generator_affine(const Ring& ring)
{
    return make_avx_g2_affine(ring, make_scalar_g2_generator_affine());
}

inline AvxProjectiveFp make_avx_g1_projective(const AvxRing& ring, const AvxAffineFp& point)
{
    return {point.x, point.y, ring.one()};
}

inline AvxProjectiveFp2 make_avx_g2_projective(const AvxFp2Ring& fp2ring, const AvxAffineFp2& point)
{
    return {point.x, point.y, fp2ring.one()};
}

inline AvxCurveProjectiveFp make_avx_curve_g1_projective(const AvxCurveRing& ring, const AvxCurveAffineFp& point)
{
    return {point.x, point.y, ring.one()};
}

inline AvxCurveProjectiveFp2 make_avx_curve_g2_projective(const AvxCurveFp2Ring& fp2ring, const AvxCurveAffineFp2& point)
{
    return {point.x, point.y, fp2ring.one()};
}

inline AvxCurve50ProjectiveFp make_avx_curve50_g1_projective(const AvxCurve50Ring& ring, const AvxCurve50AffineFp& point)
{
    return {point.x, point.y, ring.one()};
}

inline AvxCurve50ProjectiveFp2 make_avx_curve50_g2_projective(const AvxCurve50Fp2Ring& fp2ring, const AvxCurve50AffineFp2& point)
{
    return {point.x, point.y, fp2ring.one()};
}

inline AvxProjectiveFp make_random_avx_g1_projective(const AvxRing& ring, uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return {
        ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
        ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
        ring.one(),
    };
}

inline AvxAffineFp make_random_avx_g1_affine(const AvxRing& ring, uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return {
        ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
        ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
    };
}

inline AvxProjectiveFp2 make_random_avx_g2_projective(const AvxFp2Ring& fp2ring, const AvxRing& ring, uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return {
        AvxFp2(
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())))),
        AvxFp2(
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())))),
        fp2ring.one(),
    };
}

inline AvxAffineFp2 make_random_avx_g2_affine(const AvxRing& ring, uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return {
        AvxFp2(
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())))),
        AvxFp2(
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())))),
    };
}

inline AvxCurveProjectiveFp make_random_avx_curve_g1_projective(const AvxCurveRing& ring, uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return {
        ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
        ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
        ring.one(),
    };
}

inline AvxCurveAffineFp make_random_avx_curve_g1_affine(const AvxCurveRing& ring, uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return {
        ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
        ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
    };
}

inline AvxCurveProjectiveFp2 make_random_avx_curve_g2_projective(const AvxCurveFp2Ring& fp2ring, const AvxCurveRing& ring, uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return {
        AvxCurveFp2(
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())))),
        AvxCurveFp2(
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())))),
        fp2ring.one(),
    };
}

inline AvxCurveAffineFp2 make_random_avx_curve_g2_affine(const AvxCurveRing& ring, uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return {
        AvxCurveFp2(
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())))),
        AvxCurveFp2(
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())))),
    };
}

inline AvxCurve50ProjectiveFp make_random_avx_curve50_g1_projective(const AvxCurve50Ring& ring, uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return {
        ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
        ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
        ring.one(),
    };
}

inline AvxCurve50AffineFp make_random_avx_curve50_g1_affine(const AvxCurve50Ring& ring, uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return {
        ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
        ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
    };
}

inline AvxCurve50ProjectiveFp2 make_random_avx_curve50_g2_projective(const AvxCurve50Fp2Ring& fp2ring, const AvxCurve50Ring& ring, uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return {
        AvxCurve50Fp2(
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())))),
        AvxCurve50Fp2(
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())))),
        fp2ring.one(),
    };
}

inline AvxCurve50AffineFp2 make_random_avx_curve50_g2_affine(const AvxCurve50Ring& ring, uint64_t seed)
{
    const BigInt modulus = modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return {
        AvxCurve50Fp2(
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())))),
        AvxCurve50Fp2(
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())))),
    };
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

template<class Ring>
inline std::pair<BigInt, BigInt> avx_g1_projective_to_affine_bigints(
    const Ring& ring,
    const ProjectivePoint<typename Ring::StandardElement>& point)
{
    const BigInt modulus = modulus_bigint();
    const BigInt x = ring.to_bigint(point.x);
    const BigInt y = ring.to_bigint(point.y);
    const BigInt z = ring.to_bigint(point.z);

    if (z == BigInt(0)) {
        return {BigInt(0), BigInt(0)};
    }

    const BigInt z_inv = z.mod_inverse(modulus);
    BigInt affine_x = (x * z_inv) % modulus;
    if (affine_x < 0) {
        affine_x = affine_x + modulus;
    }
    BigInt affine_y = (y * z_inv) % modulus;
    if (affine_y < 0) {
        affine_y = affine_y + modulus;
    }
    return {affine_x, affine_y};
}

template<class Ring>
inline std::pair<std::pair<BigInt, BigInt>, std::pair<BigInt, BigInt>> avx_g2_projective_to_affine_bigints(
    const Ring& ring,
    const ProjectivePoint<typename FP2Ring<Ring>::StandardElement>& point)
{
    const BigInt modulus = modulus_bigint();
    const BigInt x0 = ring.to_bigint(point.x.x());
    const BigInt x1 = ring.to_bigint(point.x.y());
    const BigInt y0 = ring.to_bigint(point.y.x());
    const BigInt y1 = ring.to_bigint(point.y.y());
    const BigInt z0 = ring.to_bigint(point.z.x());
    const BigInt z1 = ring.to_bigint(point.z.y());

    if (z0 == BigInt(0) && z1 == BigInt(0)) {
        return {{BigInt(0), BigInt(0)}, {BigInt(0), BigInt(0)}};
    }

    BigInt norm = (z0 * z0 + z1 * z1) % modulus;
    if (norm < 0) {
        norm = norm + modulus;
    }
    const BigInt norm_inv = norm.mod_inverse(modulus);

    BigInt z_inv0 = (z0 * norm_inv) % modulus;
    if (z_inv0 < 0) {
        z_inv0 = z_inv0 + modulus;
    }
    BigInt z_inv1 = (-z1 * norm_inv) % modulus;
    if (z_inv1 < 0) {
        z_inv1 = z_inv1 + modulus;
    }

    BigInt affine_x0 = (x0 * z_inv0 - x1 * z_inv1) % modulus;
    if (affine_x0 < 0) {
        affine_x0 = affine_x0 + modulus;
    }
    BigInt affine_x1 = (x0 * z_inv1 + x1 * z_inv0) % modulus;
    if (affine_x1 < 0) {
        affine_x1 = affine_x1 + modulus;
    }
    BigInt affine_y0 = (y0 * z_inv0 - y1 * z_inv1) % modulus;
    if (affine_y0 < 0) {
        affine_y0 = affine_y0 + modulus;
    }
    BigInt affine_y1 = (y0 * z_inv1 + y1 * z_inv0) % modulus;
    if (affine_y1 < 0) {
        affine_y1 = affine_y1 + modulus;
    }

    return {{affine_x0, affine_x1}, {affine_y0, affine_y1}};
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

template<class Ring>
inline typename FP12Ring<Ring>::StandardElement make_avx_fp12(const Ring& ring, const std::array<BigInt, 12>& coeffs)
{
    typename FP12Ring<Ring>::StandardElement out = {};
    for (size_t i = 0; i < coeffs.size(); ++i) {
        out[i] = ring.from_bigint(coeffs[i]);
    }
    return out;
}

template<class Ring>
inline typename FP12Ring<Ring>::StandardElement make_avx_fp12_from_seed(const Ring& ring, uint64_t seed)
{
    return make_avx_fp12(ring, seeded_fp12_coeffs(seed));
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

template<class Ring>
inline typename FP12Ring<Ring>::StandardElement make_avx_cyclotomic_fp12_from_seed(const Ring& ring, uint64_t seed)
{
    return make_avx_fp12(ring, scalar_fp12_to_bigints(make_scalar_cyclotomic_fp12_from_seed(seed)));
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

template<class Ring>
inline std::array<typename Ring::StandardElement, 6> make_avx_sparse_fp6_coeffs(const Ring& ring, const std::array<BigInt, 6>& coeffs)
{
    std::array<typename Ring::StandardElement, 6> out = {};
    for (size_t i = 0; i < out.size(); ++i) {
        out[i] = ring.from_bigint(coeffs[i]);
    }
    return out;
}

template<class Ring>
inline std::array<typename Ring::StandardElement, 6> make_avx_sparse_fp6_from_seed(const Ring& ring, uint64_t seed)
{
    return make_avx_sparse_fp6_coeffs(ring, seeded_sparse_fp6_coeffs(seed));
}

template<class Ring>
inline typename FP2Ring<Ring>::StandardElement make_avx_fp2_from_sparse6(const std::array<typename Ring::StandardElement, 6>& coeffs, size_t pair_idx)
{
    return typename FP2Ring<Ring>::StandardElement(coeffs[pair_idx * 2], coeffs[pair_idx * 2 + 1]);
}

template<class Ring, class AvxFp12Like>
inline std::array<BigInt, 12> avx_fp12_to_bigints(const Ring& ring, const AvxFp12Like& value)
{
    std::array<BigInt, 12> out;
    for (size_t i = 0; i < out.size(); ++i) {
        typename Ring::StandardElement canonical = {};
        canonical.m2 = value[i].m2;
        out[i] = ring.to_bigint(canonical);
    }
    return out;
}

}  // namespace blst_avx_test

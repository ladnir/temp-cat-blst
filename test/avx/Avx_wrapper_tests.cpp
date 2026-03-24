#define BLST_AVX_SKIP_TEST_VALUES

#include <gtest/gtest.h>

#include "blst/avx/curve/GT.h"
#include "blst/avx/curve/map_to_G.h"
#include "blst/avx/curve/point.h"
#include "blst/scalar/curve/map_to_G.h"
#include "avx_test_utils.h"

#include <array>
#include <cstring>

using namespace blst_avx_test;

namespace {

blst::avx::G1Point scalar_to_avx(const blst::scalar::G1Point& point)
{
    blst::scalar::G1AffinePoint scalar_affine = {};
    blst::scalar::G1Point::to_affine(scalar_affine, point);

    blst::avx::G1AffinePoint avx_affine = {
        blst::avx::avx_curve_detail::scalar_fp_to_ring(scalar_affine.X),
        blst::avx::avx_curve_detail::scalar_fp_to_ring(scalar_affine.Y),
    };

    blst::avx::G1Point out = {};
    blst::avx::G1Point::from_affine(out, avx_affine);
    return out;
}

blst::avx::G2Point scalar_to_avx(const blst::scalar::G2Point& point)
{
    blst::scalar::G2AffinePoint scalar_affine = {};
    blst::scalar::G2Point::to_affine(scalar_affine, point);

    blst::avx::G2AffinePoint avx_affine = {
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(scalar_affine.X),
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(scalar_affine.Y),
    };

    blst::avx::G2Point out = {};
    blst::avx::G2Point::from_affine(out, avx_affine);
    return out;
}

blst::avx::G1AffinePoint scalar_to_avx(const blst::G1AffinePoint& point)
{
    return {
        blst::avx::avx_curve_detail::scalar_fp_to_ring(point.X),
        blst::avx::avx_curve_detail::scalar_fp_to_ring(point.Y),
    };
}

blst::avx::G2AffinePoint scalar_to_avx(const blst::G2AffinePoint& point)
{
    return {
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(point.X),
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(point.Y),
    };
}

void expect_same_affine(const blst::scalar::G1Point& scalar_point, const blst::avx::G1Point& avx_point)
{
    blst::scalar::G1AffinePoint scalar_affine = {};
    blst::scalar::G1Point::to_affine(scalar_affine, scalar_point);

    blst::avx::G1AffinePoint avx_affine = {};
    blst::avx::G1Point::to_affine(avx_affine, avx_point);

    EXPECT_TRUE(blst::Fp::equal(scalar_affine.X, blst::avx::avx_curve_detail::ring_fp_to_scalar(avx_affine.X)));
    EXPECT_TRUE(blst::Fp::equal(scalar_affine.Y, blst::avx::avx_curve_detail::ring_fp_to_scalar(avx_affine.Y)));
}

void expect_same_affine(const blst::scalar::G2Point& scalar_point, const blst::avx::G2Point& avx_point)
{
    blst::scalar::G2AffinePoint scalar_affine = {};
    blst::scalar::G2Point::to_affine(scalar_affine, scalar_point);

    blst::avx::G2AffinePoint avx_affine = {};
    blst::avx::G2Point::to_affine(avx_affine, avx_point);

    const blst::Fp2 avx_x = blst::avx::avx_curve_detail::ring_fp2_to_scalar(avx_affine.X);
    const blst::Fp2 avx_y = blst::avx::avx_curve_detail::ring_fp2_to_scalar(avx_affine.Y);
    EXPECT_TRUE(blst::Fp::equal(scalar_affine.X[0], avx_x[0]));
    EXPECT_TRUE(blst::Fp::equal(scalar_affine.X[1], avx_x[1]));
    EXPECT_TRUE(blst::Fp::equal(scalar_affine.Y[0], avx_y[0]));
    EXPECT_TRUE(blst::Fp::equal(scalar_affine.Y[1], avx_y[1]));
}

TEST(AvxWrapper, G1GeneratorSerializationMatchesScalar)
{
    unsigned char scalar_bytes[96] = {};
    unsigned char avx_bytes[96] = {};

    blst::G1AffinePoint::serialize(scalar_bytes, blst::G1AffinePoint::generator());
    blst::avx::G1AffinePoint::serialize(avx_bytes, blst::avx::G1AffinePoint::generator());

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, G1DoublingMatchesScalar)
{
    blst::G1Point scalar_dbl = {};
    blst::G1Point::doubling(scalar_dbl, blst::G1Point::generator());

    blst::avx::G1Point avx_dbl = {};
    blst::avx::G1Point::doubling(avx_dbl, blst::avx::G1Point::generator());

    unsigned char scalar_bytes[96] = {};
    unsigned char avx_bytes[96] = {};
    blst::G1Point::serialize(scalar_bytes, scalar_dbl);
    blst::avx::G1Point::serialize(avx_bytes, avx_dbl);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, G1ScalarMultMatchesScalar)
{
    const std::array<blst::byte, 32> scalar = {1, 2, 3, 4, 5, 6, 7, 8};

    blst::G1Point scalar_out = {};
    blst::G1Point::mult(scalar_out, blst::G1Point::generator(), scalar.data(), 256);

    blst::avx::G1Point avx_out = {};
    blst::avx::G1Point::mult(avx_out, blst::avx::G1Point::generator(), scalar.data(), 256);

    unsigned char scalar_bytes[96] = {};
    unsigned char avx_bytes[96] = {};
    blst::G1Point::serialize(scalar_bytes, scalar_out);
    blst::avx::G1Point::serialize(avx_bytes, avx_out);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, G2GeneratorSerializationMatchesScalar)
{
    unsigned char scalar_bytes[192] = {};
    unsigned char avx_bytes[192] = {};

    blst::G2AffinePoint::serialize(scalar_bytes, blst::G2AffinePoint::generator());
    blst::avx::G2AffinePoint::serialize(avx_bytes, blst::avx::G2AffinePoint::generator());

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, G2DoublingMatchesScalar)
{
    blst::G2Point scalar_dbl = {};
    blst::G2Point::doubling(scalar_dbl, blst::G2Point::generator());

    blst::avx::G2Point avx_dbl = {};
    blst::avx::G2Point::doubling(avx_dbl, blst::avx::G2Point::generator());

    unsigned char scalar_bytes[192] = {};
    unsigned char avx_bytes[192] = {};
    blst::G2Point::serialize(scalar_bytes, scalar_dbl);
    blst::avx::G2Point::serialize(avx_bytes, avx_dbl);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, G2ScalarMultMatchesScalar)
{
    const std::array<blst::byte, 32> scalar = {9, 8, 7, 6, 5, 4, 3, 2};

    blst::G2Point scalar_out = {};
    blst::G2Point::mult(scalar_out, blst::G2Point::generator(), scalar.data(), 256);

    blst::avx::G2Point avx_out = {};
    blst::avx::G2Point::mult(avx_out, blst::avx::G2Point::generator(), scalar.data(), 256);

    unsigned char scalar_bytes[192] = {};
    unsigned char avx_bytes[192] = {};
    blst::G2Point::serialize(scalar_bytes, scalar_out);
    blst::avx::G2Point::serialize(avx_bytes, avx_out);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, GTRoundTripSerializationMatchesScalar)
{
    const blst::GT scalar_gt = make_scalar_final_pairing(
        make_scalar_g1_generator_affine(),
        make_scalar_g2_generator_affine());

    unsigned char scalar_bytes[12 * 48] = {};
    blst::GT::serialize(scalar_bytes, scalar_gt);

    blst::avx::GT avx_gt = {};
    blst::avx::GT::deserialize(avx_gt, scalar_bytes);

    unsigned char avx_bytes[12 * 48] = {};
    blst::avx::GT::serialize(avx_bytes, avx_gt);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, GTOperationsMatchScalar)
{
    const blst::GT scalar_a = make_scalar_final_pairing(
        make_scalar_g1_generator_affine(),
        make_scalar_g2_generator_affine());
    const blst::GT scalar_b = make_scalar_final_pairing(
        make_scalar_g1_doubled_affine(),
        make_scalar_g2_generator_affine());

    unsigned char scalar_a_bytes[12 * 48] = {};
    unsigned char scalar_b_bytes[12 * 48] = {};
    blst::GT::serialize(scalar_a_bytes, scalar_a);
    blst::GT::serialize(scalar_b_bytes, scalar_b);

    blst::avx::GT avx_a = {};
    blst::avx::GT avx_b = {};
    blst::avx::GT::deserialize(avx_a, scalar_a_bytes);
    blst::avx::GT::deserialize(avx_b, scalar_b_bytes);

    blst::GT scalar_mul = {};
    blst::GT::add(scalar_mul, scalar_a, scalar_b);

    blst::avx::GT avx_mul = {};
    blst::avx::GT::add(avx_mul, avx_a, avx_b);

    unsigned char scalar_bytes[12 * 48] = {};
    unsigned char avx_bytes[12 * 48] = {};
    blst::GT::serialize(scalar_bytes, scalar_mul);
    blst::avx::GT::serialize(avx_bytes, avx_mul);
    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));

    blst::GT scalar_sq = {};
    blst::GT::doubling(scalar_sq, scalar_a);

    blst::avx::GT avx_sq = {};
    blst::avx::GT::doubling(avx_sq, avx_a);

    blst::GT::serialize(scalar_bytes, scalar_sq);
    blst::avx::GT::serialize(avx_bytes, avx_sq);
    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, GTScalarMultMatchesScalar)
{
    const blst::GT scalar_a = make_scalar_final_pairing(
        make_scalar_g1_generator_affine(),
        make_scalar_g2_generator_affine());
    const std::array<blst::byte, 32> scalar = {1, 0, 2, 0, 3, 0, 4, 0};

    unsigned char scalar_a_bytes[12 * 48] = {};
    blst::GT::serialize(scalar_a_bytes, scalar_a);

    blst::avx::GT avx_a = {};
    blst::avx::GT::deserialize(avx_a, scalar_a_bytes);

    blst::GT scalar_out = {};
    blst::GT::mult(scalar_out, scalar_a, scalar.data(), 256);

    blst::avx::GT avx_out = {};
    blst::avx::GT::mult(avx_out, avx_a, scalar.data(), 256);

    unsigned char scalar_bytes[12 * 48] = {};
    unsigned char avx_bytes[12 * 48] = {};
    blst::GT::serialize(scalar_bytes, scalar_out);
    blst::avx::GT::serialize(avx_bytes, avx_out);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, GTMillerLoopMatchesScalar)
{
    const std::array<blst::G1AffinePoint, 3> scalar_ps = {
        make_scalar_g1_generator_affine(),
        make_scalar_g1_doubled_affine(),
        blst::G1AffinePoint{},
    };
    const std::array<blst::G2AffinePoint, 3> scalar_qs = {
        make_scalar_g2_generator_affine(),
        make_scalar_g2_doubled_affine(),
        blst::G2AffinePoint{},
    };

    for (size_t i = 0; i < scalar_ps.size(); ++i) {
        SCOPED_TRACE(i);

        blst::GT scalar_miller = {};
        blst::GT::miller_loop(scalar_miller, &scalar_qs[i], &scalar_ps[i]);
        blst::GT scalar_final = scalar_miller;
        blst::GT::final_mult(scalar_final, scalar_final);

        blst::avx::GT avx_miller = {};
        const blst::avx::G1AffinePoint avx_p = scalar_to_avx(scalar_ps[i]);
        const blst::avx::G2AffinePoint avx_q = scalar_to_avx(scalar_qs[i]);
        blst::avx::GT::miller_loop(avx_miller, &avx_q, &avx_p);

        blst::avx::GT avx_final = {};
        blst::avx::GT::final_mult(avx_final, avx_miller);

        unsigned char scalar_bytes[12 * 48] = {};
        unsigned char avx_bytes[12 * 48] = {};
        blst::GT::serialize(scalar_bytes, scalar_final);
        blst::avx::GT::serialize(avx_bytes, avx_final);

        EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
    }
}

TEST(AvxWrapper, GTMillerLoopNMatchesScalar)
{
    const std::array<blst::G1AffinePoint, 2> scalar_ps = {
        make_scalar_g1_generator_affine(),
        make_scalar_g1_doubled_affine(),
    };
    const std::array<blst::G2AffinePoint, 2> scalar_qs = {
        make_scalar_g2_generator_affine(),
        make_scalar_g2_doubled_affine(),
    };

    blst::GT scalar_multi = {};
    blst::GT::miller_loop_n(scalar_multi, scalar_qs.data(), scalar_ps.data(), scalar_ps.size());
    blst::GT::final_mult(scalar_multi, scalar_multi);

    const std::array<blst::avx::G1AffinePoint, 2> avx_ps = {
        scalar_to_avx(scalar_ps[0]),
        scalar_to_avx(scalar_ps[1]),
    };
    const std::array<blst::avx::G2AffinePoint, 2> avx_qs = {
        scalar_to_avx(scalar_qs[0]),
        scalar_to_avx(scalar_qs[1]),
    };

    blst::avx::GT avx_multi = {};
    blst::avx::GT::miller_loop_n(avx_multi, avx_qs.data(), avx_ps.data(), avx_ps.size());
    blst::avx::GT::final_mult(avx_multi, avx_multi);

    unsigned char scalar_bytes[12 * 48] = {};
    unsigned char avx_bytes[12 * 48] = {};
    blst::GT::serialize(scalar_bytes, scalar_multi);
    blst::avx::GT::serialize(avx_bytes, avx_multi);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, GTMillerLoopNPointerMatchesScalar)
{
    const std::array<blst::G1AffinePoint, 2> scalar_ps = {
        make_scalar_g1_generator_affine(),
        make_scalar_g1_doubled_affine(),
    };
    const std::array<blst::G2AffinePoint, 2> scalar_qs = {
        make_scalar_g2_generator_affine(),
        make_scalar_g2_doubled_affine(),
    };
    const std::array<const blst::G1AffinePoint*, 2> scalar_ps_ptrs = {&scalar_ps[0], nullptr};
    const std::array<const blst::G2AffinePoint*, 2> scalar_qs_ptrs = {&scalar_qs[0], nullptr};

    blst::GT scalar_multi = {};
    blst::GT::miller_loop_n(scalar_multi, scalar_qs_ptrs.data(), scalar_ps_ptrs.data(), scalar_ps_ptrs.size());
    blst::GT::final_mult(scalar_multi, scalar_multi);

    const std::array<blst::avx::G1AffinePoint, 2> avx_ps = {
        scalar_to_avx(scalar_ps[0]),
        scalar_to_avx(scalar_ps[1]),
    };
    const std::array<blst::avx::G2AffinePoint, 2> avx_qs = {
        scalar_to_avx(scalar_qs[0]),
        scalar_to_avx(scalar_qs[1]),
    };
    const std::array<const blst::avx::G1AffinePoint*, 2> avx_ps_ptrs = {&avx_ps[0], nullptr};
    const std::array<const blst::avx::G2AffinePoint*, 2> avx_qs_ptrs = {&avx_qs[0], nullptr};

    blst::avx::GT avx_multi = {};
    blst::avx::GT::miller_loop_n(avx_multi, avx_qs_ptrs.data(), avx_ps_ptrs.data(), avx_ps_ptrs.size());
    blst::avx::GT::final_mult(avx_multi, avx_multi);

    unsigned char scalar_bytes[12 * 48] = {};
    unsigned char avx_bytes[12 * 48] = {};
    blst::GT::serialize(scalar_bytes, scalar_multi);
    blst::avx::GT::serialize(avx_bytes, avx_multi);

    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(avx_bytes)));
}

TEST(AvxWrapper, AprimeConstantsMatchScalar)
{
    const auto& scalar_e1 = blst::scalar::get_BLST_Aprime_E1();
    const auto& scalar_e2 = blst::scalar::get_BLST_Aprime_E2();
    const auto& avx_e1 = blst::avx::get_BLST_Aprime_E1();
    const auto& avx_e2 = blst::avx::get_BLST_Aprime_E2();

    EXPECT_TRUE(blst::Fp::equal(scalar_e1, blst::avx::avx_curve_detail::ring_fp_to_scalar(avx_e1)));

    const blst::Fp2 avx_e2_scalar = blst::avx::avx_curve_detail::ring_fp2_to_scalar(avx_e2);
    EXPECT_TRUE(blst::Fp::equal(scalar_e2[0], avx_e2_scalar[0]));
    EXPECT_TRUE(blst::Fp::equal(scalar_e2[1], avx_e2_scalar[1]));
}

TEST(AvxWrapper, G1DaddWithAprimeMatchesScalar)
{
    blst::scalar::G1Point scalar_a = {};
    blst::scalar::G1Point scalar_b = {};
    blst::scalar::map_to_isogenous<blst::scalar::G1Point>(scalar_a, blst::scalar::G1AffinePoint::generator().X);
    blst::scalar::map_to_isogenous<blst::scalar::G1Point>(scalar_b, blst::scalar::G1AffinePoint::generator().Y);

    const auto& scalar_a4 = blst::scalar::get_BLST_Aprime_E1();
    const auto& avx_a4 = blst::avx::get_BLST_Aprime_E1();
    const blst::avx::G1Point avx_a = scalar_to_avx(scalar_a);
    const blst::avx::G1Point avx_b = scalar_to_avx(scalar_b);

    blst::scalar::G1Point scalar_double = {};
    blst::scalar::G1Point::dadd(scalar_double, scalar_a, scalar_a, &scalar_a4);
    blst::avx::G1Point avx_double = {};
    blst::avx::G1Point::dadd(avx_double, avx_a, avx_a, &avx_a4);
    expect_same_affine(scalar_double, avx_double);

    blst::scalar::G1Point scalar_sum = {};
    blst::scalar::G1Point::dadd(scalar_sum, scalar_a, scalar_b, &scalar_a4);
    blst::avx::G1Point avx_sum = {};
    blst::avx::G1Point::dadd(avx_sum, avx_a, avx_b, &avx_a4);
    expect_same_affine(scalar_sum, avx_sum);
}

TEST(AvxWrapper, G2DaddWithAprimeMatchesScalar)
{
    blst::scalar::G2Point scalar_a = {};
    blst::scalar::G2Point scalar_b = {};
    blst::scalar::map_to_isogenous<blst::scalar::G2Point>(scalar_a, blst::scalar::G2AffinePoint::generator().X);
    blst::scalar::map_to_isogenous<blst::scalar::G2Point>(scalar_b, blst::scalar::G2AffinePoint::generator().Y);

    const auto& scalar_a4 = blst::scalar::get_BLST_Aprime_E2();
    const auto& avx_a4 = blst::avx::get_BLST_Aprime_E2();
    const blst::avx::G2Point avx_a = scalar_to_avx(scalar_a);
    const blst::avx::G2Point avx_b = scalar_to_avx(scalar_b);

    blst::scalar::G2Point scalar_double = {};
    blst::scalar::G2Point::dadd(scalar_double, scalar_a, scalar_a, &scalar_a4);
    blst::avx::G2Point avx_double = {};
    blst::avx::G2Point::dadd(avx_double, avx_a, avx_a, &avx_a4);
    expect_same_affine(scalar_double, avx_double);

    blst::scalar::G2Point scalar_sum = {};
    blst::scalar::G2Point::dadd(scalar_sum, scalar_a, scalar_b, &scalar_a4);
    blst::avx::G2Point avx_sum = {};
    blst::avx::G2Point::dadd(avx_sum, avx_a, avx_b, &avx_a4);
    expect_same_affine(scalar_sum, avx_sum);
}

TEST(AvxWrapper, G2PsiMatchesScalar)
{
    std::array<blst::scalar::G2Point, 4> scalar_inputs = {};
    scalar_inputs[0] = blst::scalar::G2Point::generator();

    const std::array<blst::byte, 32> scalar = {0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88};
    blst::scalar::G2Point::mult(scalar_inputs[1], blst::scalar::G2Point::generator(), scalar.data(), 256);

    blst::scalar::map_to_isogenous<blst::scalar::G2Point>(
        scalar_inputs[2],
        blst::scalar::G2AffinePoint::generator().X);

    scalar_inputs[3] = blst::scalar::G2Point::identity();

    for (size_t i = 0; i < scalar_inputs.size(); ++i) {
        SCOPED_TRACE(i);

        blst::scalar::G2Point scalar_out = {};
        blst::scalar::psi(scalar_out, scalar_inputs[i]);

        blst::avx::G2Point avx_out = {};
        const blst::avx::G2Point avx_in = scalar_to_avx(scalar_inputs[i]);
        blst::avx::psi(avx_out, avx_in);

        expect_same_affine(scalar_out, avx_out);
    }
}

TEST(AvxWrapper, G1TimesMinusZMatchesScalar)
{
    std::array<blst::scalar::G1Point, 4> scalar_inputs = {};
    scalar_inputs[0] = blst::scalar::G1Point::generator();

    const std::array<blst::byte, 32> scalar = {0xab, 0xcd, 0xef, 0x01, 0x23, 0x45, 0x67, 0x89};
    blst::scalar::G1Point::mult(scalar_inputs[1], blst::scalar::G1Point::generator(), scalar.data(), 256);

    blst::scalar::map_to_isogenous<blst::scalar::G1Point>(
        scalar_inputs[2],
        blst::scalar::G1AffinePoint::generator().X);

    scalar_inputs[3] = blst::scalar::G1Point::identity();

    for (size_t i = 0; i < scalar_inputs.size(); ++i) {
        SCOPED_TRACE(i);

        blst::scalar::G1Point scalar_out = {};
        blst::scalar::times_minus_z(scalar_out, scalar_inputs[i]);

        blst::avx::G1Point avx_out = {};
        const blst::avx::G1Point avx_in = scalar_to_avx(scalar_inputs[i]);
        blst::avx::times_minus_z(avx_out, avx_in);

        expect_same_affine(scalar_out, avx_out);
    }
}

TEST(AvxWrapper, G2TimesMinusZMatchesScalar)
{
    std::array<blst::scalar::G2Point, 4> scalar_inputs = {};
    scalar_inputs[0] = blst::scalar::G2Point::generator();

    const std::array<blst::byte, 32> scalar = {0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10};
    blst::scalar::G2Point::mult(scalar_inputs[1], blst::scalar::G2Point::generator(), scalar.data(), 256);

    blst::scalar::map_to_isogenous<blst::scalar::G2Point>(
        scalar_inputs[2],
        blst::scalar::G2AffinePoint::generator().X);

    scalar_inputs[3] = blst::scalar::G2Point::identity();

    for (size_t i = 0; i < scalar_inputs.size(); ++i) {
        SCOPED_TRACE(i);

        blst::scalar::G2Point scalar_out = {};
        blst::scalar::times_minus_z(scalar_out, scalar_inputs[i]);

        blst::avx::G2Point avx_out = {};
        const blst::avx::G2Point avx_in = scalar_to_avx(scalar_inputs[i]);
        blst::avx::times_minus_z(avx_out, avx_in);

        expect_same_affine(scalar_out, avx_out);
    }
}

TEST(AvxWrapper, G1EncodeAndHashToCurveMatchScalar)
{
    static constexpr unsigned char msg[] = "avx g1 map-to-group wrapper";
    static constexpr unsigned char dst[] = "BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_TEST";
    static constexpr unsigned char aug[] = "ctx";

    blst::scalar::G1Point scalar_encode = {};
    blst::scalar::Encode_to_curve<blst::scalar::G1Point>(
        scalar_encode, msg, sizeof(msg) - 1, dst, sizeof(dst) - 1, aug, sizeof(aug) - 1);
    blst::avx::G1Point avx_encode = {};
    blst::avx::Encode_to_curve(
        avx_encode, msg, sizeof(msg) - 1, dst, sizeof(dst) - 1, aug, sizeof(aug) - 1);
    expect_same_affine(scalar_encode, avx_encode);

    blst::scalar::G1Point scalar_hash = {};
    blst::scalar::Hash_to_curve<blst::scalar::G1Point>(
        scalar_hash, msg, sizeof(msg) - 1, dst, sizeof(dst) - 1, aug, sizeof(aug) - 1);
    blst::avx::G1Point avx_hash = {};
    blst::avx::Hash_to_curve(
        avx_hash, msg, sizeof(msg) - 1, dst, sizeof(dst) - 1, aug, sizeof(aug) - 1);
    expect_same_affine(scalar_hash, avx_hash);
}

TEST(AvxWrapper, G2EncodeAndHashToCurveMatchScalar)
{
    static constexpr unsigned char msg[] = "avx g2 map-to-group wrapper";
    static constexpr unsigned char dst[] = "BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_TEST";
    static constexpr unsigned char aug[] = "ctx";

    blst::scalar::G2Point scalar_encode = {};
    blst::scalar::Encode_to_curve<blst::scalar::G2Point>(
        scalar_encode, msg, sizeof(msg) - 1, dst, sizeof(dst) - 1, aug, sizeof(aug) - 1);
    blst::avx::G2Point avx_encode = {};
    blst::avx::Encode_to_curve(
        avx_encode, msg, sizeof(msg) - 1, dst, sizeof(dst) - 1, aug, sizeof(aug) - 1);
    expect_same_affine(scalar_encode, avx_encode);

    blst::scalar::G2Point scalar_hash = {};
    blst::scalar::Hash_to_curve<blst::scalar::G2Point>(
        scalar_hash, msg, sizeof(msg) - 1, dst, sizeof(dst) - 1, aug, sizeof(aug) - 1);
    blst::avx::G2Point avx_hash = {};
    blst::avx::Hash_to_curve(
        avx_hash, msg, sizeof(msg) - 1, dst, sizeof(dst) - 1, aug, sizeof(aug) - 1);
    expect_same_affine(scalar_hash, avx_hash);
}

}  // namespace

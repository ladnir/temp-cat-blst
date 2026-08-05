#include "gtest/gtest.h"

#include "blst/curve/point.h"

#include "blst/curve/pairing.h"
#include "blst/fields/Fp12.h"
#include "blst/sign/aggregate.h"

using namespace blst;

namespace blst {
BLST_HD void blst_pairing_init(PAIRING* ctx, int hash_or_encode, const void* DST, size_t DST_len);
BLST_HD void blst_pairing_commit(PAIRING* ctx);
BLST_HD int blst_pairing_finalverify(const PAIRING* ctx, const Fp12* GTsig);
BLST_HD void blst_pairing_raw_aggregate(PAIRING* ctx, const G2AffinePoint* q, const G1AffinePoint* p);
BLST_HD void blst_pairing_get_miller_result(Fp12* out, PAIRING* ctx);
BLST_HD BLST_ERROR blst_aggregate_in_g1(G1Point& out, const G1Point* in, const unsigned char* zwire);
BLST_HD BLST_ERROR blst_aggregate_in_g2(G2Point& out, const G2Point* in, const unsigned char* zwire);
}

class GTTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize points
        g1_point = G1Point::generator();
        g2_point = G2Point::generator();

        // Convert to affine
        G1Point::to_affine(g1_affine, g1_point);
        G2Point::to_affine(g2_affine, g2_point);

        // Calculate pairing
        MillerResult miller;
        GT::miller_loop(miller, &g2_affine, &g1_affine);
        GT::final_exp(gt_element, miller);
    }

    G1Point g1_point;
    G2Point g2_point;
    G1AffinePoint g1_affine;
    G2AffinePoint g2_affine;
    GT gt_element;  // An element in GT
    GT temp;        // For temporary results
};

TEST_F(GTTest, Identity) {
    // Test identity element in GT
    auto identity = GT::identity();

    EXPECT_TRUE(GT::is_identity(identity))
        << "Identity element should be 0";

    // Adding any element to identity should leave it unchanged
    GT::add(temp, gt_element, identity);
    EXPECT_TRUE(GT::equal(temp, gt_element))
        << "a + 0 should equal a in GT";
}

TEST_F(GTTest, GeneratorMatchesCanonicalPairing) {
    const GT& generator = GT::generator();

    EXPECT_TRUE(GT::equal(generator, gt_element));
    EXPECT_FALSE(GT::is_identity(generator));
}

TEST_F(GTTest, Addition) {
    // Create another element for testing
    GT another_gt;
    G1Point g1_double;
    G1Point::doubling(g1_double, g1_point);

    G1AffinePoint g1_double_affine;
    G1Point::to_affine(g1_double_affine, g1_double);

    MillerResult another_miller;
    GT::miller_loop(another_miller, &g2_affine, &g1_double_affine);
    GT::final_exp(another_gt, another_miller);

    // Test addition
    GT sum;
    GT::add(sum, gt_element, another_gt);

    // The sum should not be equal to either of the original elements
    EXPECT_FALSE(GT::equal(sum, gt_element))
        << "Sum should not equal first element";
    EXPECT_FALSE(GT::equal(sum, another_gt))
        << "Sum should not equal second element";

    // Test commutativity: a + b = b + a
    GT sum2;
    GT::add(sum2, another_gt, gt_element);
    EXPECT_TRUE(GT::equal(sum, sum2))
        << "Addition should be commutative in GT";

    // Test associativity: (a + b) + c = a + (b + c)
    // Create a third element
    GT third_gt;
    G1Point g1_triple;
    G1Point::add(g1_triple, g1_double, g1_point);

    G1AffinePoint g1_triple_affine;
    G1Point::to_affine(g1_triple_affine, g1_triple);

    MillerResult third_miller;
    GT::miller_loop(third_miller, &g2_affine, &g1_triple_affine);
    GT::final_exp(third_gt, third_miller);

    GT temp1, temp2;
    // (a + b) + c
    GT::add(temp1, gt_element, another_gt);
    GT::add(temp1, temp1, third_gt);

    // a + (b + c)
    GT::add(temp2, another_gt, third_gt);
    GT::add(temp2, gt_element, temp2);

    EXPECT_TRUE(GT::equal(temp1, temp2))
        << "Addition should be associative in GT";
}

TEST_F(GTTest, Negation) {
    // Test negation (a + (-a) = 0)
    GT negative;
    GT::neg(negative, gt_element);

    GT sum;
    GT::add(sum, gt_element, negative);

    auto identity = GT::identity();

    EXPECT_TRUE(GT::equal(sum, identity))
        << "Element added to its negative should equal identity";
}

TEST_F(GTTest, Doubling) {
    // Test that doubling equals adding to self
    GT doubled;
    GT::doubling(doubled, gt_element);

    GT sum;
    GT::add(sum, gt_element, gt_element);

    EXPECT_TRUE(GT::equal(doubled, sum))
        << "2a should equal a + a";

    // Test multiple doubling (to test repeated doubling optimization if available)
    GT repeated_double;
    GT::doubling(repeated_double, gt_element);
    GT::doubling(repeated_double, repeated_double);

    // This should equal 4 * gt_element
    GT quadruple;
    GT::add(quadruple, gt_element, gt_element);
    GT::add(quadruple, quadruple, quadruple);

    EXPECT_TRUE(GT::equal(repeated_double, quadruple))
        << "2(2a) should equal 4a";
}

TEST_F(GTTest, ScalarMultiplication) {
    // Test scalar multiplication using small scalars
    GT result, expected;

    // Test 0 * e = 0
    byte zero[4] = { 0, 0, 0, 0 };
    GT::mult(result, gt_element, zero, sizeof(zero) * 8);
    GT expected_identity = GT::identity();
    EXPECT_TRUE(GT::equal(result, expected_identity))
        << "0 * e should equal 0";

    // Test 1 * e = e
    byte one[4] = { 1, 0, 0, 0 };
    GT::mult(result, gt_element, one, sizeof(one) * 8);
    EXPECT_TRUE(GT::equal(result, gt_element))
        << "1 * e should equal e";

    // Test 2 * e = e + e
    byte two[4] = { 2, 0, 0, 0 };
    GT::mult(result, gt_element, two, sizeof(two) * 8);
    GT::doubling(expected, gt_element);
    EXPECT_TRUE(GT::equal(result, expected))
        << "2 * e should equal e + e";

    // Test 3 * e = e + e + e
    byte three[4] = { 3, 0, 0, 0 };
    GT::mult(result, gt_element, three, sizeof(three) * 8);
    GT::add(expected, expected, gt_element);
    EXPECT_TRUE(GT::equal(result, expected))
        << "3 * e should equal e + e + e";
}

#if !BLST_DEFAULT_AVX
// More detailed diagnostic
TEST_F(GTTest, DiagnosticFinalExp) {
    // Manually trace through final_exp
    MillerResult miller_result;
    GT::miller_loop(miller_result, &g2_affine, &g1_affine);

    Fp12 f = miller_result.data;
    Fp12 y1, y2, ret;

    // Easy part: f^(p^6-1) * f^(p^2+1)
    vec_copy(&y1, &f, sizeof(y1));
    Fp12::conjugate(y1);                    // y1 = f^(p^6) = conjugate(f)
    Fp12::inverse(y2, f);                   // y2 = f^(-1)
    Fp12::mul(ret, y1, y2);                 // ret = f^(p^6-1)

    std::cout << "After f^(p^6-1) - is_cyclotomic: " << Fp12::is_cyclotomic(ret) << std::endl;

    Fp12::frobenius_map(y2, ret, 2);        // y2 = ret^(p^2)
    Fp12::mul(ret, ret, y2);                // ret = ret^(p^2+1) = f^((p^6-1)(p^2+1))

    std::cout << "After easy part (p^6-1)(p^2+1) - is_cyclotomic: " << Fp12::is_cyclotomic(ret) << std::endl;

    // At this point, ret should be cyclotomic!
    // If it's not, there's a bug in conjugate, inverse, mul, or frobenius_map

    EXPECT_TRUE(Fp12::is_cyclotomic(ret))
        << "Element should be cyclotomic after easy part of final_exp";

    // Now continue with hard part using cyclotomic_sqr
    if (Fp12::is_cyclotomic(ret)) {
        Fp12 y0;
        Fp12::cyclotomic_sqr(y0, ret);

        // Check if cyclotomic_sqr preserves cyclotomic property
        std::cout << "After cyclotomic_sqr - is_cyclotomic: " << Fp12::is_cyclotomic(y0) << std::endl;

        // Compare with regular sqr
        Fp12 y0_regular;
        Fp12::sqr(y0_regular, ret);

        bool sqr_match = Fp12::equal(y0, y0_regular);
        std::cout << "cyclotomic_sqr == sqr: " << sqr_match << std::endl;

        EXPECT_TRUE(sqr_match) << "cyclotomic_sqr should equal sqr for cyclotomic elements";
    }
}

// Add this diagnostic test to understand the issue
TEST_F(GTTest, DiagnosticCyclotomic) {
    // Check if the pairing result is cyclotomic after final_exp
    MillerResult miller_result;
    GT::miller_loop(miller_result, &g2_affine, &g1_affine);

    // Before final_exp - should NOT be cyclotomic
    bool before_cyclotomic = Fp12::is_cyclotomic(miller_result.data);
    std::cout << "Before final_exp - is_cyclotomic: " << before_cyclotomic << std::endl;

    // After final_exp - SHOULD be cyclotomic
    GT final_result;
    GT::final_exp(final_result, miller_result);
    bool after_cyclotomic = Fp12::is_cyclotomic(final_result.data);
    std::cout << "After final_exp - is_cyclotomic: " << after_cyclotomic << std::endl;

    // Check if final_result equals gt_element (from SetUp)
    bool equals_setup = GT::equal(final_result, gt_element);
    std::cout << "Equals gt_element from SetUp: " << equals_setup << std::endl;

    // Test that cyclotomic_sqr and regular sqr give same result on cyclotomic element
    if (after_cyclotomic) {
        Fp12 cyclo_sqr_result, regular_sqr_result;
        Fp12::cyclotomic_sqr(cyclo_sqr_result, final_result.data);
        Fp12::sqr(regular_sqr_result, final_result.data);

        bool sqr_match = Fp12::equal(cyclo_sqr_result, regular_sqr_result);
        std::cout << "cyclotomic_sqr == sqr: " << sqr_match << std::endl;

        EXPECT_TRUE(sqr_match) << "cyclotomic_sqr should equal sqr for cyclotomic elements";
    }

    EXPECT_TRUE(after_cyclotomic) << "Element should be cyclotomic after final_exp";
}
#endif
//TEST_F(GTTest, Serialization) {
//    // Test serialization/deserialization of GT elements
//    Fp12 original = { 0 };  // Copy the gt_element
//    memcpy(original, gt_element, sizeof(Fp12));
//
//    // Serialize to bytes
//    unsigned char serialized[12 * 48];  // Fp12 consists of 12 base field elements
//    blst_fp12_serialize(serialized, original);
//
//    // Deserialize back
//    Fp12 deserialized;
//    EXPECT_TRUE(blst_fp12_deserialize(deserialized, serialized))
//        << "Deserialization should succeed";
//
//    // The original and deserialized values should be equal
//    EXPECT_TRUE(GT::equal(original, deserialized))
//        << "Deserialized value should equal original";
//}

TEST_F(GTTest, BilinearityProperty) {
    // Test the bilinearity property: e(a*P, Q) = e(P, a*Q) = a*e(P,Q)

    // Create scaled points
    G1Point g1_double;
    G2Point g2_double;
    G1Point::doubling(g1_double, g1_point);
    G2Point::doubling(g2_double, g2_point);

    G1AffinePoint g1_double_affine;
    G2AffinePoint g2_double_affine;
    G1Point::to_affine(g1_double_affine, g1_double);
    G2Point::to_affine(g2_double_affine, g2_double);

    // e(2P, Q)
    MillerResult miller_2P_Q;
    GT e_2P_Q;
    GT::miller_loop(miller_2P_Q, &g2_affine, &g1_double_affine);
    GT::final_exp(e_2P_Q, miller_2P_Q);

    // e(P, 2Q)
    MillerResult miller_P_2Q;
    GT e_P_2Q;
    GT::miller_loop(miller_P_2Q, &g2_double_affine, &g1_affine);
    GT::final_exp(e_P_2Q, miller_P_2Q);

    // 2 * e(P,Q)
    GT doubled_e_PQ;
    GT::doubling(doubled_e_PQ, gt_element);

    // All three should be equal due to bilinearity
    EXPECT_TRUE(GT::equal(e_2P_Q, e_P_2Q))
        << "e(2P, Q) should equal e(P, 2Q)";
    EXPECT_TRUE(GT::equal(e_2P_Q, doubled_e_PQ))
        << "e(2P, Q) should equal 2*e(P,Q)";
}



// NEW TESTS BELOW

TEST_F(GTTest, GroupMembership) {
    // Test that the pairing result is in the correct group
    EXPECT_TRUE(GT::in_group(gt_element))
        << "Pairing result should be in the GT group";

    // Identity should be in the group
    GT identity = GT::identity();
    EXPECT_TRUE(GT::in_group(identity))
        << "Identity element should be in the GT group";

    // Doubling, addition, negation should preserve group membership
    GT doubled;
    GT::doubling(doubled, gt_element);
    EXPECT_TRUE(GT::in_group(doubled))
        << "Doubling should preserve group membership";

    GT negative;
    GT::neg(negative, gt_element);
    EXPECT_TRUE(GT::in_group(negative))
        << "Negation should preserve group membership";
}

TEST_F(GTTest, MultiPairing) {
    // Set up multiple points for testing
    G1Point g1_double;
    G1Point::doubling(g1_double, g1_point);

    G1AffinePoint g1_points[2];
    G1Point::to_affine(g1_points[0], g1_point);
    G1Point::to_affine(g1_points[1], g1_double);

    G2AffinePoint g2_points[2];
    G2Point::to_affine(g2_points[0], g2_point);
    G2Point::to_affine(g2_points[1], g2_point);

    // Create arrays of pointers as required by GT::miller_loop_n
    const G1AffinePoint* g1_ptrs[3] = { &g1_points[0], &g1_points[1], nullptr };
    const G2AffinePoint* g2_ptrs[3] = { &g2_points[0], &g2_points[1], nullptr };

    // Calculate e(P1,Q1) and e(P2,Q2) individually
    MillerResult e1, e2, product, multi_miller;
    GT sum, multi;
    GT::miller_loop(e1, &g2_points[0], &g1_points[0]);
    GT::miller_loop(e2, &g2_points[1], &g1_points[1]);

    // Calculate their sum
    MillerResult::mul(product, e1, e2);
    GT::final_exp(sum, product);

    // Calculate using multi-pairing
    GT::miller_loop_n(multi_miller, g2_ptrs, g1_ptrs, 2);
    GT::final_exp(multi, multi_miller);

    // They should match
    EXPECT_TRUE(GT::equal(sum, multi))
        << "Multi-pairing result should match individual pairings";
}

TEST_F(GTTest, FinalVerify) {
    // Test finalverify function (equality after final exponentiation)

    // Create two different pairings that should be equal after final_exp
    G1Point g1_double;
    G1Point::doubling(g1_double, g1_point);

    // Create a "half" point in G2 by multiplying by inverse of 2 (if available)
    // For test purposes, we'll just use a different point
    G2Point g2_another;
    G2Point::add(g2_another, g2_point, g2_point);  // Just a different point

    G1AffinePoint g1_double_affine;
    G2AffinePoint g2_half_affine;
    G1Point::to_affine(g1_double_affine, g1_double);
    G2Point::to_affine(g2_half_affine, g2_another);

    // Calculate two different pairings
    MillerResult pair1, pair2;
    GT::miller_loop(pair1, &g2_affine, &g1_double_affine);
    GT::miller_loop(pair2, &g2_half_affine, &g1_affine);

    // Apply final exponentiation
    GT final1, final2;
    GT::final_exp(final1, pair1);
    GT::final_exp(final2, pair2);

    // Direct equality check
    bool direct_equal = GT::equal(final1, final2);

    // Using finalverify
    bool finalverify_equal = blst_fp12_finalverify(pair1, pair2);

    // Verify that finalverify gives same result as equality check
    // (in general they might differ due to optimizations, but for our simple test they should match)
    EXPECT_EQ(direct_equal, finalverify_equal)
        << "finalverify should agree with is_equal in basic cases";
}

#if !BLST_DEFAULT_AVX
TEST_F(GTTest, PairingWithPrecomputation) {
    // Test using precomputed lines for pairing
    Fp6 lines[68];  // Size from the API

    // Precompute lines for pairing with g2_point
    GT::precompute_lines(lines, g2_affine);

    // Compute pairing using precomputed lines
    MillerResult pre_miller;
    GT pre_result;
    GT::miller_loop_lines(pre_miller, lines, g1_affine);
    GT::final_exp(pre_result, pre_miller);

    // Verify against standard pairing computation
    EXPECT_TRUE(GT::equal(gt_element, pre_result))
        << "Pairing with precomputation should match direct computation";

    // Verify that we can reuse the precomputed lines with a different point
    G1Point g1_double;
    G1Point::doubling(g1_double, g1_point);

    G1AffinePoint g1_double_affine;
    G1Point::to_affine(g1_double_affine, g1_double);

    MillerResult direct_double_miller;
    GT direct_double_result;
    GT::miller_loop(direct_double_miller, &g2_affine, &g1_double_affine);
    GT::final_exp(direct_double_result, direct_double_miller);

    MillerResult pre_double_miller;
    GT pre_double_result;
    GT::miller_loop_lines(pre_double_miller, lines, g1_double_affine);
    GT::final_exp(pre_double_result, pre_double_miller);

    EXPECT_TRUE(GT::equal(direct_double_result, pre_double_result))
        << "Precomputed pairing should work with different G1 points";
}
#endif

TEST_F(GTTest, MillerLoopOnly) {
    // Test Miller loop without final exponentiation
    MillerResult miller_result1, miller_result2;

    // Compute Miller loop for same pairing twice
    GT::miller_loop(miller_result1, &g2_affine, &g1_affine);
    GT::miller_loop(miller_result2, &g2_affine, &g1_affine);

    // Final exponentiation of identical Miller results must be deterministic.
    GT final_result, final_result2;
    GT::final_exp(final_result, miller_result1);
    GT::final_exp(final_result2, miller_result2);
    EXPECT_TRUE(GT::equal(final_result, final_result2));
    EXPECT_TRUE(GT::in_group(final_result))
        << "Result should be in GT after final exponentiation";
}

TEST_F(GTTest, PairingIdentityProperties) {
    // Test pairing with identity elements

    // Create identity point in G1 (point at infinity)
    G1Point g1_identity;
    memset(&g1_identity, 0, sizeof(g1_identity));
    // Set z coordinate to 0 to make it point at infinity
    //g1_identity.z[0] = 0;

    G1AffinePoint g1_identity_affine;
    G1Point::to_affine(g1_identity_affine, g1_identity);

    // Pairing with identity should give identity in GT
    MillerResult identity_miller;
    GT pairing_with_identity;
    GT::miller_loop(identity_miller, &g2_affine, &g1_identity_affine);
    GT::final_exp(pairing_with_identity, identity_miller);

    EXPECT_TRUE(GT::is_identity(pairing_with_identity))
        << "Pairing with identity point should give identity in GT";
}

TEST_F(GTTest, AggregateHelpersMatchDirectAddition) {
    G1Point g1_double;
    G1Point::doubling(g1_double, g1_point);
    G1AffinePoint g1_double_affine;
    G1Point::to_affine(g1_double_affine, g1_double);

    unsigned char g1_bytes[96] = {};
    unsigned char g1_double_bytes[96] = {};
    G1AffinePoint::serialize(g1_bytes, g1_affine);
    G1AffinePoint::serialize(g1_double_bytes, g1_double_affine);

    G1Point g1_agg = {};
    ASSERT_EQ(BLST_SUCCESS, blst_aggregate_in_g1(g1_agg, nullptr, g1_bytes));
    ASSERT_EQ(BLST_SUCCESS, blst_aggregate_in_g1(g1_agg, &g1_agg, g1_double_bytes));

    G1Point g1_expected = {};
    G1Point::add(g1_expected, g1_point, g1_double);
    EXPECT_TRUE(G1Point::equal(g1_agg, g1_expected));

    G2Point g2_double;
    G2Point::doubling(g2_double, g2_point);
    G2AffinePoint g2_double_affine;
    G2Point::to_affine(g2_double_affine, g2_double);

    unsigned char g2_bytes[192] = {};
    unsigned char g2_double_bytes[192] = {};
    G2AffinePoint::serialize(g2_bytes, g2_affine);
    G2AffinePoint::serialize(g2_double_bytes, g2_double_affine);

    G2Point g2_agg = {};
    ASSERT_EQ(BLST_SUCCESS, blst_aggregate_in_g2(g2_agg, nullptr, g2_bytes));
    ASSERT_EQ(BLST_SUCCESS, blst_aggregate_in_g2(g2_agg, &g2_agg, g2_double_bytes));

    G2Point g2_expected = {};
    G2Point::add(g2_expected, g2_point, g2_double);
    EXPECT_TRUE(G2Point::equal(g2_agg, g2_expected));
}

TEST_F(GTTest, PairingAggregationRoundTripMatchesDirectMiller) {
    PAIRING ctx = {};
    blst_pairing_init(&ctx, 0, nullptr, 0);
    blst_pairing_raw_aggregate(&ctx, &g2_affine, &g1_affine);
    blst_pairing_commit(&ctx);

    MillerResult direct = {};
    GT::miller_loop(direct, &g2_affine, &g1_affine);

    Fp12 aggregated_fp12 = {};
    blst_pairing_get_miller_result(&aggregated_fp12, &ctx);
    unsigned char aggregated_bytes[12 * 48] = {};
    scalar::GT::serialize(aggregated_bytes, scalar::GT(aggregated_fp12));

    EXPECT_TRUE(blst_pairing_finalverify(&ctx, &aggregated_fp12));

    GT direct_final = {};
    GT::final_exp(direct_final, direct);
    GT aggregated_final = {};
    MillerResult aggregated_miller = {};
#if BLST_DEFAULT_AVX
    aggregated_miller.data = blst::avx::avx_gt_detail::from_scalar(aggregated_fp12);
#else
    aggregated_miller.data = aggregated_fp12;
#endif
    GT::final_exp(aggregated_final, aggregated_miller);
    EXPECT_TRUE(GT::equal(aggregated_final, direct_final));
}

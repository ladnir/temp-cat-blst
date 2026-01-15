#include "gtest/gtest.h"

#include "blst/curve/point.h"

#include "blst/curve/pairing.h"
#include "blst/fields/Fp12.h"
#include "blst/sign/aggregate.h"

using namespace blst;

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
        GT::miller_loop(gt_element, &g2_affine, &g1_affine);
        GT::final_mult(gt_element, gt_element);
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

TEST_F(GTTest, Addition) {
    // Create another element for testing
    GT another_gt;
    G1Point g1_double;
    G1Point::doubling(g1_double, g1_point);

    G1AffinePoint g1_double_affine;
    G1Point::to_affine(g1_double_affine, g1_double);

    GT::miller_loop(another_gt, &g2_affine, &g1_double_affine);
    GT::final_mult(another_gt, another_gt);

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

    GT::miller_loop(third_gt, &g2_affine, &g1_triple_affine);
    GT::final_mult(third_gt, third_gt);

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

    auto identity = Fp12::one();

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
    Fp12 expected_identity = { 0 };
    expected_identity = Fp12::one();
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

// More detailed diagnostic
TEST_F(GTTest, DiagnosticFinalExp) {
    // Manually trace through final_exp
    GT miller_result;
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
    GT miller_result;
    GT::miller_loop(miller_result, &g2_affine, &g1_affine);

    // Before final_exp - should NOT be cyclotomic
    bool before_cyclotomic = Fp12::is_cyclotomic(miller_result.data);
    std::cout << "Before final_exp - is_cyclotomic: " << before_cyclotomic << std::endl;

    // After final_exp - SHOULD be cyclotomic
    GT final_result;
    GT::final_mult(final_result, miller_result);
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
    GT e_2P_Q;
    GT::miller_loop(e_2P_Q, &g2_affine, &g1_double_affine);
    GT::final_mult(e_2P_Q, e_2P_Q);

    // e(P, 2Q)
    GT e_P_2Q;
    GT::miller_loop(e_P_2Q, &g2_double_affine, &g1_affine);
    GT::final_mult(e_P_2Q, e_P_2Q);

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
    GT identity;
    identity = Fp12::one();
    EXPECT_TRUE(GT::in_group(identity))
        << "Identity element should be in the GT group";

    // Test that a non-GT element (before final_exp) is not in the group
    GT not_finalized;
    GT::miller_loop(not_finalized, &g2_affine, &g1_affine);
    EXPECT_FALSE(GT::in_group(not_finalized))
        << "Miller loop result without final_exp should not be in GT group";

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
    GT e1, e2, sum, multi;
    GT::miller_loop(e1, &g2_points[0], &g1_points[0]);
    GT::miller_loop(e2, &g2_points[1], &g1_points[1]);

    // Calculate their sum
    GT::add(sum, e1, e2);
    GT::final_mult(sum, sum);

    // Calculate using multi-pairing
    GT::miller_loop_n(multi, g2_ptrs, g1_ptrs, 2);
    GT::final_mult(multi, multi);

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
    GT pair1, pair2;
    GT::miller_loop(pair1, &g2_affine, &g1_double_affine);
    GT::miller_loop(pair2, &g2_half_affine, &g1_affine);

    // They should not be equal before final_exp
    EXPECT_FALSE(GT::equal(pair1, pair2))
        << "Different pairings should not be equal before final_exp";

    // Apply final exponentiation
    GT final1, final2;
    GT::final_mult(final1, pair1);
    GT::final_mult(final2, pair2);

    // Direct equality check
    bool direct_equal = GT::equal(final1, final2);

    // Using finalverify
    bool finalverify_equal = blst_fp12_finalverify(final1, final2);

    // Verify that finalverify gives same result as equality check
    // (in general they might differ due to optimizations, but for our simple test they should match)
    EXPECT_EQ(direct_equal, finalverify_equal)
        << "finalverify should agree with is_equal in basic cases";
}

TEST_F(GTTest, PairingWithPrecomputation) {
    // Test using precomputed lines for pairing
    Fp6 lines[68];  // Size from the API

    // Precompute lines for pairing with g2_point
    GT::precompute_lines(lines, g2_affine);

    // Compute pairing using precomputed lines
    GT pre_result;
    GT::miller_loop_lines(pre_result, lines, g1_affine);
    GT::final_mult(pre_result, pre_result);

    // Verify against standard pairing computation
    EXPECT_TRUE(GT::equal(gt_element, pre_result))
        << "Pairing with precomputation should match direct computation";

    // Verify that we can reuse the precomputed lines with a different point
    G1Point g1_double;
    G1Point::doubling(g1_double, g1_point);

    G1AffinePoint g1_double_affine;
    G1Point::to_affine(g1_double_affine, g1_double);

    GT direct_double_result;
    GT::miller_loop(direct_double_result, &g2_affine, &g1_double_affine);
    GT::final_mult(direct_double_result, direct_double_result);

    GT pre_double_result;
    GT::miller_loop_lines(pre_double_result, lines, g1_double_affine);
    GT::final_mult(pre_double_result, pre_double_result);

    EXPECT_TRUE(GT::equal(direct_double_result, pre_double_result))
        << "Precomputed pairing should work with different G1 points";
}

TEST_F(GTTest, MillerLoopOnly) {
    // Test Miller loop without final exponentiation
    GT miller_result1, miller_result2;

    // Compute Miller loop for same pairing twice
    GT::miller_loop(miller_result1, &g2_affine, &g1_affine);
    GT::miller_loop(miller_result2, &g2_affine, &g1_affine);

    // Should be equal
    EXPECT_TRUE(GT::equal(miller_result1, miller_result2))
        << "Miller loop should be deterministic";

    // Should not be in GT group before final exp
    EXPECT_FALSE(GT::in_group(miller_result1))
        << "Miller loop result should not be in GT before final exponentiation";

    // Apply final exp and verify it's now in GT
    GT final_result;
    GT::final_mult(final_result, miller_result1);
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
    GT pairing_with_identity;
    GT::miller_loop(pairing_with_identity, &g2_affine, &g1_identity_affine);
    GT::final_mult(pairing_with_identity, pairing_with_identity);

    EXPECT_TRUE(GT::is_identity(pairing_with_identity))
        << "Pairing with identity point should give identity in GT";
}
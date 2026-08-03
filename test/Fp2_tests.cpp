#include <gtest/gtest.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "blst/fields/Fp2.h"
using namespace blst;

// Helper: Print Fp2 element for debugging
void print_fp2(const Fp2& a, const char* name) {
    printf("%s:\n", name);
    printf("  Re: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)a[0][i]);
    }
    printf("\n  Im: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)a[1][i]);
    }
    printf("\n");
}




TEST(Fp2Test, ZeroAndOne) {
    Fp2 zero, one;

    // Initialize zero and one
    vec_zero(&zero, sizeof(zero));
    vec_zero(&one, sizeof(one));
    vec_copy(one[0].l, Fp::one().l, sizeof(Fp));  // Re(1) = 1

    // Check zero properties
    Fp2 result;
    Fp2::add(result, zero, zero);
    EXPECT_TRUE(vec_is_equal(&result, &zero, sizeof(Fp2))) << "0 + 0 = 0";

    Fp2::mul(result, zero, one);
    EXPECT_TRUE(vec_is_equal(&result, &zero, sizeof(Fp2))) << "0 * 1 = 0";

    // Check one properties
    Fp2::mul(result, one, one);
    EXPECT_TRUE(vec_is_equal(&result, &one, sizeof(Fp2))) << "1 * 1 = 1";

    Fp2 temp;
    temp = one;
    Fp2::add(result, temp, temp);

    // 1 + 1 should be 2 (in Montgomery representation)
    Fp2 expected_two;
    vec_zero(&expected_two, sizeof(expected_two));
    Fp::lshift(expected_two[0], one[0], 1);  // Re(2) = 2

    EXPECT_TRUE(vec_is_equal(&result, &expected_two, sizeof(Fp2))) << "1 + 1 = 2";
}

TEST(Fp2Test, OneMultiplication) {
    Fp2 one = Fp2::one(), result;

    // Initialize `one` to Fp2::one()
    //vec_zero(&one, sizeof(one));
    //vec_copy(one[0], Fp::one(), sizeof(Fp));  // Re(1) = 1
    //vec_zero(one[1], sizeof(Fp));                  // Im(1) = 0

    // Compute result = one * one
    Fp2::mul(result, one, one);

    // Verify that result == one
    EXPECT_TRUE(vec_is_equal(&result, &one, sizeof(Fp2)))
        << "Fp2::one() * Fp2::one() should equal Fp2::one()";
}

TEST(Fp2Test, Addition) {
    Fp2 a, b, c, expected;

    // Initialize with some values
    // a = 2 + 3i (in Montgomery representation)
    vec_zero(&a, sizeof(a));
    Fp::lshift(a[0], Fp::one(), 1);     // Re(a) = 2
    Fp::mul_by_3(a[1], Fp::one());      // Im(a) = 3

    // b = 4 + 5i (in Montgomery representation)
    vec_zero(&b, sizeof(b));
    Fp::lshift(b[0], Fp::one(), 2);     // Re(b) = 4
    Fp::add(b[1], Fp::one(), b[0]);     // Im(b) = 1 + 4 = 5

    // c = a + b should be 6 + 8i
    Fp2::add(c, a, b);

    vec_zero(&expected, sizeof(expected));
    // Calculate Re(expected) = 6
    Fp::mul_by_3(expected[0], Fp::one());       // Re(expected) = 3
    Fp::lshift(expected[0], expected[0], 1);         // Re(expected) = 6
    // Calculate Im(expected) = 8
    Fp::lshift(expected[1], Fp::one(), 3);      // Im(expected) = 8

    EXPECT_TRUE(vec_is_equal(&c, &expected, sizeof(Fp2)))
        << "(2+3i) + (4+5i) = (6+8i)";

    // Test commutativity
    Fp2::add(c, b, a);
    EXPECT_TRUE(vec_is_equal(&c, &expected, sizeof(Fp2)))
        << "(4+5i) + (2+3i) = (6+8i)";
}

TEST(Fp2Test, Subtraction) {
    Fp2 a, b, c, expected;

    // Initialize with some values
    // a = 7 + 4i (in Montgomery representation)
    vec_zero(&a, sizeof(a));
    Fp::mul_by_3(a[0], Fp::one());      // Re(a) = 3
    Fp::add(a[0], a[0], Fp::one());     // Re(a) = 3+1 = 4
    Fp::add(a[0], a[0], Fp::one());     // Re(a) = 4+1 = 5
    Fp::add(a[0], a[0], Fp::one());     // Re(a) = 5+1 = 6
    Fp::add(a[0], a[0], Fp::one());     // Re(a) = 6+1 = 7
    Fp::lshift(a[1], Fp::one(), 2);     // Im(a) = 4

    // b = 2 + 3i (in Montgomery representation)
    vec_zero(&b, sizeof(b));
    Fp::lshift(b[0], Fp::one(), 1);     // Re(b) = 2
    Fp::mul_by_3(b[1], Fp::one());      // Im(b) = 3

    // c = a - b should be 5 + 1i
    Fp2::sub(c, a, b);

    vec_zero(&expected, sizeof(expected));
    Fp::add(expected[0], Fp::one(), Fp::one());  // Re(expected) = 2
    Fp::add(expected[0], expected[0], Fp::one());     // Re(expected) = 2+1 = 3
    Fp::add(expected[0], expected[0], Fp::one());     // Re(expected) = 3+1 = 4
    Fp::add(expected[0], expected[0], Fp::one());     // Re(expected) = 4+1 = 5
    vec_copy(expected[1].l, Fp::one().l, sizeof(Fp));  // Im(expected) = 1

    EXPECT_TRUE(vec_is_equal(&c, &expected, sizeof(Fp2)))
        << "(7+4i) - (2+3i) = (5+1i)";

    // c = b - a should be -5 - 1i
    Fp2::sub(c, b, a);

    // Negate expected to get -5 - 1i
    Fp2::cneg(expected, expected, 1);

    EXPECT_TRUE(vec_is_equal(&c, &expected, sizeof(Fp2)))
        << "(2+3i) - (7+4i) = (-5-1i)";
}

TEST(Fp2Test, Multiplication) {
    Fp2 a, b, c, expected;

    // a = 2 + 3i
    vec_zero(&a, sizeof(a));
    Fp::lshift(a[0], Fp::one(), 1);     // Re(a) = 2
    Fp::mul_by_3(a[1], Fp::one());      // Im(a) = 3

    // b = 4 + 5i
    vec_zero(&b, sizeof(b));
    Fp::lshift(b[0], Fp::one(), 2);     // Re(b) = 4
    Fp::add(b[1], Fp::one(), b[0]);     // Im(b) = 1 + 4 = 5

    // c = a * b should be (2*4 - 3*5) + (2*5 + 3*4)i = (8 - 15) + (10 + 12)i = -7 + 22i
    Fp2::mul(c, a, b);

    vec_zero(&expected, sizeof(expected));
    // Re(expected) = -7
    Fp::mul_by_3(expected[0], Fp::one());      // 3
    Fp::lshift(expected[0], expected[0], 1);        // 6
    Fp::add(expected[0], expected[0], Fp::one()); // 7
    Fp::cneg(expected[0], expected[0], 1);           // -7

    // Im(expected) = 22
    Fp temp;
    Fp::lshift(temp, Fp::one(), 1);             // 2
    Fp::mul_by_3(expected[1], Fp::one());       // 3
    Fp::lshift(expected[1], expected[1], 3);         // 24
    Fp::sub(expected[1], expected[1], temp);         // 22

    EXPECT_TRUE(vec_is_equal(&c, &expected, sizeof(Fp2)))
        << "(2+3i) * (4+5i) = (-7+22i)";

    // Test commutativity
    Fp2::mul(c, b, a);
    EXPECT_TRUE(vec_is_equal(&c, &expected, sizeof(Fp2)))
        << "(4+5i) * (2+3i) = (-7+22i)";
}

TEST(Fp2Test, SquareAndMultiply) {
    Fp2 a, a_squared1, a_squared2;

    // a = 2 + 3i
    vec_zero(&a, sizeof(a));
    Fp::lshift(a[0], Fp::one(), 1);     // Re(a) = 2
    Fp::mul_by_3(a[1], Fp::one());      // Im(a) = 3

    // Compute a² using multiplication
    Fp2::mul(a_squared1, a, a);

    // Compute a² using squaring
    Fp2::sqr(a_squared2, a);

    // Results should be equal: (2+3i)² = 4-9 + 12i = -5 + 12i
    EXPECT_TRUE(vec_is_equal(&a_squared1, &a_squared2, sizeof(Fp2)))
        << "(2+3i)² computed both ways should be equal";
}

TEST(Fp2Test, Inversion) {
    Fp2 a, a_inv, result;

    // a = 2 + 3i
    vec_zero(&a, sizeof(a));
    Fp::lshift(a[0], Fp::one(), 1);     // Re(a) = 2
    Fp::mul_by_3(a[1], Fp::one());      // Im(a) = 3

    // Compute 1/a
    Fp2::reciprocal(a_inv, a);

    // Check that a * (1/a) = 1
    Fp2::mul(result, a, a_inv);

    // Should equal 1 (real part 1, imaginary part 0)
    Fp2 one;
    vec_zero(&one, sizeof(one));
    vec_copy(one[0].l, Fp::one().l, sizeof(Fp));

    EXPECT_TRUE(vec_is_equal(&result, &one, sizeof(Fp2)))
        << "a * (1/a) should equal 1";
}

TEST(Fp2Test, Frobenius) {
    Fp2 a, frob_a, expected;

    // a = 2 + 3i
    vec_zero(&a, sizeof(a));
    Fp::lshift(a[0], Fp::one(), 1);     // Re(a) = 2
    Fp::mul_by_3(a[1], Fp::one());      // Im(a) = 3

    // Frobenius map in Fp2 is complex conjugation: (a+bi) -> (a-bi)
    // expected = 2 - 3i
    expected = a;
    Fp::cneg(expected[1], expected[1], 1);  // Negate imaginary part

    // Apply Frobenius
    Fp2::conjugate(frob_a, a);

    EXPECT_TRUE(vec_is_equal(&frob_a, &expected, sizeof(Fp2)))
        << "Frobenius(2+3i) should be 2-3i";
}

TEST(Fp2Test, SquareRoot) {
    Fp2 a, root, root_squared;

    // Let's try with a perfect square: a = (2+3i)²
    Fp2 temp;
    vec_zero(&temp, sizeof(temp));
    Fp::lshift(temp[0], Fp::one(), 1);    // Re(temp) = 2
    Fp::mul_by_3(temp[1], Fp::one());     // Im(temp) = 3

    Fp2::sqr(a, temp);  // a = (2+3i)²

    // Find square root
    EXPECT_TRUE(Fp2::sqrt(root, a)) << "Square root should succeed for a perfect square";

    // Verify by squaring the root
    Fp2::sqr(root_squared, root);

    // Should get back our original value
    EXPECT_TRUE(vec_is_equal(&root_squared, &a, sizeof(Fp2)))
        << "sqrt(a)² should equal a";

    // Verify that either root or -root equals our original value
    bool matches_orig = vec_is_equal(&root, &temp, sizeof(Fp2));
    Fp2::cneg(root, root, 1);
    bool matches_neg = vec_is_equal(&root, &temp, sizeof(Fp2));

    EXPECT_TRUE(matches_orig || matches_neg)
        << "Either sqrt(a) or -sqrt(a) should equal the original number";
}

TEST(Fp2Test, Serialization) {
    Fp2 a, b;
    unsigned char serialized[96]; // 2 * 48 bytes

    // a = 2 + 3i
    vec_zero(&a, sizeof(a));
    Fp::lshift(a[0], Fp::one(), 1);     // Re(a) = 2
    Fp::mul_by_3(a[1], Fp::one());      // Im(a) = 3

    // Serialize to bytes: convert from Montgomery to standard form
    Fp temp_re, temp_im;
    Fp::from_mont(temp_re, a[0]);  // Convert real part from Montgomery
    Fp::from_mont(temp_im, a[1]);  // Convert imaginary part from Montgomery

    // Copy to serialized buffer
    memcpy(serialized, temp_re.l, 48);
    memcpy(serialized + 48, temp_im.l, 48);

    // Deserialize back: copy from buffer then convert to Montgomery form
    Fp raw_re, raw_im;
    memcpy(raw_re.l, serialized, 48);
    memcpy(raw_im.l, serialized + 48, 48);

    // Convert to Montgomery form by multiplying by R^2
    Fp::mul(b[0], raw_re, get_BLST_RR());
    Fp::mul(b[1], raw_im, get_BLST_RR());

    // Should get back our original value
    EXPECT_TRUE(vec_is_equal(&a, &b, sizeof(Fp2)))
        << "Deserialized value should equal original";
}

TEST(Fp2Test, ZeroBitLShiftIsIdentity) {
    Fp2 value = Fp2::one(), result;
    Fp::mul_by_3(value[1], value[0]);

    Fp2::lshift(result, value, 0);
    EXPECT_TRUE(vec_is_equal(&result, &value, sizeof(value)));

    Fp2::lshift(value, value, 0);
    EXPECT_TRUE(vec_is_equal(&result, &value, sizeof(value)));
}

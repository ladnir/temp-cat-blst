#include <gtest/gtest.h>

#include <array>
#include <cstddef>

#include "blst/scalar/curve/point.h"

namespace {

using blst::byte;
using blst::bool_t;
using blst::scalar::G1Point;
using blst::scalar::G2Point;

template<typename Point>
Point generator_multiple(byte scalar)
{
    Point out = {};
    Point::mult(out, Point::generator(), &scalar, 8);
    return out;
}

template<typename Point>
void expect_all_mult_window_modes()
{
    constexpr size_t npoints = 3;
    constexpr int window_size = 5;
    constexpr size_t table_width = size_t{1} << (window_size - 1);

    std::array<Point, npoints> points = {
        generator_multiple<Point>(1),
        generator_multiple<Point>(2),
        generator_multiple<Point>(3),
    };
    const Point* compact_points[] = {points.data(), nullptr};

    const byte scalar_storage[npoints] = {3, 5, 7};
    const byte* compact_scalars[] = {scalar_storage, nullptr};
    const Point expected = generator_multiple<Point>(34);

    Point allocated = {};
    ASSERT_TRUE((blst::scalar::mult_window_n<Point, window_size>(
        allocated,
        compact_points,
        npoints,
        compact_scalars,
        8,
        nullptr)));
    EXPECT_TRUE(Point::equal(expected, allocated));

    Point table[npoints][table_width] = {};
    Point caller_populated = {};
    ASSERT_TRUE((blst::scalar::mult_window_n<Point, window_size>(
        caller_populated,
        compact_points,
        npoints,
        compact_scalars,
        8,
        table)));
    EXPECT_TRUE(Point::equal(expected, caller_populated));

    Point precomputed = {};
    ASSERT_TRUE((blst::scalar::mult_window_n<Point, window_size>(
        precomputed,
        nullptr,
        npoints,
        compact_scalars,
        8,
        table)));
    EXPECT_TRUE(Point::equal(expected, precomputed));
}

template<typename Point>
void expect_rejected(bool_t status, const Point& before, const Point& after)
{
    EXPECT_FALSE(status);
    EXPECT_TRUE(Point::equal(before, after));
}

struct FailingAllocator
{
    BLST_HD static void* allocate(size_t)
    {
        return nullptr;
    }

    BLST_HD static void deallocate(void*)
    {
    }
};

TEST(ScalarMultWindow, G1SupportsAllocatedAndCallerTables)
{
    expect_all_mult_window_modes<G1Point>();
}

TEST(ScalarMultWindow, G2SupportsAllocatedAndCallerTables)
{
    expect_all_mult_window_modes<G2Point>();
}

TEST(ScalarMultWindow, RejectsInvalidInputsWithoutChangingOutput)
{
    constexpr int window_size = 5;
    constexpr size_t table_width = size_t{1} << (window_size - 1);

    const G1Point point = G1Point::generator();
    const G1Point* points[] = {&point};
    const G1Point* null_first_point[] = {nullptr};
    const byte scalar = 3;
    const byte* scalars[] = {&scalar};
    const byte* null_first_scalar[] = {nullptr};
    G1Point table[1][table_width] = {};
    const G1Point before = generator_multiple<G1Point>(9);
    G1Point after = before;

    expect_rejected(
        blst::scalar::mult_window_n<G1Point, window_size>(
            after, points, 0, scalars, 8, table),
        before,
        after);

    after = before;
    expect_rejected(
        blst::scalar::mult_window_n<G1Point, window_size>(
            after, points, 1, scalars, 0, table),
        before,
        after);

    after = before;
    expect_rejected(
        blst::scalar::mult_window_n<G1Point, window_size>(
            after, points, 1, nullptr, 8, table),
        before,
        after);

    after = before;
    expect_rejected(
        blst::scalar::mult_window_n<G1Point, window_size>(
            after, points, 1, null_first_scalar, 8, table),
        before,
        after);

    after = before;
    expect_rejected(
        blst::scalar::mult_window_n<G1Point, window_size>(
            after, nullptr, 1, scalars, 8, nullptr),
        before,
        after);

    after = before;
    expect_rejected(
        blst::scalar::mult_window_n<G1Point, window_size>(
            after, null_first_point, 1, scalars, 8, table),
        before,
        after);

    after = before;
    expect_rejected(
        blst::scalar::mult_window_n<G1Point, window_size>(
            after, points, 1, scalars, (size_t)-1, table),
        before,
        after);
}

TEST(ScalarMultWindow, RejectsAllocationOverflowAndFailure)
{
    constexpr int window_size = 5;
    constexpr size_t table_width = size_t{1} << (window_size - 1);

    const G1Point point = G1Point::generator();
    const G1Point* points[] = {&point};
    const byte scalar = 3;
    const byte* scalars[] = {&scalar};
    const G1Point before = generator_multiple<G1Point>(9);
    G1Point after = before;

    const size_t table_row_size = sizeof(G1Point[table_width]);
    const size_t overflowing_npoints = (size_t)-1 / table_row_size + 1;
    expect_rejected(
        blst::scalar::mult_window_n<G1Point, window_size>(
            after, points, overflowing_npoints, scalars, 8, nullptr),
        before,
        after);

    after = before;
    expect_rejected(
        blst::scalar::mult_window_n<G1Point, window_size, FailingAllocator>(
            after, points, 1, scalars, 8, nullptr),
        before,
        after);
}

}  // namespace

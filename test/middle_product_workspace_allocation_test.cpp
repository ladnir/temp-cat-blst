// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs

#include "blst/config.h"
#include "blst/fft/middle_product.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <new>
#include <span>
#include <type_traits>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

bool fail_allocation = false;

void* allocate_unaligned(std::size_t size)
{
    if (fail_allocation) {
        throw std::bad_alloc();
    }
    size = std::max(size, std::size_t{1});
    if (void* const allocation = std::malloc(size)) {
        return allocation;
    }
    throw std::bad_alloc();
}

void* allocate_aligned(std::size_t size, const std::size_t alignment)
{
    if (fail_allocation) {
        throw std::bad_alloc();
    }
    size = std::max(size, std::size_t{1});
#if defined(_MSC_VER)
    if (void* const allocation = _aligned_malloc(size, alignment)) {
        return allocation;
    }
#else
    void* allocation = nullptr;
    if (posix_memalign(&allocation, alignment, size) == 0) {
        return allocation;
    }
#endif
    throw std::bad_alloc();
}

void free_aligned(void* const allocation) noexcept
{
#if defined(_MSC_VER)
    _aligned_free(allocation);
#else
    std::free(allocation);
#endif
}

template<class Point>
Point point_multiple(const std::uint64_t scalar)
{
    blst::Fr multiplier = {};
    blst::Fr::from_uint64(multiplier, scalar);
    blst::Fr non_montgomery = {};
    blst::Fr::from_mont(non_montgomery, multiplier);
    Point result = {};
    Point::mult(
        result,
        Point::generator(),
        reinterpret_cast<const blst::byte*>(&non_montgomery),
        256);
    return result;
}

template<class Context>
bool run_native_probe(
    const std::span<const typename Context::G1Point> a,
    const std::span<const typename Context::G2Point> b)
{
    const std::size_t n = a.size();
    const auto precomputed =
        blst::detail::precompute_middle_product_b_impl<Context>(
            b, n, nullptr);
    blst::detail::MiddleProductWorkspaceImpl<Context> workspace(
        precomputed.fft_size);

    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        std::vector<blst::GT> expected(n, blst::GT::identity());
        std::vector<blst::GT> actual(n, blst::GT::identity());
        std::vector<blst::GT> ranged_actual(n, blst::GT::identity());
        blst::detail::middle_product_impl<Context>(
            a,
            precomputed,
            std::span<blst::GT>(expected),
            workspace,
            nullptr);
        fail_allocation = true;
        try {
            blst::detail::middle_product_impl<Context>(
                a,
                precomputed,
                std::span<blst::GT>(actual),
                workspace,
                nullptr);
            blst::detail::middle_product_ranged_serial_impl<Context>(
                a,
                precomputed,
                std::span<blst::GT>(ranged_actual),
                workspace);
        } catch (...) {
            fail_allocation = false;
            return false;
        }
        fail_allocation = false;
        for (std::size_t i = 0; i < n; ++i) {
            if (!blst::GT::equal(expected[i], actual[i]) ||
                !blst::GT::equal(expected[i], ranged_actual[i])) {
                return false;
            }
        }
        return true;
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    } else {
        std::vector<blst::avx::GT> expected(n);
        std::vector<blst::avx::GT> actual(n);
        std::vector<blst::avx::GT> ranged_actual(n);
        blst::detail::middle_product_impl<Context>(
            a,
            precomputed,
            std::span<blst::avx::GT>(expected),
            workspace,
            nullptr);
        fail_allocation = true;
        try {
            blst::detail::middle_product_impl<Context>(
                a,
                precomputed,
                std::span<blst::avx::GT>(actual),
                workspace,
                nullptr);
            blst::detail::middle_product_ranged_serial_impl<Context>(
                a,
                precomputed,
                std::span<blst::avx::GT>(ranged_actual),
                workspace);
        } catch (...) {
            fail_allocation = false;
            return false;
        }
        fail_allocation = false;
        for (std::size_t i = 0; i < n; ++i) {
            if (!blst::avx::GT::equal(expected[i], actual[i]) ||
                !blst::avx::GT::equal(expected[i], ranged_actual[i])) {
                return false;
            }
        }
        return true;
#endif
    }
}

} // namespace

void* operator new(const std::size_t size)
{
    return allocate_unaligned(size);
}

void* operator new[](const std::size_t size)
{
    return allocate_unaligned(size);
}

void operator delete(void* const allocation) noexcept
{
    std::free(allocation);
}

void operator delete[](void* const allocation) noexcept
{
    std::free(allocation);
}

void operator delete(void* const allocation, std::size_t) noexcept
{
    std::free(allocation);
}

void operator delete[](void* const allocation, std::size_t) noexcept
{
    std::free(allocation);
}

void* operator new(const std::size_t size, const std::align_val_t alignment)
{
    return allocate_aligned(size, static_cast<std::size_t>(alignment));
}

void* operator new[](const std::size_t size, const std::align_val_t alignment)
{
    return allocate_aligned(size, static_cast<std::size_t>(alignment));
}

void operator delete(
    void* const allocation,
    const std::align_val_t) noexcept
{
    free_aligned(allocation);
}

void operator delete[](
    void* const allocation,
    const std::align_val_t) noexcept
{
    free_aligned(allocation);
}

void operator delete(
    void* const allocation,
    std::size_t,
    const std::align_val_t) noexcept
{
    free_aligned(allocation);
}

void operator delete[](
    void* const allocation,
    std::size_t,
    const std::align_val_t) noexcept
{
    free_aligned(allocation);
}

int main()
{
    constexpr std::size_t n = 17;
    std::vector<blst::scalar::G1Point> scalar_a;
    std::vector<blst::scalar::G2Point> scalar_b;
    scalar_a.reserve(n);
    scalar_b.reserve(2 * n);
    for (std::size_t i = 0; i < n; ++i) {
        scalar_a.push_back(point_multiple<blst::scalar::G1Point>(i + 1));
    }
    for (std::size_t i = 0; i < 2 * n; ++i) {
        scalar_b.push_back(point_multiple<blst::scalar::G2Point>(i + 3));
    }
    if (!run_native_probe<blst::MiddleProductScalarContext>(
            scalar_a, scalar_b)) {
        return 1;
    }

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    std::vector<blst::avx::G1Point> avx_a;
    std::vector<blst::avx::G2Point> avx_b;
    avx_a.reserve(n);
    avx_b.reserve(2 * n);
    for (const blst::scalar::G1Point& point : scalar_a) {
        avx_a.push_back(blst::avx::avx_map_to_G_detail::from_scalar(point));
    }
    for (const blst::scalar::G2Point& point : scalar_b) {
        avx_b.push_back(blst::avx::avx_map_to_G_detail::from_scalar(point));
    }
    if (!run_native_probe<blst::MiddleProductAvxContext>(avx_a, avx_b)) {
        return 2;
    }
#endif

    return 0;
}

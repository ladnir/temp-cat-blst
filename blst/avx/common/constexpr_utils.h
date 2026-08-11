#pragma once
#include "../../common/defines.h"
#include <bit>
#include <cstdint>
#include <type_traits>

#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#endif

namespace blst::avx {

#if defined(_MSC_VER) && !defined(__clang__)

struct uint128_t {
    uint64_t lo;
    uint64_t hi;

    constexpr uint128_t() : lo(0), hi(0) {}
    constexpr uint128_t(uint64_t value) : lo(value), hi(0) {}
    constexpr uint128_t(uint64_t high, uint64_t low) : lo(low), hi(high) {}

    constexpr operator uint64_t() const { return lo; }
};

constexpr BLST_FORCEINLINE uint128_t operator+(uint128_t a, uint128_t b) {
    const uint64_t low = a.lo + b.lo;
    const uint64_t carry = low < a.lo ? 1u : 0u;
    return uint128_t(a.hi + b.hi + carry, low);
}

constexpr BLST_FORCEINLINE uint128_t& operator+=(uint128_t& a, uint128_t b) {
    a = a + b;
    return a;
}

constexpr BLST_FORCEINLINE uint128_t operator|(uint128_t a, uint128_t b) {
    return uint128_t(a.hi | b.hi, a.lo | b.lo);
}

constexpr BLST_FORCEINLINE uint128_t& operator|=(uint128_t& a, uint128_t b) {
    a = a | b;
    return a;
}

constexpr BLST_FORCEINLINE uint128_t operator&(uint128_t a, uint64_t mask) {
    return uint128_t(0, a.lo & mask);
}

constexpr BLST_FORCEINLINE uint128_t operator<<(uint128_t value, int shift) {
    if (shift <= 0) {
        return value;
    }
    if (shift >= 128) {
        return uint128_t();
    }
    if (shift >= 64) {
        return uint128_t(value.lo << (shift - 64), 0);
    }
    return uint128_t(
        (value.hi << shift) | (value.lo >> (64 - shift)),
        value.lo << shift
    );
}

constexpr BLST_FORCEINLINE uint128_t operator>>(uint128_t value, int shift) {
    if (shift <= 0) {
        return value;
    }
    if (shift >= 128) {
        return uint128_t();
    }
    if (shift >= 64) {
        return uint128_t(0, value.hi >> (shift - 64));
    }
    return uint128_t(
        value.hi >> shift,
        (value.lo >> shift) | (value.hi << (64 - shift))
    );
}

constexpr BLST_FORCEINLINE uint128_t& operator>>=(uint128_t& value, int shift) {
    value = value >> shift;
    return value;
}

constexpr BLST_FORCEINLINE uint128_t mul_u64_u64(uint64_t a, uint64_t b) {
    if (!std::is_constant_evaluated()) {
        uint64_t high = 0;
        const uint64_t low = _umul128(a, b, &high);
        return uint128_t(high, low);
    }

    const uint64_t a0 = a & 0xffffffffULL;
    const uint64_t a1 = a >> 32;
    const uint64_t b0 = b & 0xffffffffULL;
    const uint64_t b1 = b >> 32;

    const uint64_t p0 = a0 * b0;
    const uint64_t p1 = a0 * b1;
    const uint64_t p2 = a1 * b0;
    const uint64_t p3 = a1 * b1;

    const uint64_t carry = (p0 >> 32) + (p1 & 0xffffffffULL) + (p2 & 0xffffffffULL);
    const uint64_t high = p3 + (p1 >> 32) + (p2 >> 32) + (carry >> 32);
    const uint64_t low = (carry << 32) | (p0 & 0xffffffffULL);
    return uint128_t(high, low);
}

constexpr BLST_FORCEINLINE uint128_t mul_u128_u64(uint128_t a, uint64_t b) {
    const uint128_t low = mul_u64_u64(a.lo, b);
    const uint128_t high = mul_u64_u64(a.hi, b);
    return uint128_t(low.hi + high.lo, low.lo);
}

constexpr BLST_FORCEINLINE uint128_t operator*(uint128_t a, uint64_t b) {
    return mul_u128_u64(a, b);
}

constexpr BLST_FORCEINLINE uint128_t operator*(uint64_t a, uint128_t b) {
    return mul_u128_u64(b, a);
}

constexpr BLST_FORCEINLINE uint128_t operator*(uint128_t a, uint128_t b) {
    if (a.hi == 0 && b.hi == 0) {
        return mul_u64_u64(a.lo, b.lo);
    }
    if (a.hi == 0) {
        return mul_u128_u64(b, a.lo);
    }
    if (b.hi == 0) {
        return mul_u128_u64(a, b.lo);
    }

    const uint128_t low = mul_u64_u64(a.lo, b.lo);
    const uint128_t cross1 = mul_u64_u64(a.hi, b.lo);
    const uint128_t cross2 = mul_u64_u64(a.lo, b.hi);
    return uint128_t(low.hi + cross1.lo + cross2.lo, low.lo);
}

// MSVC does not provide a native signed 128-bit integer. Keep a complete
// two's-complement carrier for bound arithmetic and Debug reconstruction
// rather than silently narrowing those calculations to int64_t.
struct int128_t {
    uint64_t lo;
    uint64_t hi;

    constexpr int128_t() : lo(0), hi(0) {}
    constexpr int128_t(int64_t value)
        : lo(static_cast<uint64_t>(value))
        , hi(value < 0 ? UINT64_MAX : 0)
    {}
    constexpr int128_t(uint64_t value) : lo(value), hi(0) {}
    constexpr int128_t(uint64_t high, uint64_t low) : lo(low), hi(high) {}

    constexpr explicit operator int64_t() const
    {
        return std::bit_cast<int64_t>(lo);
    }
};

static_assert(sizeof(uint128_t) == 16);
static_assert(sizeof(int128_t) == 16);

constexpr BLST_FORCEINLINE bool operator==(int128_t a, int128_t b)
{
    return a.hi == b.hi && a.lo == b.lo;
}

constexpr BLST_FORCEINLINE bool operator<(int128_t a, int128_t b)
{
    const bool a_negative = (a.hi >> 63) != 0;
    const bool b_negative = (b.hi >> 63) != 0;
    if (a_negative != b_negative) {
        return a_negative;
    }
    return a.hi != b.hi ? a.hi < b.hi : a.lo < b.lo;
}

constexpr BLST_FORCEINLINE bool operator>(int128_t a, int128_t b)
{
    return b < a;
}

constexpr BLST_FORCEINLINE bool operator<=(int128_t a, int128_t b)
{
    return !(b < a);
}

constexpr BLST_FORCEINLINE bool operator>=(int128_t a, int128_t b)
{
    return !(a < b);
}

constexpr BLST_FORCEINLINE int128_t operator+(int128_t a, int128_t b)
{
    const uint64_t low = a.lo + b.lo;
    const uint64_t carry = low < a.lo ? 1u : 0u;
    return int128_t(a.hi + b.hi + carry, low);
}

constexpr BLST_FORCEINLINE int128_t operator-(int128_t value)
{
    const uint64_t low = ~value.lo + 1;
    return int128_t(~value.hi + (low == 0 ? 1u : 0u), low);
}

constexpr BLST_FORCEINLINE int128_t operator-(int128_t a, int128_t b)
{
    return a + (-b);
}

constexpr BLST_FORCEINLINE int128_t operator<<(int128_t value, int shift)
{
    if (shift <= 0) {
        return value;
    }
    if (shift >= 128) {
        return int128_t();
    }
    if (shift >= 64) {
        return int128_t(value.lo << (shift - 64), 0);
    }
    return int128_t(
        (value.hi << shift) | (value.lo >> (64 - shift)),
        value.lo << shift);
}

constexpr BLST_FORCEINLINE int128_t operator*(int128_t a, int128_t b)
{
    // The low 128 bits of a two's-complement product are the unsigned low
    // product plus both cross products modulo 2^128.
    const uint128_t low = mul_u64_u64(a.lo, b.lo);
    return int128_t(
        low.hi + a.hi * b.lo + a.lo * b.hi,
        low.lo);
}

constexpr uint64_t add_mod_u64(uint64_t x, uint64_t y, uint64_t modulus) {
    return x >= modulus - y ? x - (modulus - y) : x + y;
}

constexpr uint64_t operator%(uint128_t value, uint64_t modulus) {
    BLST_ASSERT(modulus != 0);
    if (modulus == 0) {
        return 0;
    }

    uint64_t result = 0;
    for (int bit = 127; bit >= 0; --bit) {
        result = add_mod_u64(result, result, modulus);
        const uint64_t next =
            bit >= 64 ? ((value.hi >> (bit - 64)) & 1ULL) : ((value.lo >> bit) & 1ULL);
        result = add_mod_u64(result, next, modulus);
    }
    return result;
}
#else
typedef __int128_t int128_t;
typedef unsigned __int128 uint128_t;
#endif

// Zero is the fail-closed Release result for an invalid modulus.
constexpr uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m) {
    BLST_ASSERT(m != 0);
    if (m == 0) {
        return 0;
    }
    return ((uint128_t)a * (uint128_t)b) % m;
}

}  // namespace blst::avx

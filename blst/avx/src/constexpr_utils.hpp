#pragma once
#include <bit>
#include <cstdint>
#include <array>

#if defined(_MSC_VER) && !defined(__clang__)
using int128_t = int64_t;

struct uint128_t {
    uint64_t lo;
    uint64_t hi;

    constexpr uint128_t() : lo(0), hi(0) {}
    constexpr uint128_t(uint64_t value) : lo(value), hi(0) {}
    constexpr uint128_t(uint64_t high, uint64_t low) : lo(low), hi(high) {}

    constexpr operator uint64_t() const { return lo; }
};

constexpr uint128_t operator+(uint128_t a, uint128_t b) {
    const uint64_t low = a.lo + b.lo;
    const uint64_t carry = low < a.lo ? 1u : 0u;
    return uint128_t(a.hi + b.hi + carry, low);
}

constexpr uint128_t& operator+=(uint128_t& a, uint128_t b) {
    a = a + b;
    return a;
}

constexpr uint128_t operator|(uint128_t a, uint128_t b) {
    return uint128_t(a.hi | b.hi, a.lo | b.lo);
}

constexpr uint128_t& operator|=(uint128_t& a, uint128_t b) {
    a = a | b;
    return a;
}

constexpr uint128_t operator&(uint128_t a, uint64_t mask) {
    return uint128_t(0, a.lo & mask);
}

constexpr uint128_t operator<<(uint128_t value, int shift) {
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

constexpr uint128_t operator>>(uint128_t value, int shift) {
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

constexpr uint128_t& operator>>=(uint128_t& value, int shift) {
    value = value >> shift;
    return value;
}

constexpr uint128_t mul_u64_u64(uint64_t a, uint64_t b) {
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

constexpr uint128_t mul_u128_u64(uint128_t a, uint64_t b) {
    const uint128_t low = mul_u64_u64(a.lo, b);
    const uint128_t high = mul_u64_u64(a.hi, b);
    return uint128_t(low.hi + high.lo, low.lo);
}

constexpr uint128_t operator*(uint128_t a, uint64_t b) {
    return mul_u128_u64(a, b);
}

constexpr uint128_t operator*(uint64_t a, uint128_t b) {
    return mul_u128_u64(b, a);
}

constexpr uint128_t operator*(uint128_t a, uint128_t b) {
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

constexpr uint64_t add_mod_u64(uint64_t x, uint64_t y, uint64_t modulus) {
    return x >= modulus - y ? x - (modulus - y) : x + y;
}

constexpr uint64_t operator%(uint128_t value, uint64_t modulus) {
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

constexpr int64_t ceil_div(int64_t a, int64_t b) {
    return (a + b - 1) / b;
}

constexpr int64_t get_mask(int64_t shift) {
    return (1ull << shift) - 1;
}

constexpr int ceil_log2(int x) {
    if (x == 0) return 0;
    if (x == 1) return 0;
    return static_cast<int>(std::bit_width(static_cast<unsigned int>(x - 1)));
}


constexpr uint64_t mod_inverse(uint64_t a, uint64_t m) {
    uint64_t m0 = m, t = 0, q = 0;
    int64_t x0 = 0, x1 = 1;  // Use signed for intermediate calculations
    a = a % m;  // Reduce a modulo m first
    if (a == 0) return 0;  // No inverse if a is 0 mod m
    while (a > 1) {
        if (m == 0) return 0;  // No inverse exists (gcd != 1)
        q = a / m;
        t = m;
        m = a % m;
        a = t;
        int64_t t_x = x0;
        x0 = x1 - static_cast<int64_t>(q) * x0;
        x1 = t_x;
    }
    // Make sure result is positive
    if (x1 < 0) x1 += static_cast<int64_t>(m0);
    return static_cast<uint64_t>(x1);
}

constexpr uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m) {
    return ((uint128_t)a * (uint128_t)b) % m;
}

// Helper to convert hex char to value
constexpr uint8_t hex_char_to_value(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return 0;
}

template<int bytes>
constexpr std::array<uint64_t, 6> read_384bit(const char(&lit)[bytes]) {
    static_assert(bytes >= 2, "Must be at least 2 bytes");
    
    // Parse hex string logically (matching blst/fp12.hpp behavior)
    // Hex string is big-endian (most significant digit first), e.g., "0x1234..."
    // We need to parse it and store as little-endian 64-bit limbs
    // Strategy: Parse 6 limbs of 16 hex digits each, reversing order
    
    // Count hex digits (skip "0x")
    int hex_digits = bytes - 2;
    if (hex_digits > 96) hex_digits = 96; // Max 384 bits = 96 hex digits
    
    // Pad hex digits to 96 (add leading zeros)
    std::array<char, 96> padded_hex = {0};
    int padding = 96 - hex_digits;
    for (int i = 0; i < padding; i++) {
        padded_hex[i] = '0';
    }
    for (int i = 0; i < hex_digits; i++) {
        padded_hex[padding + i] = lit[2 + i];
    }
    
    // Parse 6 limbs of 16 hex digits each
    // Hex string is big-endian, so first 16 digits are most significant -> go to result[5]
    // Last 16 digits are least significant -> go to result[0]
    std::array<uint64_t, 6> result = {0};
    for (int i = 0; i < 6; i++) {
        int hex_offset = (5 - i) * 16; // Reverse order: limb 0 gets last 16 digits, limb 5 gets first 16 digits
        uint64_t limb_value = 0;
        for (int j = 0; j < 16; j++) {
            char hex_char = padded_hex[hex_offset + j];
            uint8_t hex_value = hex_char_to_value(hex_char);
            limb_value = (limb_value << 4) | hex_value;
        }
        result[i] = limb_value;
    }
    
    return result;
}

template<int element_bits>
constexpr std::array<uint64_t, 8> parse_chunks(const std::array<uint64_t, 6> &chunks) {
    std::array<uint64_t, 8> result;
    uint64_t mask = get_mask(element_bits);
    uint128_t current = uint128_t(chunks[0]) | (uint128_t(chunks[1]) << 64);
    int bits_in_current = 128;
    int next_chunk_index = 2;
    for (int i = 0; i < 8; i++) {
        result[i] = current & mask;
        current >>= element_bits;
        bits_in_current -= element_bits;
        if ((bits_in_current < 64) && (next_chunk_index < 6)) {
            current |= uint128_t(chunks[next_chunk_index]) << bits_in_current;
            next_chunk_index++;
            bits_in_current += 64;
        }
    }
    return result;
}

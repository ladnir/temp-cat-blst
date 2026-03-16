#pragma once
#include <cstdint>
#include <array>

typedef __int128_t int128_t;
typedef unsigned __int128 uint128_t;

constexpr int64_t ceil_div(int64_t a, int64_t b) {
    return (a + b - 1) / b;
}

constexpr int64_t get_mask(int64_t shift) {
    return (1ull << shift) - 1;
}

constexpr int ceil_log2(int x) {
    if (x == 0) return 0;
    if (x == 1) return 0;
    return 32 - __builtin_clz(x - 1);
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
    uint128_t current = chunks[0] | (uint128_t)chunks[1] << 64;
    int bits_in_current = 128;
    int next_chunk_index = 2;
    for (int i = 0; i < 8; i++) {
        result[i] = current & mask;
        current >>= element_bits;
        bits_in_current -= element_bits;
        if ((bits_in_current < 64) && (next_chunk_index < 6)) {
            current |= (uint128_t)chunks[next_chunk_index] << bits_in_current;
            next_chunk_index++;
            bits_in_current += 64;
        }
    }
    return result;
}
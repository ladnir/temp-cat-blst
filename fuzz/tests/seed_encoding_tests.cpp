// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/support/seed_encoding.h"

#include <array>

int main()
{
    using namespace catblst::fuzz;
    constexpr std::array hex{
        std::byte{'h'}, std::byte{'e'}, std::byte{'x'}, std::byte{':'},
        std::byte{'0'}, std::byte{'0'}, std::byte{'f'}, std::byte{'f'},
        std::byte{'\n'}};
    const auto decoded_hex = decode_hex_seed<2>(hex);
    if (!decoded_hex || decoded_hex->size != 2 ||
        decoded_hex->bytes[0] != std::byte{0} ||
        decoded_hex->bytes[1] != std::byte{0xff})
        return 1;

    constexpr std::array base64{
        std::byte{'b'}, std::byte{'6'}, std::byte{'4'}, std::byte{':'},
        std::byte{'A'}, std::byte{'Q'}, std::byte{'I'}, std::byte{'='}};
    const auto decoded_base64 = decode_base64_seed<2>(base64);
    if (!decoded_base64 || decoded_base64->size != 2 ||
        decoded_base64->bytes[0] != std::byte{1} ||
        decoded_base64->bytes[1] != std::byte{2})
        return 2;

    constexpr std::array text{
        std::byte{'t'}, std::byte{'e'}, std::byte{'x'}, std::byte{'t'},
        std::byte{':'}, std::byte{'o'}, std::byte{'k'}, std::byte{'\n'}};
    const auto payload = text_seed_payload(text);
    if (payload.size() != 2 || payload[0] != std::byte{'o'} ||
        payload[1] != std::byte{'k'})
        return 3;

    constexpr std::array malformed{
        std::byte{'h'}, std::byte{'e'}, std::byte{'x'}, std::byte{':'},
        std::byte{'0'}, std::byte{'g'}};
    if (decode_hex_seed<1>(malformed))
        return 4;
    return 0;
}

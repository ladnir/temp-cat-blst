// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter.h"
#include "fuzz/support/driver.h"
#include "fuzz/support/seed_encoding.h"

#include <algorithm>

namespace catblst::fuzz::target
{
    void run(std::span<const std::byte> input)
    {
        adapter::run_bytes_case(input.first(
            std::min<std::size_t>(input.size(), 1024)));
        if (const auto decoded = decode_hex_seed<1024>(input))
            adapter::run_bytes_case(decoded->span());
        const auto text = text_seed_payload(input);
        if (!text.empty())
            adapter::run_bytes_case(text);
    }
}

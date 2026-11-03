// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "test/fuzz/internal/adapter.h"
#include "test/fuzz/support/driver.h"
#include "test/fuzz/support/seed_encoding.h"

#include <algorithm>

namespace catblst::fuzz::target
{
    void run(std::span<const std::byte> input)
    {
        if (const auto decoded = decode_hex_seed<1024>(input))
        {
            adapter::run_bytes_case(decoded->span());
            return;
        }
        if (has_prefix(input, "text:"))
        {
            adapter::run_bytes_case(text_seed_payload(input));
            return;
        }
        adapter::run_bytes_case(input.first(
            std::min<std::size_t>(input.size(), 1024)));
    }
}

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
        adapter::run_decode_case(input.first(
            std::min<std::size_t>(input.size(), 192)));
        if (const auto decoded = decode_hex_seed<192>(input))
            adapter::run_decode_case(decoded->span());
        if (const auto decoded = decode_base64_seed<192>(input))
            adapter::run_decode_case(decoded->span());
    }
}

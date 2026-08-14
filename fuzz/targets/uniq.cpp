// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter.h"
#include "fuzz/support/driver.h"
#include "fuzz/support/seed_encoding.h"

namespace catblst::fuzz::target
{
    void run(std::span<const std::byte> input)
    {
        adapter::run_uniq_case(input);
        if (const auto decoded = decode_hex_seed<4097>(input))
            adapter::run_uniq_case(decoded->span());
    }
}


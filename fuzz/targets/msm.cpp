// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter.h"
#include "fuzz/support/driver.h"
#include "fuzz/support/seed_encoding.h"

namespace catblst::fuzz::target
{
    void run(std::span<const std::byte> input)
    {
        if (const auto decoded = decode_hex_seed<64>(input))
        {
            adapter::run_msm_case(decoded->span());
            return;
        }
        adapter::run_msm_case(input);
    }
}

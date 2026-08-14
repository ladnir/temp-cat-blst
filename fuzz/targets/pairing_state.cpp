// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter.h"
#include "fuzz/support/driver.h"

namespace catblst::fuzz::target
{
    void run(std::span<const std::byte> input)
    {
        adapter::run_pairing_state_case(input);
    }
}


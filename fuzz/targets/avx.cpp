// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter.h"

namespace catblst::fuzz::target
{
    void run(std::span<const std::byte> input)
    {
        adapter::run_avx_case(input);
    }
}

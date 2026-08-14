// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/support/driver.h"

namespace catblst::fuzz::target
{
    void run(std::span<const std::byte> input);
}

extern "C" int LLVMFuzzerTestOneInput(
    const unsigned char* data,
    std::size_t size)
{
    catblst::fuzz::target::run(std::as_bytes(std::span(data, size)));
    return 0;
}


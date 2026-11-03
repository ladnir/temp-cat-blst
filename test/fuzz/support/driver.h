// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <span>

namespace catblst::fuzz
{
    using Driver = void (*)(std::span<const std::byte>);
}


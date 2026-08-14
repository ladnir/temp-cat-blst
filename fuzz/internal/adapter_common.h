// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "bindings/blst_scalar.h"
#include "fuzz/internal/adapter.h"
#include "fuzz/internal/adapter_access.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <source_location>
#include <span>
#include <type_traits>

namespace catblst::fuzz::adapter::detail
{
    inline void require(
        bool condition,
        const std::source_location location =
            std::source_location::current()) noexcept
    {
        if (!condition)
        {
            std::fprintf(
                stderr,
                "fuzz property failed at %s:%u\n",
                location.file_name(),
                location.line());
            std::abort();
        }
    }

    template<class Native, std::size_t Size, class Tag>
    [[nodiscard]] Native native(
        const DomainValue<OpaqueStorage<Size>, Tag>& input) noexcept
    {
        return internal::AdapterAccess::unwrap<Native>(input);
    }

    template<class T>
    [[nodiscard]] bool object_equal(const T& left, const T& right) noexcept
    {
        static_assert(std::is_trivially_copyable_v<T>);
        return std::memcmp(&left, &right, sizeof(T)) == 0;
    }

    template<std::size_t Size>
    [[nodiscard]] std::array<blst_byte, Size> fixed_bytes(
        std::span<const std::byte> input) noexcept
    {
        std::array<blst_byte, Size> result{};
        const std::size_t copied = std::min(input.size(), result.size());
        if (copied != 0)
            std::memcpy(result.data(), input.data(), copied);
        return result;
    }
}

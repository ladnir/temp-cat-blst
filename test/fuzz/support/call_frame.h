// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <span>
#include <type_traits>
#include <vector>

namespace catblst::fuzz
{
    enum class StorageShape : std::uint8_t
    {
        disjoint,
        exact_alias,
        null_zero_length,
        invalid_size,
        short_scratch,
        rejected_overlap,
    };

    struct CallRanges final
    {
        std::span<std::byte> first;
        std::span<std::byte> second;
    };

    /**
     * Owns every byte range supplied to one C call.
     *
     * Builders publish pointers only after all vectors reach their final size.
     * A campaign selects only shapes documented for the target function.
     */
    class CallFrame final
    {
      public:
        [[nodiscard]] static CallFrame disjoint() noexcept
        {
            return CallFrame(StorageShape::disjoint);
        }

        [[nodiscard]] static CallFrame exact_alias() noexcept
        {
            return CallFrame(StorageShape::exact_alias);
        }

        [[nodiscard]] static CallFrame null_zero_length() noexcept
        {
            return CallFrame(StorageShape::null_zero_length);
        }

        [[nodiscard]] static CallFrame invalid_size() noexcept
        {
            return CallFrame(StorageShape::invalid_size);
        }

        [[nodiscard]] static CallFrame short_scratch() noexcept
        {
            return CallFrame(StorageShape::short_scratch);
        }

        [[nodiscard]] static CallFrame rejected_overlap() noexcept
        {
            return CallFrame(StorageShape::rejected_overlap);
        }

        [[nodiscard]] StorageShape shape() const noexcept { return shape_; }

        [[nodiscard]] std::span<std::byte> add_region(std::size_t size)
        {
            regions_.emplace_back(size);
            return regions_.back();
        }

        /**
         * Builds two readable/writable ranges with the selected relationship.
         *
         * Exact aliasing requires equal nonzero sizes. Rejected overlap means
         * a one-byte overlap and therefore requires two nonempty ranges. Short
         * scratch treats the second requested size as the documented minimum.
         */
        [[nodiscard]] CallRanges two_ranges(
            std::size_t first_size,
            std::size_t second_size)
        {
            switch (shape_)
            {
            case StorageShape::disjoint:
            case StorageShape::invalid_size:
                return {add_region(first_size), add_region(second_size)};
            case StorageShape::exact_alias:
                if (first_size == 0 || first_size != second_size)
                    return {};
                {
                    auto storage = add_region(first_size);
                    return {storage, storage};
                }
            case StorageShape::null_zero_length:
                return {};
            case StorageShape::short_scratch:
                return {
                    add_region(first_size),
                    add_region(second_size == 0 ? 0 : second_size - 1)};
            case StorageShape::rejected_overlap:
                if (first_size == 0 || second_size == 0)
                    return {};
                if (first_size >
                    std::numeric_limits<std::size_t>::max() -
                        (second_size - 1))
                    return {};
                {
                    auto storage = add_region(first_size + second_size - 1);
                    return {
                        storage.first(first_size),
                        storage.subspan(first_size - 1, second_size)};
                }
            }
            return {};
        }

        [[nodiscard]] static bool byte_size(
            std::size_t count,
            std::size_t element_size,
            std::size_t& result) noexcept
        {
            if (element_size != 0 &&
                count > std::numeric_limits<std::size_t>::max() / element_size)
            {
                result = 0;
                return false;
            }
            result = count * element_size;
            return true;
        }

      private:
        explicit CallFrame(StorageShape shape) noexcept : shape_(shape) {}

        StorageShape shape_;
        std::deque<std::vector<std::byte>> regions_;
    };
}

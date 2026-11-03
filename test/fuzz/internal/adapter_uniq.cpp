// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "test/fuzz/internal/adapter_common.h"

#include <new>

namespace catblst::fuzz::adapter
{
    void run_uniq_case(std::span<const std::byte> input)
    {
        using detail::require;
        const std::size_t capacity = input.empty()
            ? 0
            : std::to_integer<unsigned char>(input.front()) % 65U;
        const std::size_t size = ::blst_uniq_sizeof(capacity);
        const std::size_t alignment = ::blst_uniq_alignof();
        require(alignment != 0);
        if (capacity == 0)
        {
            require(size == 0);
            require(!::blst_uniq_init(nullptr, 0));
            return;
        }
        require(size != 0);
        void* memory = ::operator new(size, std::align_val_t(alignment));
        const auto cleanup = [&]() noexcept {
            ::operator delete(memory, std::align_val_t(alignment));
        };
        try
        {
            auto* tree = static_cast<blst_uniq*>(memory);
            require(::blst_uniq_init(tree, capacity));
            input = input.subspan(1).first(
                std::min<std::size_t>(input.size() - 1, 4096));
            const blst_byte* message = input.empty()
                ? nullptr
                : static_cast<const blst_byte*>(static_cast<const void*>(input.data()));
            require(::blst_uniq_insert(tree, message, input.size()) ==
                    BLST_UNIQ_INSERTED);
            require(::blst_uniq_insert(tree, message, input.size()) ==
                    BLST_UNIQ_DUPLICATE);

            require(::blst_uniq_init(tree, capacity));
            std::array<std::array<blst_byte, 64>, 64> retained{};
            for (std::size_t index = 0; index < capacity; ++index)
            {
                retained[index][0] = static_cast<blst_byte>(index + 1);
                for (std::size_t byte = 1; byte <= index % 63; ++byte)
                    retained[index][byte] = static_cast<blst_byte>(index + byte);
                require(::blst_uniq_insert(
                            tree, retained[index].data(), (index % 63) + 1) ==
                        BLST_UNIQ_INSERTED);
            }
            const std::array<blst_byte, 2> extra{0xfe, 0xed};
            require(::blst_uniq_insert(tree, extra.data(), extra.size()) ==
                    BLST_UNIQ_CAPACITY_EXCEEDED);
        }
        catch (...)
        {
            cleanup();
            throw;
        }
        cleanup();
    }
}

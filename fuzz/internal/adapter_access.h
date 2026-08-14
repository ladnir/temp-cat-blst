// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "fuzz/support/domain.h"
#include "fuzz/support/storage.h"

#include <cstring>
#include <type_traits>
#include <utility>

namespace catblst::fuzz::internal
{
    /** Representation bridge owned exclusively by the fuzz adapter. */
    struct AdapterAccess final
    {
        template<class Native, std::size_t Size, class Tag>
        [[nodiscard]] static DomainValue<OpaqueStorage<Size>, Tag> wrap(
            const Native& native) noexcept
        {
            static_assert(sizeof(Native) == Size);
            static_assert(std::is_trivially_copyable_v<Native>);
            OpaqueStorage<Size> storage;
            std::memcpy(storage.bytes.data(), &native, Size);
            return DomainValue<OpaqueStorage<Size>, Tag>(std::move(storage));
        }

        template<class Native, std::size_t Size, class Tag>
        [[nodiscard]] static Native unwrap(
            const DomainValue<OpaqueStorage<Size>, Tag>& value) noexcept
        {
            static_assert(sizeof(Native) == Size);
            static_assert(std::is_trivially_copyable_v<Native>);
            Native native;
            std::memcpy(&native, value.storage_.bytes.data(), Size);
            return native;
        }

        template<class Storage, class Tag>
        [[nodiscard]] static const Storage& storage(
            const DomainValue<Storage, Tag>& value) noexcept
        {
            return value.storage_;
        }

        template<class Storage, class Tag>
        [[nodiscard]] static DomainValue<Storage, Tag> make(Storage storage)
        {
            return DomainValue<Storage, Tag>(std::move(storage));
        }

        template<class Storage, class Tag>
        [[nodiscard]] static Storage take(
            DomainValue<Storage, Tag>&& value) noexcept(
                std::is_nothrow_move_constructible_v<Storage>)
        {
            return std::move(value.storage_);
        }
    };
}

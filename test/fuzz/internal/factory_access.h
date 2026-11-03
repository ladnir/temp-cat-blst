// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "test/fuzz/support/domain.h"

#include <utility>

namespace catblst::fuzz::internal
{
    /** Private construction and observation authority for reviewed factories. */
    struct FactoryAccess final
    {
        template<class Storage, class Tag>
        [[nodiscard]] static DomainValue<Storage, Tag> make(Storage storage)
        {
            return DomainValue<Storage, Tag>(std::move(storage));
        }

        template<class Storage, class Tag>
        [[nodiscard]] static const Storage& storage(
            const DomainValue<Storage, Tag>& value) noexcept
        {
            return value.storage_;
        }
    };
}

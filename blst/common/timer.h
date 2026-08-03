// Copyright Supranational LLC
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <chrono>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>

namespace blst {

class TimingStats {
public:
    TimingStats() noexcept = default;
    ~TimingStats() noexcept = default;

    TimingStats(const TimingStats&) = delete;
    TimingStats& operator=(const TimingStats&) = delete;
    TimingStats(TimingStats&&) = delete;
    TimingStats& operator=(TimingStats&&) = delete;

    void add(std::string_view key, std::uint64_t ns) noexcept
    {
        try {
            if (!totals_)
                totals_.emplace();
            totals_->try_emplace(std::string(key), 0).first->second += ns;
        } catch (...) {
            // Timing is diagnostic-only and must never affect product work.
        }
    }

    [[nodiscard]] std::uint64_t get_ns(std::string_view key) const noexcept
    {
        if (!totals_)
            return 0;

        try {
            const auto it = totals_->find(key);
            if (it != totals_->end())
                return it->second;
        } catch (...) {
            // A failed diagnostic lookup is indistinguishable from no sample.
        }
        return 0;
    }

    void clear() noexcept
    {
        totals_.reset();
    }

private:
    struct StringHash {
        using is_transparent = void;

        std::size_t operator()(std::string_view value) const noexcept
        {
            return std::hash<std::string_view>{}(value);
        }
    };

    struct StringEqual {
        using is_transparent = void;

        bool operator()(std::string_view lhs, std::string_view rhs) const noexcept
        {
            return lhs == rhs;
        }
    };

    using Totals = std::unordered_map<
        std::string, std::uint64_t, StringHash, StringEqual>;
    std::optional<Totals> totals_;
};

class ScopedTimer {
public:
    ScopedTimer(TimingStats* stats, std::string_view key) noexcept
        : stats_(stats), key_(key),
          start_(stats == nullptr ? Clock::time_point{} : Clock::now())
    {
    }

    ~ScopedTimer() noexcept
    {
        if (stats_ == nullptr)
            return;

        const auto end = Clock::now();
        const auto ns = static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_).count());
        stats_->add(key_, ns);
    }

    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;
    ScopedTimer(ScopedTimer&&) = delete;
    ScopedTimer& operator=(ScopedTimer&&) = delete;

private:
    using Clock = std::chrono::steady_clock;

    TimingStats* stats_;
    std::string_view key_;
    Clock::time_point start_;
};

} // namespace blst

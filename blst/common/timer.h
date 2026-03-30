#pragma once

#include <chrono>
#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>

namespace blst {

class TimingStats {
public:
    void add(std::string_view key, std::uint64_t ns)
    {
        totals_[std::string(key)] += ns;
    }

    std::uint64_t get_ns(std::string_view key) const
    {
        auto it = totals_.find(std::string(key));
        if (it == totals_.end()) {
            return 0;
        }
        return it->second;
    }

    void clear()
    {
        totals_.clear();
    }

private:
    std::unordered_map<std::string, std::uint64_t> totals_;
};

class ScopedTimer {
public:
    ScopedTimer(TimingStats* stats, std::string_view key)
        : stats_(stats), key_(key), start_(Clock::now())
    {
    }

    ~ScopedTimer()
    {
        if (stats_ == nullptr) {
            return;
        }
        auto end = Clock::now();
        auto ns = static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_).count());
        stats_->add(key_, ns);
    }

private:
    using Clock = std::chrono::steady_clock;

    TimingStats* stats_;
    std::string_view key_;
    Clock::time_point start_;
};

} // namespace blst

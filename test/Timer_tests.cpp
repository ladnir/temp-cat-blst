// Copyright 2026 Peter Rindal.
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <string_view>
#include <type_traits>
#include <utility>

#include "blst/common/timer.h"

static_assert(noexcept(
    std::declval<blst::TimingStats&>().add("sample", 1)));
static_assert(noexcept(
    std::declval<blst::TimingStats&>().add("sample", ".suffix", 1)));
static_assert(noexcept(
    std::declval<const blst::TimingStats&>().get_ns("sample")));
static_assert(noexcept(std::declval<blst::TimingStats&>().clear()));
static_assert(std::is_nothrow_default_constructible_v<blst::TimingStats>);
static_assert(std::is_nothrow_destructible_v<blst::TimingStats>);
static_assert(std::is_nothrow_constructible_v<
    blst::ScopedTimer, blst::TimingStats*, std::string_view>);
static_assert(std::is_nothrow_destructible_v<blst::ScopedTimer>);

TEST(Timer, AccumulatesAndClearsWithoutARequiredSink)
{
    blst::TimingStats timings;
    timings.add("sample", 4);
    timings.add("sample", 6);
    timings.add("sample", ".suffix", 7);
    EXPECT_EQ(timings.get_ns("sample"), 10);
    EXPECT_EQ(timings.get_ns("sample.suffix"), 7);
    EXPECT_EQ(timings.get_ns("missing"), 0);

    {
        blst::ScopedTimer disabled(nullptr, "disabled");
    }

    timings.clear();
    EXPECT_EQ(timings.get_ns("sample"), 0);
}

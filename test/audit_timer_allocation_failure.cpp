// Copyright 2026 Peter Rindal.
// SPDX-License-Identifier: Apache-2.0

#include <cstdlib>
#include <new>

#include "blst/common/timer.h"

namespace {

bool fail_allocation = false;

}

void* operator new(std::size_t size)
{
    if (fail_allocation)
        throw std::bad_alloc();
    if (void* const allocation = std::malloc(size))
        return allocation;
    throw std::bad_alloc();
}

void operator delete(void* allocation) noexcept
{
    std::free(allocation);
}

void operator delete(void* allocation, std::size_t) noexcept
{
    std::free(allocation);
}

int main()
{
    constexpr char key[] =
        "timer-allocation-failure-probe-key-longer-than-small-string-storage";
    blst::TimingStats timings;
    timings.add(key, 7);

    fail_allocation = true;
    timings.add(key, 9);
    return timings.get_ns(key) == 7 ? 0 : 1;
}

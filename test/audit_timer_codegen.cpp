// Copyright 2026 Peter Rindal.
// SPDX-License-Identifier: Apache-2.0

#include "blst/common/timer.h"

#if defined(_MSC_VER)
# define AUDIT_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define AUDIT_NOINLINE __attribute__((noinline))
#else
# define AUDIT_NOINLINE
#endif

extern "C" AUDIT_NOINLINE void audit_disabled_scoped_timer() noexcept
{
    blst::ScopedTimer timer(nullptr, "disabled");
}

#undef AUDIT_NOINLINE

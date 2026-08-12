/*
 * Copyright Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

// Compatibility shim for the historical misspelling. New CatBLST code must
// include primitive_root.h. Frozen BTx/BTE callers retain this path unchanged.
#include "blst/fft/primitive_root.h"

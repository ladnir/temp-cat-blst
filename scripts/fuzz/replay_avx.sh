#!/usr/bin/env bash
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 BUILD_DIR" >&2
    exit 2
fi

cmake --build "$1" --target blst_fuzz_avx --parallel 2
ctest --test-dir "$1" -R blst_fuzz_avx_replay --output-on-failure -j1

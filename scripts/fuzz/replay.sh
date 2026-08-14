#!/usr/bin/env bash
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 BUILD_DIR" >&2
    exit 2
fi

build_dir=$1
cmake --build "$build_dir" \
    --target blst_fuzz_phase0 blst_fuzz_phase1 blst_fuzz_phase2 blst_fuzz_phase3 \
    --parallel 2
ctest --test-dir "$build_dir" -L fuzz --output-on-failure -j1

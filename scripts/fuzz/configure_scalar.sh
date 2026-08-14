#!/usr/bin/env bash
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
    echo "usage: $0 BUILD_DIR [Debug|RelWithDebInfo] [coverage]" >&2
    exit 2
fi

build_dir=$1
build_type=${2:-RelWithDebInfo}
coverage=${3:-off}
coverage_flag=OFF
if [[ "$coverage" == coverage ]]; then
    coverage_flag=ON
elif [[ "$coverage" != off ]]; then
    echo "coverage mode must be 'coverage' or 'off'" >&2
    exit 2
fi

cmake -S . -B "$build_dir" -G Ninja \
    -DCMAKE_BUILD_TYPE="$build_type" \
    -DCMAKE_C_COMPILER=clang \
    -DCMAKE_CXX_COMPILER=clang++ \
    -DBLST_FUZZ=ON \
    -DBLST_FUZZ_AVX=OFF \
    -DBLST_FUZZ_COVERAGE="$coverage_flag" \
    -DBLST_ASAN=ON \
    -DBLST_UBSAN=ON \
    -DBLST_ASM=OFF \
    -DBLST_ADX=OFF \
    -DBLST_DEFAULT_ASM=OFF \
    -DBLST_AVX512=OFF \
    -DBLST_AVX512_PORTABLE=OFF \
    -DBLST_DEFAULT_AVX=OFF \
    -DBLST_CUDA=OFF \
    -DBLST_CUDA_MSM=OFF \
    -DBLST_TEST=OFF \
    -DBLST_BENCHMARK=OFF

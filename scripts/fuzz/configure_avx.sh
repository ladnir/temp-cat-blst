#!/usr/bin/env bash
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

if [[ $# -ne 2 || ( $2 != portable && $2 != native ) ]]; then
    echo "usage: $0 BUILD_DIR portable|native" >&2
    exit 2
fi

build_dir=$1
mode=$2
native=OFF
portable=OFF
if [[ "$mode" == native ]]; then
    if [[ ${BLST_NATIVE_AVX_FUZZ_APPROVED:-} != YES ]]; then
        echo "set BLST_NATIVE_AVX_FUZZ_APPROVED=YES before native AVX execution" >&2
        exit 2
    fi
    native=ON
else
    portable=ON
fi

cmake -S . -B "$build_dir" -G Ninja \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_C_COMPILER=clang \
    -DCMAKE_CXX_COMPILER=clang++ \
    -DBLST_FUZZ=ON \
    -DBLST_FUZZ_AVX=ON \
    -DBLST_ASAN=ON \
    -DBLST_UBSAN=ON \
    -DBLST_ASM=OFF \
    -DBLST_ADX=OFF \
    -DBLST_DEFAULT_ASM=OFF \
    -DBLST_AVX512="$native" \
    -DBLST_AVX512_PORTABLE="$portable" \
    -DBLST_DEFAULT_AVX=ON \
    -DBLST_CUDA=OFF \
    -DBLST_TEST=OFF \
    -DBLST_BENCHMARK=OFF

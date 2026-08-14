#!/usr/bin/env bash
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 BUILD_DIR" >&2
    exit 2
fi

build_dir=$1
coverage_dir="$build_dir/coverage"
mkdir -p "$coverage_dir"

targets=(
    blst_fuzz_decode
    blst_fuzz_bytes
    blst_fuzz_fields
    blst_fuzz_groups
    blst_fuzz_protocol
    blst_fuzz_pairing_state
    blst_fuzz_uniq
    blst_fuzz_fft_middle
    blst_fuzz_msm
)

cmake --build "$build_dir" \
    --target blst_fuzz_phase1 blst_fuzz_phase2 blst_fuzz_phase3 \
    --parallel 2

profiles=()
objects=()
for target in "${targets[@]}"; do
    profile="$coverage_dir/$target.profraw"
    LLVM_PROFILE_FILE="$profile" "$build_dir/fuzz/$target" \
        -runs=0 \
        -rss_limit_mb=2048 \
        "fuzz/corpus/$target"
    profiles+=("$profile")
    if [[ "$target" != "${targets[0]}" ]]; then
        objects+=("-object" "$build_dir/fuzz/$target")
    fi
done

llvm-profdata merge -sparse "${profiles[@]}" -o "$coverage_dir/corpus.profdata"
llvm-cov report \
    "$build_dir/fuzz/${targets[0]}" \
    "${objects[@]}" \
    -instr-profile="$coverage_dir/corpus.profdata" \
    -ignore-filename-regex='(fuzz/|test/|_deps/)' \
    > "$coverage_dir/report.txt"
llvm-cov show \
    "$build_dir/fuzz/${targets[0]}" \
    "${objects[@]}" \
    -instr-profile="$coverage_dir/corpus.profdata" \
    -ignore-filename-regex='(fuzz/|test/|_deps/)' \
    -format=html \
    -output-dir="$coverage_dir/html"

echo "coverage report: $coverage_dir/report.txt"

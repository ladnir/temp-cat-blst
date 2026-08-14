#!/usr/bin/env bash
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/target_limits.sh"

if [[ $# -ne 1 ]]; then
    echo "usage: $0 BUILD_DIR" >&2
    exit 2
fi

build_dir=$1
coverage_dir="$build_dir/coverage"
mkdir -p "$coverage_dir"
rm -f "$coverage_dir"/*.profraw "$coverage_dir"/*.profdata \
    "$coverage_dir/report.txt"
rm -rf "$coverage_dir/html"

find_llvm_tool() {
    local name=$1
    local candidate
    for candidate in "$name" "$name-18" "$name-19" "$name-20" "$name-21"; do
        if command -v "$candidate" >/dev/null 2>&1; then
            command -v "$candidate"
            return 0
        fi
    done
    echo "required LLVM tool not found: $name" >&2
    return 1
}

llvm_profdata=$(find_llvm_tool llvm-profdata)
llvm_cov=$(find_llvm_tool llvm-cov)

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

: > "$coverage_dir/report.txt"
for target in "${targets[@]}"; do
    blst_fuzz_target_limits "$target"
    profile="$coverage_dir/$target.profraw"
    LLVM_PROFILE_FILE="$profile" "$build_dir/fuzz/$target" \
        -runs=0 \
        -max_len="$BLST_FUZZ_MAX_LEN" \
        -timeout="$BLST_FUZZ_TIMEOUT" \
        -rss_limit_mb=2048 \
        "fuzz/corpus/$target"
    profdata="$coverage_dir/$target.profdata"
    "$llvm_profdata" merge -sparse "$profile" -o "$profdata"
    python3 scripts/fuzz/check_coverage_paths.py \
        --target "$target" \
        --llvm-cov "$llvm_cov" \
        --binary "$build_dir/fuzz/$target" \
        --profile "$profdata"
    echo "===== $target =====" >> "$coverage_dir/report.txt"
    "$llvm_cov" report \
        "$build_dir/fuzz/$target" \
        -instr-profile="$profdata" \
        -ignore-filename-regex='(fuzz/|test/|_deps/)' \
        >> "$coverage_dir/report.txt"
    "$llvm_cov" show \
        "$build_dir/fuzz/$target" \
        -instr-profile="$profdata" \
        -ignore-filename-regex='(fuzz/|test/|_deps/)' \
        -format=html \
        -output-dir="$coverage_dir/html/$target"
done

echo "coverage report: $coverage_dir/report.txt"

#!/usr/bin/env bash
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

if [[ ${BLST_FUZZ_DISCOVERY_APPROVED:-} != YES ]]; then
    echo "set BLST_FUZZ_DISCOVERY_APPROVED=YES for an approved long run" >&2
    exit 2
fi
if [[ $# -gt 2 ]]; then
    echo "usage: $0 [SECONDS] [TARGET]" >&2
    exit 2
fi

readonly seconds=${1:-900}
if ! [[ "$seconds" =~ ^[1-9][0-9]*$ ]]; then
    echo "SECONDS must be a positive integer" >&2
    exit 2
fi

readonly toolchain=nightly-2026-01-15
readonly cargo_fuzz_version=0.13.2
readonly fuzz_root=bindings/rust/fuzz
readonly lock_dir=out/fuzz-rust/.single-fuzz-worker.lock
readonly default_targets=(
    blst_rust_fuzz_decode
    blst_rust_fuzz_state
    blst_rust_fuzz_bls_api
)
if [[ $# -eq 2 ]]; then
    targets=("$2")
else
    targets=("${default_targets[@]}")
fi
case "${targets[0]}" in
    blst_rust_fuzz_decode|blst_rust_fuzz_state|blst_rust_fuzz_bls_api) ;;
    *) echo "unknown Rust fuzz target: ${targets[0]}" >&2; exit 2 ;;
esac

export CARGO_BUILD_JOBS=2
export CMAKE_BUILD_PARALLEL_LEVEL=2

version=$(cargo fuzz --version)
if [[ "$version" != "cargo-fuzz $cargo_fuzz_version" ]]; then
    echo "cargo-fuzz $cargo_fuzz_version is required; found: $version" >&2
    exit 2
fi

mkdir -p "$(dirname "$lock_dir")"
if ! mkdir "$lock_dir" 2>/dev/null; then
    echo "another fuzz runner owns $lock_dir" >&2
    exit 1
fi
trap 'rmdir "$lock_dir"' EXIT

export CFLAGS="-fsanitize=address -fno-omit-frame-pointer"
export CXXFLAGS="-fsanitize=address -fno-omit-frame-pointer"
export BLST_FUZZ_NATIVE_SANITIZER=address
export CARGO_TARGET_DIR="$fuzz_root/target/address"

for target in "${targets[@]}"; do
    if [[ "$target" == blst_rust_fuzz_decode ]]; then
        timeout_seconds=2
        max_len=192
    elif [[ "$target" == blst_rust_fuzz_bls_api ]]; then
        timeout_seconds=10
        max_len=512
    else
        timeout_seconds=10
        max_len=4096
    fi
    artifact_dir="out/fuzz-rust/artifacts/$target"
    working_corpus="out/fuzz-rust/discovery-corpus/$target"
    mkdir -p "$artifact_dir"
    mkdir -p "$working_corpus"
    cp -n "$fuzz_root/corpus/$target"/* "$working_corpus"/
    cargo "+$toolchain" fuzz run \
        --fuzz-dir "$fuzz_root" \
        --sanitizer address \
        "$target" "$working_corpus" -- \
        -max_total_time="$seconds" \
        -timeout="$timeout_seconds" \
        -rss_limit_mb=2048 \
        -max_len="$max_len" \
        -artifact_prefix="$artifact_dir/" 2>&1 | tee "$artifact_dir/discovery.log"
done

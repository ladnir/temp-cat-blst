#!/usr/bin/env bash
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

readonly toolchain=nightly-2026-01-15
readonly cargo_fuzz_version=0.13.2
readonly manifest=bindings/rust/fuzz/Cargo.toml
readonly fuzz_root=bindings/rust/fuzz
readonly lock_dir=out/fuzz-rust/.single-fuzz-worker.lock
readonly targets=(
    blst_rust_fuzz_decode
    blst_rust_fuzz_state
    blst_rust_fuzz_bls_api
)

profiles=("$@")
if [[ ${#profiles[@]} -eq 0 ]]; then
    profiles=(address undefined)
fi
for profile in "${profiles[@]}"; do
    if [[ "$profile" != address && "$profile" != undefined ]]; then
        echo "usage: $0 [address|undefined ...]" >&2
        exit 2
    fi
done

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

for sanitizer in "${profiles[@]}"; do
    case "$sanitizer" in
        address)
            native_sanitizer=address
            cargo_sanitizer=address
            ;;
        undefined)
            native_sanitizer=undefined
            cargo_sanitizer=none
            ;;
    esac
    export CFLAGS="-fsanitize=$native_sanitizer -fno-omit-frame-pointer"
    export CXXFLAGS="-fsanitize=$native_sanitizer -fno-omit-frame-pointer"
    export BLST_FUZZ_NATIVE_SANITIZER="$native_sanitizer"
    export CARGO_TARGET_DIR="$fuzz_root/target/$sanitizer"

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
        artifact_dir="out/fuzz-rust/artifacts/$target/$sanitizer"
        mkdir -p "$artifact_dir"
        cargo "+$toolchain" fuzz run \
            --fuzz-dir "$fuzz_root" \
            --sanitizer "$cargo_sanitizer" \
            "$target" "$fuzz_root/corpus/$target" -- \
            -runs=0 \
            -timeout="$timeout_seconds" \
            -rss_limit_mb=2048 \
            -max_len="$max_len" \
            -artifact_prefix="$artifact_dir/"
    done
done

#!/usr/bin/env bash
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

if [[ ${BLST_NATIVE_AVX_FUZZ_APPROVED:-} != YES ]]; then
    echo "set BLST_NATIVE_AVX_FUZZ_APPROVED=YES before native AVX discovery" >&2
    exit 2
fi
if [[ ${BLST_FUZZ_DISCOVERY_APPROVED:-} != YES ]]; then
    echo "set BLST_FUZZ_DISCOVERY_APPROVED=YES for an approved long run" >&2
    exit 2
fi
if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "usage: $0 BUILD_DIR [SECONDS]" >&2
    exit 2
fi

build_dir=$1
seconds=${2:-900}
if ! [[ "$seconds" =~ ^[1-9][0-9]*$ ]]; then
    echo "SECONDS must be a positive integer" >&2
    exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/target_limits.sh"
target=blst_fuzz_avx
blst_fuzz_target_limits "$target"

lock_dir="$build_dir/.single-fuzz-worker.lock"
if ! mkdir "$lock_dir" 2>/dev/null; then
    echo "another fuzz runner owns $lock_dir" >&2
    exit 1
fi
trap 'rmdir "$lock_dir"' EXIT

cmake --build "$build_dir" --target "$target" --parallel 2
artifact_dir="$build_dir/artifacts/$target"
working_corpus="$build_dir/discovery-corpus/$target"
mkdir -p "$artifact_dir"
mkdir -p "$working_corpus"
cp -n "fuzz/corpus/$target"/* "$working_corpus"/
"$build_dir/fuzz/$target" \
    -max_total_time="$seconds" \
    -max_len="$BLST_FUZZ_MAX_LEN" \
    -timeout="$BLST_FUZZ_TIMEOUT" \
    -rss_limit_mb=2048 \
    -artifact_prefix="$artifact_dir/" \
    "$working_corpus" 2>&1 | tee "$artifact_dir/discovery.log"

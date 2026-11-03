#!/usr/bin/env bash
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/target_limits.sh"

if [[ ${BLST_FUZZ_DISCOVERY_APPROVED:-} != YES ]]; then
    echo "set BLST_FUZZ_DISCOVERY_APPROVED=YES for an approved long run" >&2
    exit 2
fi
if [[ $# -lt 1 ]]; then
    echo "usage: $0 BUILD_DIR [SECONDS] [TARGET ...]" >&2
    exit 2
fi

build_dir=$1
seconds=${2:-900}
if ! [[ "$seconds" =~ ^[1-9][0-9]*$ ]]; then
    echo "SECONDS must be a positive integer" >&2
    exit 2
fi
shift $(( $# >= 2 ? 2 : 1 ))

default_targets=(
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
targets=("${@:-}")
if [[ ${#targets[@]} -eq 0 || -z ${targets[0]} ]]; then
    targets=("${default_targets[@]}")
fi

is_known_target() {
    local candidate=$1
    local known
    for known in "${default_targets[@]}"; do
        [[ "$candidate" == "$known" ]] && return 0
    done
    return 1
}

lock_dir="$build_dir/.single-fuzz-worker.lock"
if ! mkdir "$lock_dir" 2>/dev/null; then
    echo "another fuzz runner owns $lock_dir" >&2
    exit 1
fi
trap 'rmdir "$lock_dir"' EXIT

for target in "${targets[@]}"; do
    if ! is_known_target "$target"; then
        echo "unknown scalar fuzz target: $target" >&2
        exit 2
    fi
    cmake --build "$build_dir" --target "$target" --parallel 2
    blst_fuzz_target_limits "$target"
    artifact_dir="$build_dir/artifacts/$target"
    working_corpus="$build_dir/discovery-corpus/$target"
    mkdir -p "$artifact_dir"
    mkdir -p "$working_corpus"
    cp -n "test/fuzz/corpus/$target"/* "$working_corpus"/
    "$build_dir/test/fuzz/$target" \
        -jobs=1 \
        -workers=1 \
        -max_total_time="$seconds" \
        -max_len="$BLST_FUZZ_MAX_LEN" \
        -timeout="$BLST_FUZZ_TIMEOUT" \
        -rss_limit_mb=2048 \
        -artifact_prefix="$artifact_dir/" \
        "$working_corpus" 2>&1 | tee "$artifact_dir/discovery.log"
done

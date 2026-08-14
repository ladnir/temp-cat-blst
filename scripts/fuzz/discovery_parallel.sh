#!/usr/bin/env bash
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/target_limits.sh"

if [[ ${BLST_FUZZ_DISCOVERY_APPROVED:-} != YES ||
      ${BLST_FUZZ_PARALLEL_APPROVED:-} != YES ]]; then
    echo "set both BLST_FUZZ_DISCOVERY_APPROVED=YES and BLST_FUZZ_PARALLEL_APPROVED=YES" >&2
    exit 2
fi
if [[ $# -ne 3 ]]; then
    echo "usage: $0 BUILD_DIR SECONDS TOTAL_WORKERS" >&2
    exit 2
fi

repo_root=$(pwd -P)
build_dir=$(realpath "$1")
seconds=$2
total_workers=$3
if ! [[ "$seconds" =~ ^[1-9][0-9]*$ && "$total_workers" =~ ^[1-9][0-9]*$ ]]; then
    echo "SECONDS and TOTAL_WORKERS must be positive integers" >&2
    exit 2
fi

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

cpu_count=$(nproc)
if (( total_workers < ${#targets[@]} )); then
    echo "TOTAL_WORKERS must be at least ${#targets[@]} to run every target concurrently" >&2
    exit 2
fi
if (( total_workers > cpu_count )); then
    echo "TOTAL_WORKERS=$total_workers exceeds the detected $cpu_count logical CPUs" >&2
    exit 2
fi

# Give every target one worker. Additional workers go first to the campaigns
# with the most expensive arithmetic and state oracles.
worker_counts=(1 1 1 1 1 1 1 1 1)
priority=(2 3 4 5 7 8 0 1 6)
remaining=$((total_workers - ${#targets[@]}))
for ((i = 0; remaining > 0; ++i, --remaining)); do
    index=${priority[$((i % ${#priority[@]}))]}
    worker_counts[$index]=$((worker_counts[$index] + 1))
done

# This compatibility lock is shared with the sequential scalar and AVX
# runners. Its historical name describes the old policy, not this campaign.
lock_dir="$build_dir/.single-fuzz-worker.lock"
if ! mkdir "$lock_dir" 2>/dev/null; then
    echo "another fuzz runner owns $lock_dir" >&2
    exit 1
fi

pids=()
cleanup() {
    local status=$?
    trap - EXIT INT TERM
    if (( ${#pids[@]} != 0 )); then
        kill "${pids[@]}" 2>/dev/null || true
        wait "${pids[@]}" 2>/dev/null || true
    fi
    rmdir "$lock_dir"
    exit "$status"
}
trap cleanup EXIT INT TERM

cmake --build "$build_dir" --target "${targets[@]}" --parallel 2

campaign=$(date -u +%Y%m%dT%H%M%SZ)-$$
hard_seconds=$((seconds + 120))
run_root="$build_dir/parallel-work/$campaign"
mkdir -p "$run_root"

for i in "${!targets[@]}"; do
    target=${targets[$i]}
    worker_count=${worker_counts[$i]}
    blst_fuzz_target_limits "$target"
    artifact_dir="$build_dir/artifacts/$target"
    working_corpus="$build_dir/discovery-corpus/$target"
    target_run_dir="$run_root/$target"
    log="$artifact_dir/parallel-$campaign.log"
    mkdir -p "$artifact_dir" "$working_corpus" "$target_run_dir"
    cp --update=none "$repo_root/fuzz/corpus/$target"/* "$working_corpus"/

    echo "START $target workers=$worker_count log=$log"
    (
        cd "$target_run_dir"
        if (( worker_count == 1 )); then
            exec timeout --signal=TERM --kill-after=15s "$hard_seconds" \
                "$build_dir/fuzz/$target" \
                -max_total_time="$seconds" \
                -max_len="$BLST_FUZZ_MAX_LEN" \
                -timeout="$BLST_FUZZ_TIMEOUT" \
                -rss_limit_mb=2048 \
                -artifact_prefix="$artifact_dir/" \
                "$working_corpus"
        fi
        exec timeout --signal=TERM --kill-after=15s "$hard_seconds" \
            "$build_dir/fuzz/$target" \
            -jobs="$worker_count" \
            -workers="$worker_count" \
            -max_total_time="$seconds" \
            -max_len="$BLST_FUZZ_MAX_LEN" \
            -timeout="$BLST_FUZZ_TIMEOUT" \
            -rss_limit_mb=2048 \
            -artifact_prefix="$artifact_dir/" \
            "$working_corpus"
    ) >"$log" 2>&1 &
    pids+=("$!")
done

status=0
for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then
        echo "DONE ${targets[$i]} workers=${worker_counts[$i]}"
    else
        echo "FAILED ${targets[$i]} workers=${worker_counts[$i]}" >&2
        status=1
    fi
done

pids=()
exit "$status"

#!/usr/bin/env bash

set -euo pipefail

readonly expected_clang_version=18.1.3
readonly checks='-*,bugprone-*,performance-*,portability-*,-bugprone-easily-swappable-parameters,-bugprone-reserved-identifier'

usage()
{
    cat <<'EOF'
Usage: run_clang_tidy_local.sh [scalar|avx|all]

Optional environment variables:
  BLST_CLANG_TIDY_JOBS          clang-tidy workers (default: 2)
  BLST_CLANG_TIDY_WORK_PARENT   WSL-native work parent
                                (default: /var/tmp/catblst-clang-tidy)
EOF
}

mode=${1:-all}
case "$mode" in
    scalar|avx|all)
        ;;
    -h|--help)
        usage
        exit 0
        ;;
    *)
        usage >&2
        exit 2
        ;;
esac

if [[ $(uname -s) != Linux || $(uname -m) != x86_64 ]]; then
    echo 'This runner requires an x86_64 Linux environment such as WSL2.' >&2
    exit 1
fi

for tool in clang clang++ clang-tidy run-clang-tidy cmake ninja; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Required tool not found: $tool" >&2
        exit 1
    fi
done

clang_version=$(clang-tidy --version)
if ! grep -Fq "LLVM version $expected_clang_version" <<<"$clang_version"; then
    printf 'Expected clang-tidy %s, but received:\n%s\n' \
        "$expected_clang_version" "$clang_version" >&2
    exit 1
fi

repo_root=$(git rev-parse --show-toplevel)
commit=$(git -C "$repo_root" rev-parse HEAD)
short_commit=${commit:0:12}
timestamp=$(date -u +%Y%m%dT%H%M%SZ)
run_id="$short_commit-$timestamp"
jobs=${BLST_CLANG_TIDY_JOBS:-2}
work_parent=${BLST_CLANG_TIDY_WORK_PARENT:-/var/tmp/catblst-clang-tidy}
work_root="$work_parent/$run_id"
results_root="$repo_root/out/audit/clang-tidy/$run_id"
source_root="$work_root/source"
source_list="$work_root/source-files.nul"
source_archive="$work_root/source-snapshot.tar"

if [[ -e $work_root || -e $results_root ]]; then
    echo "Refusing to overwrite an existing clang-tidy run: $run_id" >&2
    exit 1
fi

mkdir -p "$work_root" "$results_root" "$source_root"
exec > >(tee "$results_root/run.log") 2>&1

(
    cd "$repo_root"
    git ls-files --cached --others --exclude-standard -z |
        while IFS= read -r -d '' path; do
            if [[ -f $path || -L $path ]]; then
                printf '%s\0' "$path"
            fi
        done >"$source_list"
    tar --null --no-recursion --files-from="$source_list" \
        -cf "$source_archive"
)

source_sha256=$(sha256sum "$source_archive" | cut -d ' ' -f 1)
tar -xf "$source_archive" -C "$source_root"
rm -f -- "$source_archive"
tr '\0' '\n' <"$source_list" >"$results_root/source-files.txt"
git -C "$repo_root" status --porcelain=v1 --untracked-files=all \
    >"$results_root/working-tree-status.txt"

cat >"$results_root/run-metadata.txt" <<EOF
commit=$commit
timestamp_utc=$timestamp
mode=$mode
clang_tidy_version=$expected_clang_version
checks=$checks
jobs=$jobs
work_root=$work_root
source_root=$source_root
source_snapshot_sha256=$source_sha256
EOF

printf '%s\n' "$clang_version" >"$results_root/clang-tidy-version.txt"
cmake --version >"$results_root/cmake-version.txt"

run_profile()
{
    local profile=$1
    local build_root="$work_root/$profile"
    local profile_results="$results_root/$profile"
    local profile_flags=()

    case "$profile" in
        scalar)
            profile_flags=(
                -DBLST_AVX512=OFF
                -DBLST_AVX512_PORTABLE=OFF
                -DBLST_DEFAULT_AVX=OFF
            )
            ;;
        avx-portable)
            profile_flags=(
                -DBLST_AVX512=OFF
                -DBLST_AVX512_PORTABLE=ON
                -DBLST_DEFAULT_AVX=ON
            )
            ;;
    esac

    mkdir -p "$profile_results/fixes"
    cmake -S "$source_root" -B "$build_root" \
        -G Ninja \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_C_COMPILER=clang \
        -DCMAKE_CXX_COMPILER=clang++ \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        -DBLST_BUILD_SHARED=OFF \
        -DBLST_PORTABLE=ON \
        -DBLST_ASM=OFF \
        -DBLST_ADX=OFF \
        -DBLST_TEST=OFF \
        -DBLST_BENCHMARK=OFF \
        -DBLST_CUDA=OFF \
        -DBLST_CUDA_MSM=OFF \
        -DBLST_PARALLEL_MSM=OFF \
        -DBLST_CRNS_BATCH_VROOM_VENDOR=OFF \
        "${profile_flags[@]}" \
        2>&1 | tee "$profile_results/configure.log"

    run-clang-tidy \
        -p "$build_root" \
        -j "$jobs" \
        -quiet \
        -checks="$checks" \
        -header-filter="^$source_root/" \
        -export-fixes="$profile_results/fixes" \
        "^$source_root/blst/.*" \
        2>&1 | tee "$profile_results/clang-tidy.log"
}

case "$mode" in
    scalar)
        run_profile scalar
        ;;
    avx)
        run_profile avx-portable
        ;;
    all)
        run_profile scalar
        run_profile avx-portable
        ;;
esac

printf '\nclang-tidy results: %s\n' "$results_root"

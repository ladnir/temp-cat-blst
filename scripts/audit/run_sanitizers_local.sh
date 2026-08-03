#!/usr/bin/env bash

set -euo pipefail

readonly expected_clang_version=18.1.3
readonly test_filter='BackendAlias.*:BindingHeader.*:Blst256Shift.*:BlstTemplateDispatch.*:CApi.*:CompressedInfinity.*:Conformance.*:FftApi.*:FrMontgomery.*:MiddleProductCApi.*:InterfaceTests.*:ScalarMultWindow.*-BindingHeader.NativeFrBatchAffineG1DifferentialDistributions:BindingHeader.NativeFrBatchAffineG2MatchesPippenger:BindingHeader.SelectedNativeFrG1MsmHandlesCapacityAndIdentity:BindingHeader.SelectedNativeFrG2MsmHandlesBatchRange:CApi.ScalarG1WbitsPrecomputeCrossesScratchThresholdWithoutAllocation:CApi.ScalarG2WbitsPrecomputeCrossesScratchThresholdWithoutAllocation:InterfaceTests.Div3LimbsRandomFuzz'

usage()
{
    cat <<'EOF'
Usage: run_sanitizers_local.sh [asan|ubsan|all]

Optional environment variables:
  BLST_SANITIZER_JOBS          build workers (default: 2)
  BLST_SANITIZER_WORK_PARENT   Linux-native work parent
                               (default: /var/tmp/catblst-sanitizers)
EOF
}

mode=${1:-all}
case "$mode" in
    asan|ubsan|all)
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

for tool in clang clang++ cmake ninja git tar sha256sum; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Required tool not found: $tool" >&2
        exit 1
    fi
done

clang_version=$(clang --version)
if ! grep -Fq "clang version $expected_clang_version" <<<"$clang_version"; then
    printf 'Expected Clang %s, but received:\n%s\n' \
        "$expected_clang_version" "$clang_version" >&2
    exit 1
fi

repo_root=$(git rev-parse --show-toplevel)
commit=$(git -C "$repo_root" rev-parse HEAD)
short_commit=${commit:0:12}
timestamp=$(date -u +%Y%m%dT%H%M%SZ)
run_id="$short_commit-$timestamp"
jobs=${BLST_SANITIZER_JOBS:-2}
work_parent=${BLST_SANITIZER_WORK_PARENT:-/var/tmp/catblst-sanitizers}
work_root="$work_parent/$run_id"
results_root="$repo_root/out/audit/sanitizers/$run_id"
source_root="$work_root/source"
source_list="$work_root/source-files.nul"
source_archive="$work_root/source-snapshot.tar"

if [[ -e $work_root || -e $results_root ]]; then
    echo "Refusing to overwrite an existing sanitizer run: $run_id" >&2
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
clang_version=$expected_clang_version
build_type=Debug
jobs=$jobs
test_filter=$test_filter
work_root=$work_root
source_root=$source_root
source_snapshot_sha256=$source_sha256
EOF

printf '%s\n' "$clang_version" >"$results_root/clang-version.txt"
cmake --version >"$results_root/cmake-version.txt"

symbolizer=$(command -v llvm-symbolizer-18 || command -v llvm-symbolizer || true)

run_profile()
{
    local profile=$1
    local build_root="$work_root/$profile"
    local profile_results="$results_root/$profile"
    local asan=OFF
    local ubsan=OFF

    case "$profile" in
        asan)
            asan=ON
            ;;
        ubsan)
            ubsan=ON
            ;;
    esac

    mkdir -p "$profile_results"
    cmake -S "$source_root" -B "$build_root" \
        -G Ninja \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_C_COMPILER=clang \
        -DCMAKE_CXX_COMPILER=clang++ \
        -DBLST_BUILD_SHARED=OFF \
        -DBLST_PORTABLE=ON \
        -DBLST_ASM=OFF \
        -DBLST_ADX=OFF \
        -DBLST_DEFAULT_ASM=OFF \
        -DBLST_AVX512=OFF \
        -DBLST_AVX512_PORTABLE=OFF \
        -DBLST_DEFAULT_AVX=OFF \
        -DBLST_TEST=ON \
        -DBLST_BENCHMARK=OFF \
        -DBLST_CUDA=OFF \
        -DBLST_CUDA_MSM=OFF \
        -DBLST_PARALLEL_MSM=OFF \
        -DBLST_CRNS_BATCH_VROOM_VENDOR=OFF \
        -DBLST_ASAN="$asan" \
        -DBLST_UBSAN="$ubsan" \
        2>&1 | tee "$profile_results/configure.log"

    cmake --build "$build_root" --target blst_sanitizer_boundary_tests \
        --parallel "$jobs" \
        2>&1 | tee "$profile_results/build.log"

    local test_executable="$build_root/test/blst_sanitizer_boundary_tests"
    if [[ ! -x $test_executable ]]; then
        echo "Sanitizer test executable not found: $test_executable" >&2
        exit 1
    fi

    case "$profile" in
        asan)
            ASAN_OPTIONS='detect_leaks=1:halt_on_error=1:strict_string_checks=1:check_initialization_order=1' \
            LSAN_OPTIONS='exitcode=23' \
            ASAN_SYMBOLIZER_PATH="$symbolizer" \
                "$test_executable" --gtest_filter="$test_filter" \
                2>&1 | tee "$profile_results/tests.log"
            ;;
        ubsan)
            UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1' \
                "$test_executable" --gtest_filter="$test_filter" \
                2>&1 | tee "$profile_results/tests.log"
            ;;
    esac

    if ! grep -Eq '^\[  PASSED  \] [1-9][0-9]* tests?\.$' \
        "$profile_results/tests.log"; then
        echo "Profile '$profile' did not report a nonzero passing test count." >&2
        exit 1
    fi
}

case "$mode" in
    asan)
        run_profile asan
        ;;
    ubsan)
        run_profile ubsan
        ;;
    all)
        run_profile asan
        run_profile ubsan
        ;;
esac

printf '\nSanitizer results: %s\n' "$results_root"

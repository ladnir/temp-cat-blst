#!/usr/bin/env bash

set -euo pipefail

repo_root=${BLST_CODEQL_SOURCE_ROOT:-}
if [[ -z $repo_root ]]; then
    repo_root=$(git rev-parse --show-toplevel)
fi

if [[ ! -f $repo_root/CMakeLists.txt ]]; then
    echo "Invalid CodeQL source root: $repo_root" >&2
    exit 1
fi
build_root=$(mktemp -d /var/tmp/catblst-codeql-build.XXXXXX)
build_jobs=${BLST_CODEQL_BUILD_JOBS:-2}
profile=${BLST_CODEQL_PROFILE:-}

case "$profile" in
    scalar|avx-portable)
        ;;
    *)
        echo 'BLST_CODEQL_PROFILE must be scalar or avx-portable.' >&2
        exit 2
        ;;
esac

cleanup()
{
    case "$build_root" in
        /var/tmp/catblst-codeql-build.*)
            rm -rf -- "$build_root"
            ;;
        *)
            printf 'Refusing to remove unexpected build directory: %s\n' \
                "$build_root" >&2
            return 1
            ;;
    esac
}
trap cleanup EXIT

common_flags=(
    -G Ninja
    -DCMAKE_BUILD_TYPE=Debug
    -DBLST_BUILD_SHARED=OFF
    -DBLST_PORTABLE=ON
    -DBLST_ASM=OFF
    -DBLST_ADX=OFF
    -DBLST_TEST=OFF
    -DBLST_BENCHMARK=OFF
    -DBLST_CUDA=OFF
    -DBLST_CUDA_MSM=OFF
    -DBLST_PARALLEL_MSM=OFF
    -DBLST_CRNS_BATCH_VROOM_VENDOR=OFF
)

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

cmake -S "$repo_root" -B "$build_root/$profile" \
    "${common_flags[@]}" \
    "${profile_flags[@]}"
cmake --build "$build_root/$profile" \
    --target blst_static \
    --parallel "$build_jobs"

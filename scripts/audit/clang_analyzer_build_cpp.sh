#!/usr/bin/env bash

set -euo pipefail

repo_root=${BLST_CLANG_ANALYZER_SOURCE_ROOT:-}
profile=${BLST_CLANG_ANALYZER_PROFILE:-}

if [[ ! -f $repo_root/CMakeLists.txt ]]; then
    echo "Invalid Clang analyzer source root: $repo_root" >&2
    exit 1
fi

case "$profile" in
    scalar|avx-portable)
        ;;
    *)
        echo 'BLST_CLANG_ANALYZER_PROFILE must be scalar or avx-portable.' >&2
        exit 2
        ;;
esac

build_root=$(mktemp -d /var/tmp/catblst-clang-analyzer-build.XXXXXX)
build_jobs=${BLST_CLANG_ANALYZER_BUILD_JOBS:-2}

cleanup()
{
    case "$build_root" in
        /var/tmp/catblst-clang-analyzer-build.*)
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
    -DCMAKE_C_COMPILER="${CC:-clang}"
    -DCMAKE_CXX_COMPILER="${CXX:-clang++}"
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
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

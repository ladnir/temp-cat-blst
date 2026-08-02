#!/usr/bin/env bash

set -euo pipefail

readonly expected_clang_version=18.1.3

usage()
{
    cat <<'EOF'
Usage: run_clang_analyzer_local.sh [scalar|avx|all]

Optional environment variables:
  BLST_CLANG_ANALYZER_BUILD_JOBS   CMake build workers (default: 2)
  BLST_CLANG_ANALYZER_WORK_PARENT  WSL-native work parent
                                   (default: /var/tmp/catblst-clang-analyzer)
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

for tool in clang clang++ scan-build cmake ninja; do
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
work_parent=${BLST_CLANG_ANALYZER_WORK_PARENT:-/var/tmp/catblst-clang-analyzer}
work_root="$work_parent/$run_id"
results_root="$repo_root/out/audit/clang-analyzer/$run_id"
source_root="$work_root/source"
source_list="$work_root/source-files.nul"
source_archive="$work_root/source-snapshot.tar"

if [[ -e $work_root || -e $results_root ]]; then
    echo "Refusing to overwrite an existing analyzer run: $run_id" >&2
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
build_jobs=${BLST_CLANG_ANALYZER_BUILD_JOBS:-2}
work_root=$work_root
source_root=$source_root
source_snapshot_sha256=$source_sha256
EOF

printf '%s\n' "$clang_version" >"$results_root/clang-version.txt"
scan-build --help >"$results_root/scan-build-help.txt"
cmake --version >"$results_root/cmake-version.txt"

cd "$source_root"

run_profile()
{
    local profile=$1
    local result_name=$2

    scan-build \
        -analyze-headers \
        -sarif \
        --use-cc clang \
        --use-c++ clang++ \
        -o "$results_root/$result_name" \
        env \
        BLST_CLANG_ANALYZER_SOURCE_ROOT="$source_root" \
        BLST_CLANG_ANALYZER_PROFILE="$profile" \
        bash scripts/audit/clang_analyzer_build_cpp.sh
}

case "$mode" in
    scalar)
        run_profile scalar scalar
        ;;
    avx)
        run_profile avx-portable avx-portable
        ;;
    all)
        run_profile scalar scalar
        run_profile avx-portable avx-portable
        ;;
esac

printf '\nClang analyzer results: %s\n' "$results_root"

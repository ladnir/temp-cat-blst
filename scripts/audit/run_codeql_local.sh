#!/usr/bin/env bash

set -euo pipefail

readonly expected_codeql_version=2.26.2
readonly cpp_query_suite='codeql/cpp-queries@1.8.0:codeql-suites/cpp-security-and-quality.qls'
readonly rust_query_suite='codeql/rust-queries@0.1.39:codeql-suites/rust-security-and-quality.qls'

usage()
{
    cat <<'EOF'
Usage: CODEQL_EXE=/path/to/codeql run_codeql_local.sh [cpp|rust|all]

Optional environment variables:
  BLST_CODEQL_THREADS       CodeQL worker count (default: 2)
  BLST_CODEQL_BUILD_JOBS    CMake build worker count (default: 2)
  BLST_CODEQL_WORK_PARENT   WSL-native database parent
                            (default: /var/tmp/catblst-codeql)
EOF
}

mode=${1:-all}
case "$mode" in
    cpp|rust|all)
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

if [[ -z ${CODEQL_EXE:-} || ! -x $CODEQL_EXE ]]; then
    echo 'CODEQL_EXE must name the executable from the verified CodeQL bundle.' >&2
    exit 1
fi

codeql_version=$($CODEQL_EXE version 2>&1)
if ! grep -Fq "CodeQL command-line toolchain release $expected_codeql_version" \
    <<<"$codeql_version"; then
    printf 'Expected CodeQL %s, but received:\n%s\n' \
        "$expected_codeql_version" "$codeql_version" >&2
    exit 1
fi

repo_root=$(git rev-parse --show-toplevel)
commit=$(git -C "$repo_root" rev-parse HEAD)
short_commit=${commit:0:12}
timestamp=$(date -u +%Y%m%dT%H%M%SZ)
run_id="$short_commit-$timestamp"
threads=${BLST_CODEQL_THREADS:-2}
work_parent=${BLST_CODEQL_WORK_PARENT:-/var/tmp/catblst-codeql}
work_root="$work_parent/$run_id"
results_root="$repo_root/out/audit/codeql/$run_id"
source_root="$work_root/source"
source_list="$work_root/source-files.nul"
source_archive="$work_root/source-snapshot.tar"

if [[ -e $work_root || -e $results_root ]]; then
    echo "Refusing to overwrite an existing CodeQL run: $run_id" >&2
    exit 1
fi

mkdir -p "$work_root" "$results_root" "$source_root"
exec > >(tee "$results_root/run.log") 2>&1

# Build-free extraction otherwise discovers ignored compiler output under
# out/ and bindings/rust/target. Snapshot the current non-ignored checkout so
# local edits remain visible without allowing stale generated files into the
# database.
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
codeql_version=$expected_codeql_version
cpp_query_suite=$cpp_query_suite
rust_query_suite=$rust_query_suite
threads=$threads
build_jobs=${BLST_CODEQL_BUILD_JOBS:-2}
work_root=$work_root
source_root=$source_root
source_snapshot_sha256=$source_sha256
EOF

printf '%s\n' "$codeql_version" >"$results_root/codeql-version.txt"
cmake --version >"$results_root/cmake-version.txt"
c++ --version >"$results_root/cxx-version.txt"
if command -v rustc >/dev/null 2>&1; then
    rustc --version --verbose >"$results_root/rust-version.txt"
fi

cd "$source_root"

run_cpp()
{
    "$CODEQL_EXE" database create "$work_root/cpp-db" \
        --language=c-cpp \
        --source-root="$source_root" \
        --threads="$threads" \
        --command="BLST_CODEQL_SOURCE_ROOT='$source_root' bash scripts/audit/codeql_build_cpp.sh"

    "$CODEQL_EXE" database analyze "$work_root/cpp-db" \
        "$cpp_query_suite" \
        --format=sarifv2.1.0 \
        --sarif-category=c-cpp \
        --sarif-add-baseline-file-info \
        --threads="$threads" \
        --output="$results_root/cpp.sarif"
}

run_rust()
{
    "$CODEQL_EXE" database create "$work_root/rust-db" \
        --language=rust \
        --build-mode=none \
        --source-root="$source_root" \
        --threads="$threads"

    "$CODEQL_EXE" database analyze "$work_root/rust-db" \
        "$rust_query_suite" \
        --format=sarifv2.1.0 \
        --sarif-category=rust \
        --sarif-add-baseline-file-info \
        --threads="$threads" \
        --output="$results_root/rust.sarif"
}

case "$mode" in
    cpp)
        run_cpp
        ;;
    rust)
        run_rust
        ;;
    all)
        run_cpp
        run_rust
        ;;
esac

printf '\nCodeQL results: %s\n' "$results_root"
printf 'CodeQL databases: %s\n' "$work_root"
printf 'Retain the databases until extraction diagnostics are triaged.\n'

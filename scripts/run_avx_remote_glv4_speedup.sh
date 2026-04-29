#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PLINK_EXE="${PLINK_EXE:-plink}"
PLINK_KEY="${PLINK_KEY:-/path/to/key}"
REMOTE_HOST="${REMOTE_HOST:-18.188.231.170}"
REMOTE_USER="${REMOTE_USER:-ec2-user}"
REMOTE_HOSTKEY="${REMOTE_HOSTKEY:-ssh-ed25519 255 SHA256:/qd0qT37U/NAXgZT1kN6kSHaiyhoqQUXrlQBQXTa4Go}"
REMOTE_BASE="${REMOTE_BASE:-/tmp/avx-glv4}"
CLANG_ROOT="${CLANG_ROOT:-/opt/llvm-21/bin}"
PIN_CORE="${PIN_CORE:-0}"
REPETITIONS="${REPETITIONS:-5}"
MIN_TIME="${MIN_TIME:-0.05s}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
REMOTE_RUN_DIR="${REMOTE_BASE}/${RUN_ID}"
LOCAL_REPORT_DIR="${ROOT_DIR}/out/avx-glv4/${RUN_ID}"

resolve_windows_path()
{
    local candidate="$1"
    if [[ -e "${candidate}" ]]; then
        printf '%s\n' "${candidate}"
        return 0
    fi

    case "${candidate}" in
        /c/*)
            local alt="/mnt/${candidate#/}"
            if [[ -e "${alt}" ]]; then
                printf '%s\n' "${alt}"
                return 0
            fi
            ;;
        /mnt/c/*)
            local alt="/${candidate#/mnt/}"
            if [[ -e "${alt}" ]]; then
                printf '%s\n' "${alt}"
                return 0
            fi
            ;;
    esac

    return 1
}

to_windows_path()
{
    local candidate="$1"
    case "${candidate}" in
        /mnt/?/*)
            local drive_letter="${candidate:5:1}"
            local remainder="${candidate:7}"
            remainder="${remainder//\//\\}"
            printf '%s:\\%s\n' "${drive_letter^^}" "${remainder}"
            return 0
            ;;
        /?/*)
            local drive_letter="${candidate:1:1}"
            local remainder="${candidate:3}"
            remainder="${remainder//\//\\}"
            printf '%s:\\%s\n' "${drive_letter^^}" "${remainder}"
            return 0
            ;;
    esac

    printf '%s\n' "${candidate}"
}

PLINK_EXE="$(resolve_windows_path "${PLINK_EXE}")" || {
    echo "missing plink executable: ${PLINK_EXE}" >&2
    exit 1
}
PLINK_KEY="$(resolve_windows_path "${PLINK_KEY}")" || {
    echo "missing plink key: ${PLINK_KEY}" >&2
    exit 1
}
PLINK_KEY_WIN="$(to_windows_path "${PLINK_KEY}")"

mkdir -p "${LOCAL_REPORT_DIR}"

remote()
{
    "${PLINK_EXE}" -batch -hostkey "${REMOTE_HOSTKEY}" -i "${PLINK_KEY_WIN}" "${REMOTE_USER}@${REMOTE_HOST}" "$@"
}

upload()
{
    local src="$1"
    local dst="$2"
    cat "${src}" | "${PLINK_EXE}" -batch -hostkey "${REMOTE_HOSTKEY}" -i "${PLINK_KEY_WIN}" "${REMOTE_USER}@${REMOTE_HOST}" "cat > '${dst}'"
}

run_remote_json()
{
    local name="$1"
    local command="$2"
    local output_path="${LOCAL_REPORT_DIR}/${name}.json"
    remote "${command}" > "${output_path}"
    echo "${output_path}"
}

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

SNAPSHOT_ARCHIVE="${TMP_DIR}/blst_current.tgz"
SNAPSHOT_DIR="${TMP_DIR}/blst_snapshot"
mkdir -p "${SNAPSHOT_DIR}"

git -C "${ROOT_DIR}" archive HEAD | tar -xf - -C "${SNAPSHOT_DIR}"

while IFS= read -r -d '' path; do
    [[ "${path}" == "nul" ]] && continue
    if [[ -f "${ROOT_DIR}/${path}" ]]; then
        mkdir -p "${SNAPSHOT_DIR}/$(dirname "${path}")"
        cp "${ROOT_DIR}/${path}" "${SNAPSHOT_DIR}/${path}"
    else
        rm -f "${SNAPSHOT_DIR:?}/${path}"
    fi
done < <(git -C "${ROOT_DIR}" diff --name-only -z HEAD)

while IFS= read -r -d '' path; do
    [[ "${path}" == "nul" ]] && continue
    mkdir -p "${SNAPSHOT_DIR}/$(dirname "${path}")"
    cp "${ROOT_DIR}/${path}" "${SNAPSHOT_DIR}/${path}"
done < <(git -C "${ROOT_DIR}" ls-files --others --exclude-standard -z)

tar -czf "${SNAPSHOT_ARCHIVE}" -C "${SNAPSHOT_DIR}" .

echo "remote host: ${REMOTE_USER}@${REMOTE_HOST}"
echo "remote run dir: ${REMOTE_RUN_DIR}"
echo "local report dir: ${LOCAL_REPORT_DIR}"

remote "mkdir -p '${REMOTE_RUN_DIR}/archives'"
upload "${SNAPSHOT_ARCHIVE}" "${REMOTE_RUN_DIR}/archives/blst_current.tgz"
remote "rm -rf '${REMOTE_RUN_DIR}/blst_current' && mkdir -p '${REMOTE_RUN_DIR}/blst_current' && tar -xzf '${REMOTE_RUN_DIR}/archives/blst_current.tgz' -C '${REMOTE_RUN_DIR}/blst_current'"

remote "cd '${REMOTE_RUN_DIR}/blst_current' && cmake -S . -B out/glv4-clang21 -G Ninja -DCMAKE_C_COMPILER='${CLANG_ROOT}/clang' -DCMAKE_CXX_COMPILER='${CLANG_ROOT}/clang++' -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS_RELEASE='-O3 -march=native' -DCMAKE_CXX_FLAGS_RELEASE='-O3 -march=native' -DBLST_PORTABLE=OFF -DBLST_ASM=ON -DBLST_ADX=ON -DBLST_AVX512=ON -DBLST_DEFAULT_AVX=ON -DBLST_AVX512_PORTABLE=OFF -DBLST_BENCHMARK=ON -DBLST_TEST=ON -DBLST_CUDA=OFF >/dev/null && cmake --build out/glv4-clang21 --target blst_benchmarks blst_tests -j2 >/dev/null"

remote "cd '${REMOTE_RUN_DIR}/blst_current' && ./out/glv4-clang21/test/blst_tests --gtest_filter='AvxGtGlv4.*' >/dev/null"

SCALAR_BASIC_JSON="$(run_remote_json "scalar_basic" "cd '${REMOTE_RUN_DIR}/blst_current/out/glv4-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_GTScalarCyclotomicExp_Basic$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
SCALAR_GLV4_JSON="$(run_remote_json "scalar_glv4" "cd '${REMOTE_RUN_DIR}/blst_current/out/glv4-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_GTScalarCyclotomicExp_GLV4$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
SCALAR_GLV4_CONST_JSON="$(run_remote_json "scalar_glv4_constant" "cd '${REMOTE_RUN_DIR}/blst_current/out/glv4-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_GTScalarCyclotomicExp_GLV4_Constant$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
AVX_BASIC_JSON="$(run_remote_json "avx_basic_52" "cd '${REMOTE_RUN_DIR}/blst_current/out/glv4-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_AvxGTCyclotomicExp_Basic_52$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
AVX_GLV4_JSON="$(run_remote_json "avx_glv4_52" "cd '${REMOTE_RUN_DIR}/blst_current/out/glv4-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_AvxGTCyclotomicExp_GLV4_52$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
AVX_GLV4_CONST_JSON="$(run_remote_json "avx_glv4_constant_52" "cd '${REMOTE_RUN_DIR}/blst_current/out/glv4-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_AvxGTCyclotomicExp_GLV4_Constant_52$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"

SUMMARY_PATH="${LOCAL_REPORT_DIR}/summary.txt"
python3 - "${SCALAR_BASIC_JSON}" "${SCALAR_GLV4_JSON}" "${SCALAR_GLV4_CONST_JSON}" "${AVX_BASIC_JSON}" "${AVX_GLV4_JSON}" "${AVX_GLV4_CONST_JSON}" "${SUMMARY_PATH}" <<'PY'
import json
import pathlib
import sys

def load_median(path_str):
    data = json.loads(pathlib.Path(path_str).read_text())
    for bench in data["benchmarks"]:
        if bench.get("aggregate_name") == "median":
            return bench["name"], float(bench["cpu_time"]), bench["time_unit"]
    raise RuntimeError(f"median aggregate not found in {path_str}")

def speedup(baseline, optimized):
    return baseline / optimized

scalar_basic_name, scalar_basic, unit = load_median(sys.argv[1])
scalar_glv4_name, scalar_glv4, _ = load_median(sys.argv[2])
scalar_const_name, scalar_const, _ = load_median(sys.argv[3])
avx_basic_name, avx_basic, _ = load_median(sys.argv[4])
avx_glv4_name, avx_glv4, _ = load_median(sys.argv[5])
avx_const_name, avx_const, _ = load_median(sys.argv[6])

lines = [
    f"{scalar_basic_name}: {scalar_basic:.3f} {unit}",
    f"{scalar_glv4_name}: {scalar_glv4:.3f} {unit}",
    f"{scalar_const_name}: {scalar_const:.3f} {unit}",
    f"scalar glv4 speedup: {speedup(scalar_basic, scalar_glv4):.3f}x",
    f"scalar glv4 constant speedup: {speedup(scalar_basic, scalar_const):.3f}x",
    f"{avx_basic_name}: {avx_basic:.3f} {unit}",
    f"{avx_glv4_name}: {avx_glv4:.3f} {unit}",
    f"{avx_const_name}: {avx_const:.3f} {unit}",
    f"avx glv4 speedup: {speedup(avx_basic, avx_glv4):.3f}x",
    f"avx glv4 constant speedup: {speedup(avx_basic, avx_const):.3f}x",
    f"avx/scalar glv4 speedup ratio: {speedup(avx_basic, avx_glv4) / speedup(scalar_basic, scalar_glv4):.3f}x",
]

pathlib.Path(sys.argv[7]).write_text("\n".join(lines) + "\n")
print("\n".join(lines))
PY

echo "json reports written under ${LOCAL_REPORT_DIR}"
echo "summary written to ${SUMMARY_PATH}"

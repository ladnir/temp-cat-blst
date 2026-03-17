#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PLINK_EXE="${PLINK_EXE:-/c/Users/peter/OneDrive/tools/plink.exe}"
PLINK_KEY="${PLINK_KEY:-/c/Users/peter/OneDrive/keys/aws_2_19_26.ppk}"
REMOTE_HOST="${REMOTE_HOST:-18.188.231.170}"
REMOTE_USER="${REMOTE_USER:-ec2-user}"
REMOTE_HOSTKEY="${REMOTE_HOSTKEY:-ssh-ed25519 255 SHA256:/qd0qT37U/NAXgZT1kN6kSHaiyhoqQUXrlQBQXTa4Go}"
REMOTE_BASE="${REMOTE_BASE:-/home/ec2-user/vroom-parity}"
UPSTREAM_VROOM_REPO="${UPSTREAM_VROOM_REPO:-$ROOT_DIR/../VROOM}"
UPSTREAM_VROOM_COMMIT="${UPSTREAM_VROOM_COMMIT:-31e98bbf7f32e1345299d4276dd6a5bc44dc0b69}"
CLANG_ROOT="${CLANG_ROOT:-/opt/llvm-21/bin}"
PIN_CORE="${PIN_CORE:-0}"
REPETITIONS="${REPETITIONS:-5}"
MIN_TIME="${MIN_TIME:-0.05s}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
REMOTE_RUN_DIR="${REMOTE_BASE}/${RUN_ID}"
LOCAL_REPORT_DIR="${ROOT_DIR}/out/vroom-parity/${RUN_ID}"

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

mkdir -p "${LOCAL_REPORT_DIR}"

if [[ ! -f "${PLINK_EXE}" ]]; then
    echo "missing plink executable: ${PLINK_EXE}" >&2
    exit 1
fi

if [[ ! -d "${UPSTREAM_VROOM_REPO}" ]]; then
    echo "missing upstream VROOM repo: ${UPSTREAM_VROOM_REPO}" >&2
    exit 1
fi

PLINK_KEY_WIN="$(to_windows_path "${PLINK_KEY}")"

if ! git -C "${ROOT_DIR}" rev-parse --verify HEAD >/dev/null 2>&1; then
    echo "failed to resolve HEAD in ${ROOT_DIR}" >&2
    exit 1
fi

if ! git -C "${UPSTREAM_VROOM_REPO}" rev-parse --verify "${UPSTREAM_VROOM_COMMIT}" >/dev/null 2>&1; then
    echo "failed to resolve upstream VROOM commit ${UPSTREAM_VROOM_COMMIT}" >&2
    exit 1
fi

CURRENT_HEAD="$(git -C "${ROOT_DIR}" rev-parse HEAD)"
BLST_SNAPSHOT_MODE="committed HEAD ${CURRENT_HEAD}"
if [[ -n "$(git -C "${ROOT_DIR}" status --short)" ]]; then
    BLST_SNAPSHOT_MODE="current worktree over HEAD ${CURRENT_HEAD}"
    echo "note: local worktree is dirty; benchmarking ${BLST_SNAPSHOT_MODE}" >&2
fi

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

BLST_ARCHIVE="${TMP_DIR}/blst_current.tgz"
VROOM_ARCHIVE="${TMP_DIR}/VROOM_upstream.tgz"
BLST_SNAPSHOT_DIR="${TMP_DIR}/blst_snapshot"

mkdir -p "${BLST_SNAPSHOT_DIR}"
git -C "${ROOT_DIR}" archive HEAD | tar -xf - -C "${BLST_SNAPSHOT_DIR}"

while IFS= read -r -d '' path; do
    [[ "${path}" == "nul" ]] && continue
    if [[ -f "${ROOT_DIR}/${path}" ]]; then
        mkdir -p "${BLST_SNAPSHOT_DIR}/$(dirname "${path}")"
        cp "${ROOT_DIR}/${path}" "${BLST_SNAPSHOT_DIR}/${path}"
    else
        rm -f "${BLST_SNAPSHOT_DIR:?}/${path}"
    fi
done < <(git -C "${ROOT_DIR}" diff --name-only -z HEAD)

while IFS= read -r -d '' path; do
    [[ "${path}" == "nul" ]] && continue
    mkdir -p "${BLST_SNAPSHOT_DIR}/$(dirname "${path}")"
    cp "${ROOT_DIR}/${path}" "${BLST_SNAPSHOT_DIR}/${path}"
done < <(git -C "${ROOT_DIR}" ls-files --others --exclude-standard -z)

tar -czf "${BLST_ARCHIVE}" -C "${BLST_SNAPSHOT_DIR}" .
git -C "${UPSTREAM_VROOM_REPO}" archive --format=tar.gz --output="${VROOM_ARCHIVE}" "${UPSTREAM_VROOM_COMMIT}"

echo "remote host: ${REMOTE_USER}@${REMOTE_HOST}"
echo "blst snapshot: ${BLST_SNAPSHOT_MODE}"
echo "upstream VROOM commit: ${UPSTREAM_VROOM_COMMIT}"
echo "remote run dir: ${REMOTE_RUN_DIR}"
echo "local report dir: ${LOCAL_REPORT_DIR}"

remote "mkdir -p '${REMOTE_RUN_DIR}/archives'"
upload "${BLST_ARCHIVE}" "${REMOTE_RUN_DIR}/archives/blst_current.tgz"
upload "${VROOM_ARCHIVE}" "${REMOTE_RUN_DIR}/archives/VROOM_upstream.tgz"

remote "rm -rf '${REMOTE_RUN_DIR}/blst_current' '${REMOTE_RUN_DIR}/VROOM_upstream' && mkdir -p '${REMOTE_RUN_DIR}/blst_current' '${REMOTE_RUN_DIR}/VROOM_upstream' && tar -xzf '${REMOTE_RUN_DIR}/archives/blst_current.tgz' -C '${REMOTE_RUN_DIR}/blst_current' && tar -xzf '${REMOTE_RUN_DIR}/archives/VROOM_upstream.tgz' -C '${REMOTE_RUN_DIR}/VROOM_upstream'"

remote "cd '${REMOTE_RUN_DIR}/blst_current' && cmake -S . -B out/parity-clang21 -G Ninja -DCMAKE_C_COMPILER='${CLANG_ROOT}/clang' -DCMAKE_CXX_COMPILER='${CLANG_ROOT}/clang++' -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS_RELEASE='-O3 -march=native' -DCMAKE_CXX_FLAGS_RELEASE='-O3 -march=native' -DBLST_PORTABLE=OFF -DBLST_ASM=ON -DBLST_ADX=ON -DBLST_AVX512=ON -DBLST_DEFAULT_AVX=ON -DBLST_AVX512_PORTABLE=OFF -DBLST_BENCHMARK=ON -DBLST_TEST=OFF -DBLST_CUDA=OFF >/dev/null && cmake --build out/parity-clang21 --target blst_benchmarks -j2 >/dev/null"

remote "cd '${REMOTE_RUN_DIR}/VROOM_upstream/blst' && make clean >/dev/null 2>&1 || true && make -j2 CC='${CLANG_ROOT}/clang' CXX='${CLANG_ROOT}/clang++' AR=ar RANLIB=ranlib libblst_small.a >/dev/null"
remote "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && make clean >/dev/null 2>&1 || true && make -j2 CXX='${CLANG_ROOT}/clang++' INCLUDES='-I.. -I${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/_deps/googlebenchmark-src/include' BENCHMARK_LIBS='${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/_deps/googlebenchmark-build/src/libbenchmark.a -lpthread' bench_pairing >/dev/null"

UPSTREAM_MODMUL_JSON="$(run_remote_json "upstream_modmul" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_ModMul$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_JSON="$(run_remote_json "upstream_fp12_mul" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_Mul$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_CYCLO_JSON="$(run_remote_json "upstream_fp12_cyclotomic_sqr" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_CyclotomicSqr$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_MODMUL_JSON="$(run_remote_json "intree_modmul" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_ModMul$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_JSON="$(run_remote_json "intree_fp12_mul" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_Mul$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_CYCLO_JSON="$(run_remote_json "intree_fp12_cyclotomic_sqr" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_CyclotomicSqr$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
SCALAR_FP12_JSON="$(run_remote_json "scalar_fp12_mul" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_ScalarMulBaseline$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
SCALAR_FP12_CYCLO_JSON="$(run_remote_json "scalar_fp12_cyclotomic_sqr" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_ScalarCyclotomicSqrBaseline$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"

SUMMARY_PATH="${LOCAL_REPORT_DIR}/summary.txt"

python3 - "${UPSTREAM_MODMUL_JSON}" "${INTREE_MODMUL_JSON}" "${UPSTREAM_FP12_JSON}" "${INTREE_FP12_JSON}" "${UPSTREAM_FP12_CYCLO_JSON}" "${INTREE_FP12_CYCLO_JSON}" "${SCALAR_FP12_JSON}" "${SCALAR_FP12_CYCLO_JSON}" "${SUMMARY_PATH}" <<'PY'
import json
import pathlib
import sys

def load_median(path_str):
    data = json.loads(pathlib.Path(path_str).read_text())
    for bench in data["benchmarks"]:
        if bench.get("aggregate_name") == "median":
            return bench["name"], float(bench["cpu_time"]), bench["time_unit"]
    raise RuntimeError(f"median aggregate not found in {path_str}")

up_mod_name, up_mod, mod_unit = load_median(sys.argv[1])
in_mod_name, in_mod, _ = load_median(sys.argv[2])
up_fp12_name, up_fp12, fp12_unit = load_median(sys.argv[3])
in_fp12_name, in_fp12, _ = load_median(sys.argv[4])
up_cyclo_name, up_cyclo, cyclo_unit = load_median(sys.argv[5])
in_cyclo_name, in_cyclo, _ = load_median(sys.argv[6])
scalar_name, scalar_fp12, scalar_unit = load_median(sys.argv[7])
scalar_cyclo_name, scalar_cyclo, scalar_cyclo_unit = load_median(sys.argv[8])
summary_path = pathlib.Path(sys.argv[9])

def drift_percent(current, baseline):
    return ((current - baseline) / baseline) * 100.0

mod_drift = drift_percent(in_mod, up_mod)
fp12_drift = drift_percent(in_fp12, up_fp12)
cyclo_drift = drift_percent(in_cyclo, up_cyclo)

lines = [
    f"upstream {up_mod_name}: {up_mod:.3f} {mod_unit}",
    f"in-tree  {in_mod_name}: {in_mod:.3f} {mod_unit}",
    f"drift BM_ModMul: {mod_drift:+.3f}%",
    f"upstream {up_fp12_name}: {up_fp12:.3f} {fp12_unit}",
    f"in-tree  {in_fp12_name}: {in_fp12:.3f} {fp12_unit}",
    f"drift BM_FP12_Mul: {fp12_drift:+.3f}%",
    f"upstream {up_cyclo_name}: {up_cyclo:.3f} {cyclo_unit}",
    f"in-tree  {in_cyclo_name}: {in_cyclo:.3f} {cyclo_unit}",
    f"drift BM_FP12_CyclotomicSqr: {cyclo_drift:+.3f}%",
    f"scalar   {scalar_name}: {scalar_fp12:.3f} {scalar_unit}",
    f"scalar   {scalar_cyclo_name}: {scalar_cyclo:.3f} {scalar_cyclo_unit}",
    f"BM_ModMul parity gate: {'PASS' if abs(mod_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_Mul parity gate: {'PASS' if abs(fp12_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_CyclotomicSqr parity gate: {'PASS' if abs(cyclo_drift) <= 5.0 else 'FAIL'}",
]

summary_path.write_text("\n".join(lines) + "\n")
print("\n".join(lines))
PY

echo "json reports written under ${LOCAL_REPORT_DIR}"
echo "summary written to ${SUMMARY_PATH}"

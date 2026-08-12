#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PLINK_EXE="${PLINK_EXE:-/c/Users/peter/OneDrive/tools/plink.exe}"
PLINK_KEY="${PLINK_KEY:-/c/Users/peter/OneDrive/keys/aws_2_19_26.ppk}"
REMOTE_HOST="${REMOTE_HOST:-18.188.231.170}"
REMOTE_USER="${REMOTE_USER:-ec2-user}"
REMOTE_HOSTKEY="${REMOTE_HOSTKEY:-ssh-ed25519 255 SHA256:/qd0qT37U/NAXgZT1kN6kSHaiyhoqQUXrlQBQXTa4Go}"
REMOTE_BASE="${REMOTE_BASE:-/home/ec2-user/avx-fft}"
CLANG_ROOT="${CLANG_ROOT:-/opt/llvm-21/bin}"
PIN_CORE="${PIN_CORE:-0}"
REPETITIONS="${REPETITIONS:-5}"
MIN_TIME="${MIN_TIME:-0.05s}"
FFT_SIZES_REGEX="${FFT_SIZES_REGEX:-256|512|1024|2048}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
REMOTE_RUN_DIR="${REMOTE_BASE}/${RUN_ID}"
LOCAL_REPORT_DIR="${ROOT_DIR}/out/avx-fft/${RUN_ID}"

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

remote "cd '${REMOTE_RUN_DIR}/blst_current' && cmake -S . -B out/fft-clang21 -G Ninja -DCMAKE_C_COMPILER='${CLANG_ROOT}/clang' -DCMAKE_CXX_COMPILER='${CLANG_ROOT}/clang++' -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS_RELEASE='-O3 -march=native -DNDEBUG' -DCMAKE_CXX_FLAGS_RELEASE='-O3 -march=native -DNDEBUG' -DBLST_PORTABLE=OFF -DBLST_ASM=ON -DBLST_ADX=ON -DBLST_AVX512=ON -DBLST_DEFAULT_AVX=ON -DBLST_AVX512_PORTABLE=OFF -DBLST_BENCHMARK=ON -DBLST_TEST=ON -DBLST_CUDA=OFF >/dev/null && cmake --build out/fft-clang21 --target blst_benchmarks blst_tests -j2 >/dev/null"

remote "cd '${REMOTE_RUN_DIR}/blst_current' && ./out/fft-clang21/test/blst_tests --gtest_filter='Fft.GTAvxInverseFftMatchesScalar' >/dev/null"

SCALAR_JSON="$(run_remote_json "scalar_fft" "cd '${REMOTE_RUN_DIR}/blst_current/out/fft-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^bench_gt_ifft_cpu/(${FFT_SIZES_REGEX})/real_time$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
AVX_JSON="$(run_remote_json "avx_fft" "cd '${REMOTE_RUN_DIR}/blst_current/out/fft-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^bench_gt_ifft_avx52_precomputed_cpu/(${FFT_SIZES_REGEX})/real_time$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"

SUMMARY_PATH="${LOCAL_REPORT_DIR}/summary.txt"
python3 - "${SCALAR_JSON}" "${AVX_JSON}" "${SUMMARY_PATH}" <<'PY'
import json
import pathlib
import re
import sys

def load_medians(path_str):
    data = json.loads(pathlib.Path(path_str).read_text())
    medians = {}
    for bench in data["benchmarks"]:
        if bench.get("aggregate_name") != "median":
            continue
        medians[bench["name"]] = (float(bench["cpu_time"]), bench["time_unit"])
    return medians

scalar = load_medians(sys.argv[1])
avx = load_medians(sys.argv[2])

lines = []
for scalar_name, (scalar_time, unit) in sorted(scalar.items()):
    match = re.search(r"/(\\d+)$", scalar_name)
    if not match:
        continue
    size = match.group(1)
    avx_name = f"bench_gt_ifft_avx52_precomputed_cpu/{size}"
    if avx_name not in avx:
        continue
    avx_time, _ = avx[avx_name]
    drift = ((avx_time - scalar_time) / scalar_time) * 100.0
    speedup = scalar_time / avx_time
    lines.append(
        f"N={size}: scalar {scalar_time:.3f} {unit}, avx {avx_time:.3f} {unit}, drift {drift:+.3f}%, speedup {speedup:.3f}x"
    )

pathlib.Path(sys.argv[3]).write_text("\\n".join(lines) + "\\n")
print("\\n".join(lines))
PY

echo "json reports written under ${LOCAL_REPORT_DIR}"
echo "summary written to ${SUMMARY_PATH}"

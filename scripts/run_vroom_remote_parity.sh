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
remote "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && make -j2 CXX='${CLANG_ROOT}/clang++' INCLUDES='-I.. -I${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/_deps/googlebenchmark-src/include' BENCHMARK_LIBS='${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/_deps/googlebenchmark-build/src/libbenchmark.a -lpthread' bench_pairing_50bit >/dev/null"

UPSTREAM_MODMUL_JSON="$(run_remote_json "upstream_modmul" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_ModMul$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_INVERSE_JSON="$(run_remote_json "upstream_inverse" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_Inverse$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_1_JSON="$(run_remote_json "upstream_batch_reduce_expand_1" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_BatchReduceExpand_1$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_2_JSON="$(run_remote_json "upstream_batch_reduce_expand_2" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_BatchReduceExpand_2$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_3_JSON="$(run_remote_json "upstream_batch_reduce_expand_3" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_BatchReduceExpand_3$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_4_JSON="$(run_remote_json "upstream_batch_reduce_expand_4" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_BatchReduceExpand_4$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_5_JSON="$(run_remote_json "upstream_batch_reduce_expand_5" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_BatchReduceExpand_5$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_6_JSON="$(run_remote_json "upstream_batch_reduce_expand_6" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_BatchReduceExpand_6$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_12_JSON="$(run_remote_json "upstream_batch_reduce_expand_12" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_BatchReduceExpand_12$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G1_ADD_JSON="$(run_remote_json "upstream_g1_add" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_G1_Add$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G1_MIXED_ADD_JSON="$(run_remote_json "upstream_g1_mixed_add" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_G1_MixedAdd$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G1_DOUBLE_JSON="$(run_remote_json "upstream_g1_double" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_G1_Double$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G2_ADD_JSON="$(run_remote_json "upstream_g2_add" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_G2_Add$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G2_MIXED_ADD_JSON="$(run_remote_json "upstream_g2_mixed_add" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_G2_MixedAdd$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G2_DOUBLE_JSON="$(run_remote_json "upstream_g2_double" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_G2_Double$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G1_SCALAR_MULT_JSON="$(run_remote_json "upstream_g1_scalar_mult_glv" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_G1_ScalarMult_GLV$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G2_SCALAR_MULT_JSON="$(run_remote_json "upstream_g2_scalar_mult_gls" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_G2_ScalarMult_GLS$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_50_JSON="$(run_remote_json "upstream_fp12_mul_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP12_Mul$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_SQR_50_JSON="$(run_remote_json "upstream_fp12_sqr_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP12_Sqr$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_CYCLO_50_JSON="$(run_remote_json "upstream_fp12_cyclotomic_sqr_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP12_CyclotomicSqr$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_CONJ_50_JSON="$(run_remote_json "upstream_fp12_conjugate_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP12_Conjugate$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_FROB1_50_JSON="$(run_remote_json "upstream_fp12_frobenius_map1_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP12_FrobeniusMap1$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_FROB2_50_JSON="$(run_remote_json "upstream_fp12_frobenius_map2_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP12_FrobeniusMap2$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_FROB3_50_JSON="$(run_remote_json "upstream_fp12_frobenius_map3_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP12_FrobeniusMap3$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_FINALEXP_50_JSON="$(run_remote_json "upstream_fp12_final_exp_after_inverse_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP12_FinalExpAfterInverse$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_SPARSE_50_JSON="$(run_remote_json "upstream_fp12_mul_by_xy0z00_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP12_MulByXy0z00$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_INV_BLST_50_JSON="$(run_remote_json "upstream_fp12_inverse_blst_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP12_Inverse_BLST$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_INV_RNS_50_JSON="$(run_remote_json "upstream_fp12_inverse_rns_blst_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP12_Inverse_RNS_BLST$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_FINALEXP_BLST_50_JSON="$(run_remote_json "upstream_fp12_final_exp_blst_inverter_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP12_FinalExp_BLST_Inverter$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_FINALEXP_RNS_50_JSON="$(run_remote_json "upstream_fp12_final_exp_rns_blst_inverter_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP12_FinalExp_RNS_BLST_Inverter$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_MILLER_DOUBLE_50_JSON="$(run_remote_json "upstream_miller_double_step_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_MillerLoop_DoubleStep$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_MILLER_ADD_50_JSON="$(run_remote_json "upstream_miller_add_step_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_MillerLoop_AddStep$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_MILLER_50_JSON="$(run_remote_json "upstream_miller_loop_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_MillerLoop$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_PAIRING_BLST_50_JSON="$(run_remote_json "upstream_pairing_blst_inverter_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_Pairing_BLST_Inverter$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_PAIRING_RNS_50_JSON="$(run_remote_json "upstream_pairing_rns_blst_inverter_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_Pairing_RNS_BLST_Inverter$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G1_SCALAR_MULT_50_JSON="$(run_remote_json "upstream_g1_scalar_mult_glv_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_G1_ScalarMult_GLV_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G2_SCALAR_MULT_50_JSON="$(run_remote_json "upstream_g2_scalar_mult_gls_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_G2_ScalarMult_GLS_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_MODMUL_50_JSON="$(run_remote_json "upstream_modmul_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_ModMul$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_MODMUL_1_50_JSON="$(run_remote_json "upstream_batch_modmul_1_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchModMul_1$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_MODMUL_2_50_JSON="$(run_remote_json "upstream_batch_modmul_2_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchModMul_2$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_MODMUL_3_50_JSON="$(run_remote_json "upstream_batch_modmul_3_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchModMul_3$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_MODMUL_4_50_JSON="$(run_remote_json "upstream_batch_modmul_4_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchModMul_4$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_MODMUL_5_50_JSON="$(run_remote_json "upstream_batch_modmul_5_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchModMul_5$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_MODMUL_6_50_JSON="$(run_remote_json "upstream_batch_modmul_6_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchModMul_6$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_MODMUL_7_50_JSON="$(run_remote_json "upstream_batch_modmul_7_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchModMul_7$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_MODMUL_8_50_JSON="$(run_remote_json "upstream_batch_modmul_8_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchModMul_8$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_MODMUL_10_50_JSON="$(run_remote_json "upstream_batch_modmul_10_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchModMul_10$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_MODMUL_12_50_JSON="$(run_remote_json "upstream_batch_modmul_12_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchModMul_12$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_1_50_JSON="$(run_remote_json "upstream_batch_reduce_expand_1_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchReduceExpand_1$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_2_50_JSON="$(run_remote_json "upstream_batch_reduce_expand_2_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchReduceExpand_2$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_3_50_JSON="$(run_remote_json "upstream_batch_reduce_expand_3_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchReduceExpand_3$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_4_50_JSON="$(run_remote_json "upstream_batch_reduce_expand_4_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchReduceExpand_4$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_5_50_JSON="$(run_remote_json "upstream_batch_reduce_expand_5_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchReduceExpand_5$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_6_50_JSON="$(run_remote_json "upstream_batch_reduce_expand_6_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchReduceExpand_6$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_BATCH_REDUCE_12_50_JSON="$(run_remote_json "upstream_batch_reduce_expand_12_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_BatchReduceExpand_12$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP2_MUL_50_JSON="$(run_remote_json "upstream_fp2_mul_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP2_Mul$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP2_SQR_50_JSON="$(run_remote_json "upstream_fp2_sqr_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_FP2_Sqr$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_INVERSE_50_JSON="$(run_remote_json "upstream_inverse_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_Inverse$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G1_ADD_50_JSON="$(run_remote_json "upstream_g1_add_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_G1_Add$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G1_MIXED_ADD_50_JSON="$(run_remote_json "upstream_g1_mixed_add_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_G1_MixedAdd$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G1_DOUBLE_50_JSON="$(run_remote_json "upstream_g1_double_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_G1_Double$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G2_ADD_50_JSON="$(run_remote_json "upstream_g2_add_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_G2_Add$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G2_MIXED_ADD_50_JSON="$(run_remote_json "upstream_g2_mixed_add_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_G2_MixedAdd$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_G2_DOUBLE_50_JSON="$(run_remote_json "upstream_g2_double_50bit" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing_50bit --benchmark_filter='^BM_G2_Double$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_JSON="$(run_remote_json "upstream_fp12_mul" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_Mul$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_CYCLO_JSON="$(run_remote_json "upstream_fp12_cyclotomic_sqr" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_CyclotomicSqr$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_INV_BLST_JSON="$(run_remote_json "upstream_fp12_inverse_blst" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_Inverse_BLST$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_INV_RNS_JSON="$(run_remote_json "upstream_fp12_inverse_rns_blst" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_Inverse_RNS_BLST$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_CONJ_JSON="$(run_remote_json "upstream_fp12_conjugate" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_Conjugate$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_FROB1_JSON="$(run_remote_json "upstream_fp12_frobenius_map1" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_FrobeniusMap1$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_FROB2_JSON="$(run_remote_json "upstream_fp12_frobenius_map2" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_FrobeniusMap2$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_FROB3_JSON="$(run_remote_json "upstream_fp12_frobenius_map3" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_FrobeniusMap3$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_FINALEXP_JSON="$(run_remote_json "upstream_fp12_final_exp_after_inverse" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_FinalExpAfterInverse$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_FINALEXP_BLST_JSON="$(run_remote_json "upstream_fp12_final_exp_blst_inverter" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_FinalExp_BLST_Inverter$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_FINALEXP_RNS_JSON="$(run_remote_json "upstream_fp12_final_exp_rns_blst_inverter" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_FinalExp_RNS_BLST_Inverter$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_FP12_SPARSE_JSON="$(run_remote_json "upstream_fp12_mul_by_xy0z00" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_FP12_MulByXy0z00$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_MILLER_DOUBLE_JSON="$(run_remote_json "upstream_miller_double_step" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_MillerLoop_DoubleStep$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_MILLER_ADD_JSON="$(run_remote_json "upstream_miller_add_step" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_MillerLoop_AddStep$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_MILLER_JSON="$(run_remote_json "upstream_miller_loop" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_MillerLoop$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_PAIRING_BLST_JSON="$(run_remote_json "upstream_pairing_blst_inverter" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_Pairing_BLST_Inverter$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
UPSTREAM_PAIRING_RNS_JSON="$(run_remote_json "upstream_pairing_rns_blst_inverter" "cd '${REMOTE_RUN_DIR}/VROOM_upstream/src' && taskset -c ${PIN_CORE} ./bench_pairing --benchmark_filter='^BM_Pairing_RNS_BLST_Inverter$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_MODMUL_JSON="$(run_remote_json "intree_modmul" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_ModMul$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_INVERSE_JSON="$(run_remote_json "intree_inverse" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_Inverse$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_1_JSON="$(run_remote_json "intree_batch_reduce_expand_1" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_1$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_2_JSON="$(run_remote_json "intree_batch_reduce_expand_2" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_2$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_3_JSON="$(run_remote_json "intree_batch_reduce_expand_3" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_3$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_4_JSON="$(run_remote_json "intree_batch_reduce_expand_4" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_4$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_5_JSON="$(run_remote_json "intree_batch_reduce_expand_5" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_5$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_6_JSON="$(run_remote_json "intree_batch_reduce_expand_6" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_6$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_12_JSON="$(run_remote_json "intree_batch_reduce_expand_12" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_12$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G1_ADD_JSON="$(run_remote_json "intree_g1_add" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G1_Add$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G1_MIXED_ADD_JSON="$(run_remote_json "intree_g1_mixed_add" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G1_MixedAdd$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G1_DOUBLE_JSON="$(run_remote_json "intree_g1_double" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G1_Double$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G2_ADD_JSON="$(run_remote_json "intree_g2_add" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G2_Add$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G2_MIXED_ADD_JSON="$(run_remote_json "intree_g2_mixed_add" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G2_MixedAdd$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G2_DOUBLE_JSON="$(run_remote_json "intree_g2_double" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G2_Double$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G1_SCALAR_MULT_JSON="$(run_remote_json "intree_g1_scalar_mult_glv" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G1_ScalarMult_GLV$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G2_SCALAR_MULT_JSON="$(run_remote_json "intree_g2_scalar_mult_gls" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G2_ScalarMult_GLS$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_50_JSON="$(run_remote_json "intree_fp12_mul_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_Mul_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_SQR_50_JSON="$(run_remote_json "intree_fp12_sqr_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_Sqr_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_CYCLO_50_JSON="$(run_remote_json "intree_fp12_cyclotomic_sqr_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_CyclotomicSqr_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_CONJ_50_JSON="$(run_remote_json "intree_fp12_conjugate_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_Conjugate_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_FROB1_50_JSON="$(run_remote_json "intree_fp12_frobenius_map1_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_FrobeniusMap1_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_FROB2_50_JSON="$(run_remote_json "intree_fp12_frobenius_map2_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_FrobeniusMap2_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_FROB3_50_JSON="$(run_remote_json "intree_fp12_frobenius_map3_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_FrobeniusMap3_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_FINALEXP_50_JSON="$(run_remote_json "intree_fp12_final_exp_after_inverse_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_FinalExpAfterInverse_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_SPARSE_50_JSON="$(run_remote_json "intree_fp12_mul_by_xy0z00_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_MulByXy0z00_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_INV_BLST_50_JSON="$(run_remote_json "intree_fp12_inverse_blst_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_Inverse_BLST_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_INV_RNS_50_JSON="$(run_remote_json "intree_fp12_inverse_rns_blst_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_Inverse_RNS_BLST_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_FINALEXP_BLST_50_JSON="$(run_remote_json "intree_fp12_final_exp_blst_inverter_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_FinalExp_BLST_Inverter_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_FINALEXP_RNS_50_JSON="$(run_remote_json "intree_fp12_final_exp_rns_blst_inverter_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_FinalExp_RNS_BLST_Inverter_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_MILLER_DOUBLE_50_JSON="$(run_remote_json "intree_miller_double_step_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_MillerLoop_DoubleStep_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_MILLER_ADD_50_JSON="$(run_remote_json "intree_miller_add_step_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_MillerLoop_AddStep_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_MILLER_50_JSON="$(run_remote_json "intree_miller_loop_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_MillerLoop_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_PAIRING_BLST_50_JSON="$(run_remote_json "intree_pairing_blst_inverter_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_Pairing_BLST_Inverter_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_PAIRING_RNS_50_JSON="$(run_remote_json "intree_pairing_rns_blst_inverter_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_Pairing_RNS_BLST_Inverter_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G1_SCALAR_MULT_50_JSON="$(run_remote_json "intree_g1_scalar_mult_glv_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G1_ScalarMult_GLV_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G2_SCALAR_MULT_50_JSON="$(run_remote_json "intree_g2_scalar_mult_gls_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G2_ScalarMult_GLS_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_MODMUL_50_JSON="$(run_remote_json "intree_modmul_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_ModMul_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_MODMUL_1_50_JSON="$(run_remote_json "intree_batch_modmul_1_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchModMul_1_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_MODMUL_2_50_JSON="$(run_remote_json "intree_batch_modmul_2_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchModMul_2_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_MODMUL_3_50_JSON="$(run_remote_json "intree_batch_modmul_3_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchModMul_3_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_MODMUL_4_50_JSON="$(run_remote_json "intree_batch_modmul_4_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchModMul_4_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_MODMUL_5_50_JSON="$(run_remote_json "intree_batch_modmul_5_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchModMul_5_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_MODMUL_6_50_JSON="$(run_remote_json "intree_batch_modmul_6_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchModMul_6_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_MODMUL_7_50_JSON="$(run_remote_json "intree_batch_modmul_7_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchModMul_7_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_MODMUL_8_50_JSON="$(run_remote_json "intree_batch_modmul_8_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchModMul_8_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_MODMUL_10_50_JSON="$(run_remote_json "intree_batch_modmul_10_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchModMul_10_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_MODMUL_12_50_JSON="$(run_remote_json "intree_batch_modmul_12_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchModMul_12_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_1_50_JSON="$(run_remote_json "intree_batch_reduce_expand_1_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_1_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_2_50_JSON="$(run_remote_json "intree_batch_reduce_expand_2_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_2_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_3_50_JSON="$(run_remote_json "intree_batch_reduce_expand_3_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_3_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_4_50_JSON="$(run_remote_json "intree_batch_reduce_expand_4_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_4_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_5_50_JSON="$(run_remote_json "intree_batch_reduce_expand_5_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_5_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_6_50_JSON="$(run_remote_json "intree_batch_reduce_expand_6_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_6_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_BATCH_REDUCE_12_50_JSON="$(run_remote_json "intree_batch_reduce_expand_12_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BatchReduceExpand_12_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP2_MUL_50_JSON="$(run_remote_json "intree_fp2_mul_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP2_Mul_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP2_SQR_50_JSON="$(run_remote_json "intree_fp2_sqr_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP2_Sqr_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_INVERSE_50_JSON="$(run_remote_json "intree_inverse_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_Inverse_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G1_ADD_50_JSON="$(run_remote_json "intree_g1_add_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G1_Add_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G1_MIXED_ADD_50_JSON="$(run_remote_json "intree_g1_mixed_add_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G1_MixedAdd_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G1_DOUBLE_50_JSON="$(run_remote_json "intree_g1_double_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G1_Double_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G2_ADD_50_JSON="$(run_remote_json "intree_g2_add_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G2_Add_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G2_MIXED_ADD_50_JSON="$(run_remote_json "intree_g2_mixed_add_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G2_MixedAdd_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_G2_DOUBLE_50_JSON="$(run_remote_json "intree_g2_double_50bit" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_G2_Double_50bit$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_JSON="$(run_remote_json "intree_fp12_mul" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_Mul$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_CYCLO_JSON="$(run_remote_json "intree_fp12_cyclotomic_sqr" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_CyclotomicSqr$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_INV_BLST_JSON="$(run_remote_json "intree_fp12_inverse_blst" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_Inverse_BLST$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_INV_RNS_JSON="$(run_remote_json "intree_fp12_inverse_rns_blst" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_Inverse_RNS_BLST$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_CONJ_JSON="$(run_remote_json "intree_fp12_conjugate" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_Conjugate$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_FROB1_JSON="$(run_remote_json "intree_fp12_frobenius_map1" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_FrobeniusMap1$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_FROB2_JSON="$(run_remote_json "intree_fp12_frobenius_map2" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_FrobeniusMap2$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_FROB3_JSON="$(run_remote_json "intree_fp12_frobenius_map3" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_FrobeniusMap3$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_FINALEXP_JSON="$(run_remote_json "intree_fp12_final_exp_after_inverse" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_FinalExpAfterInverse$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_FINALEXP_BLST_JSON="$(run_remote_json "intree_fp12_final_exp_blst_inverter" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_FinalExp_BLST_Inverter$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_FINALEXP_RNS_JSON="$(run_remote_json "intree_fp12_final_exp_rns_blst_inverter" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_FinalExp_RNS_BLST_Inverter$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_FP12_SPARSE_JSON="$(run_remote_json "intree_fp12_mul_by_xy0z00" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_MulByXy0z00$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_MILLER_DOUBLE_JSON="$(run_remote_json "intree_miller_double_step" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_MillerLoop_DoubleStep$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_MILLER_ADD_JSON="$(run_remote_json "intree_miller_add_step" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_MillerLoop_AddStep$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_MILLER_JSON="$(run_remote_json "intree_miller_loop" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_MillerLoop$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_PAIRING_BLST_JSON="$(run_remote_json "intree_pairing_blst_inverter" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_Pairing_BLST_Inverter$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
INTREE_PAIRING_RNS_JSON="$(run_remote_json "intree_pairing_rns_blst_inverter" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_Pairing_RNS_BLST_Inverter$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
SCALAR_FP12_JSON="$(run_remote_json "scalar_fp12_mul" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_ScalarMulBaseline$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
SCALAR_FP12_CYCLO_JSON="$(run_remote_json "scalar_fp12_cyclotomic_sqr" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_ScalarCyclotomicSqrBaseline$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
SCALAR_FP12_SPARSE_JSON="$(run_remote_json "scalar_fp12_mul_by_xy0z00" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_FP12_ScalarMulByXy0z00Baseline$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"
SCALAR_PAIRING_JSON="$(run_remote_json "scalar_pairing_final_exp" "cd '${REMOTE_RUN_DIR}/blst_current/out/parity-clang21/benchmark' && taskset -c ${PIN_CORE} ./blst_benchmarks --benchmark_filter='^BM_BlstScalarPairingFinalExp$' --benchmark_repetitions=${REPETITIONS} --benchmark_report_aggregates_only=true --benchmark_min_time=${MIN_TIME} --benchmark_format=json")"

SUMMARY_PATH="${LOCAL_REPORT_DIR}/summary.txt"

python3 - "${UPSTREAM_MODMUL_JSON}" "${INTREE_MODMUL_JSON}" "${UPSTREAM_FP12_JSON}" "${INTREE_FP12_JSON}" "${UPSTREAM_FP12_CYCLO_JSON}" "${INTREE_FP12_CYCLO_JSON}" "${UPSTREAM_FP12_INV_BLST_JSON}" "${INTREE_FP12_INV_BLST_JSON}" "${UPSTREAM_FP12_INV_RNS_JSON}" "${INTREE_FP12_INV_RNS_JSON}" "${UPSTREAM_FP12_CONJ_JSON}" "${INTREE_FP12_CONJ_JSON}" "${UPSTREAM_FP12_FROB1_JSON}" "${INTREE_FP12_FROB1_JSON}" "${UPSTREAM_FP12_FROB2_JSON}" "${INTREE_FP12_FROB2_JSON}" "${UPSTREAM_FP12_FROB3_JSON}" "${INTREE_FP12_FROB3_JSON}" "${UPSTREAM_FP12_FINALEXP_JSON}" "${INTREE_FP12_FINALEXP_JSON}" "${UPSTREAM_FP12_FINALEXP_BLST_JSON}" "${INTREE_FP12_FINALEXP_BLST_JSON}" "${UPSTREAM_FP12_FINALEXP_RNS_JSON}" "${INTREE_FP12_FINALEXP_RNS_JSON}" "${UPSTREAM_FP12_SPARSE_JSON}" "${INTREE_FP12_SPARSE_JSON}" "${UPSTREAM_MILLER_DOUBLE_JSON}" "${INTREE_MILLER_DOUBLE_JSON}" "${UPSTREAM_MILLER_ADD_JSON}" "${INTREE_MILLER_ADD_JSON}" "${UPSTREAM_MILLER_JSON}" "${INTREE_MILLER_JSON}" "${UPSTREAM_PAIRING_BLST_JSON}" "${INTREE_PAIRING_BLST_JSON}" "${UPSTREAM_PAIRING_RNS_JSON}" "${INTREE_PAIRING_RNS_JSON}" "${SCALAR_FP12_JSON}" "${SCALAR_FP12_CYCLO_JSON}" "${SCALAR_FP12_SPARSE_JSON}" "${SCALAR_PAIRING_JSON}" "${SUMMARY_PATH}" <<'PY'
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
up_inv_blst_name, up_inv_blst, inv_blst_unit = load_median(sys.argv[7])
in_inv_blst_name, in_inv_blst, _ = load_median(sys.argv[8])
up_inv_rns_name, up_inv_rns, inv_rns_unit = load_median(sys.argv[9])
in_inv_rns_name, in_inv_rns, _ = load_median(sys.argv[10])
up_conj_name, up_conj, conj_unit = load_median(sys.argv[11])
in_conj_name, in_conj, _ = load_median(sys.argv[12])
up_frob1_name, up_frob1, frob1_unit = load_median(sys.argv[13])
in_frob1_name, in_frob1, _ = load_median(sys.argv[14])
up_frob2_name, up_frob2, frob2_unit = load_median(sys.argv[15])
in_frob2_name, in_frob2, _ = load_median(sys.argv[16])
up_frob3_name, up_frob3, frob3_unit = load_median(sys.argv[17])
in_frob3_name, in_frob3, _ = load_median(sys.argv[18])
up_finalexp_name, up_finalexp, finalexp_unit = load_median(sys.argv[19])
in_finalexp_name, in_finalexp, _ = load_median(sys.argv[20])
up_finalexp_blst_name, up_finalexp_blst, finalexp_blst_unit = load_median(sys.argv[21])
in_finalexp_blst_name, in_finalexp_blst, _ = load_median(sys.argv[22])
up_finalexp_rns_name, up_finalexp_rns, finalexp_rns_unit = load_median(sys.argv[23])
in_finalexp_rns_name, in_finalexp_rns, _ = load_median(sys.argv[24])
up_sparse_name, up_sparse, sparse_unit = load_median(sys.argv[25])
in_sparse_name, in_sparse, _ = load_median(sys.argv[26])
up_miller_double_name, up_miller_double, miller_double_unit = load_median(sys.argv[27])
in_miller_double_name, in_miller_double, _ = load_median(sys.argv[28])
up_miller_add_name, up_miller_add, miller_add_unit = load_median(sys.argv[29])
in_miller_add_name, in_miller_add, _ = load_median(sys.argv[30])
up_miller_name, up_miller, miller_unit = load_median(sys.argv[31])
in_miller_name, in_miller, _ = load_median(sys.argv[32])
up_pairing_blst_name, up_pairing_blst, pairing_blst_unit = load_median(sys.argv[33])
in_pairing_blst_name, in_pairing_blst, _ = load_median(sys.argv[34])
up_pairing_rns_name, up_pairing_rns, pairing_rns_unit = load_median(sys.argv[35])
in_pairing_rns_name, in_pairing_rns, _ = load_median(sys.argv[36])
scalar_name, scalar_fp12, scalar_unit = load_median(sys.argv[37])
scalar_cyclo_name, scalar_cyclo, scalar_cyclo_unit = load_median(sys.argv[38])
scalar_sparse_name, scalar_sparse, scalar_sparse_unit = load_median(sys.argv[39])
scalar_pairing_name, scalar_pairing, scalar_pairing_unit = load_median(sys.argv[40])
summary_path = pathlib.Path(sys.argv[41])

def drift_percent(current, baseline):
    return ((current - baseline) / baseline) * 100.0

mod_drift = drift_percent(in_mod, up_mod)
fp12_drift = drift_percent(in_fp12, up_fp12)
cyclo_drift = drift_percent(in_cyclo, up_cyclo)
inv_blst_drift = drift_percent(in_inv_blst, up_inv_blst)
inv_rns_drift = drift_percent(in_inv_rns, up_inv_rns)
conj_drift = drift_percent(in_conj, up_conj)
frob1_drift = drift_percent(in_frob1, up_frob1)
frob2_drift = drift_percent(in_frob2, up_frob2)
frob3_drift = drift_percent(in_frob3, up_frob3)
finalexp_drift = drift_percent(in_finalexp, up_finalexp)
finalexp_blst_drift = drift_percent(in_finalexp_blst, up_finalexp_blst)
finalexp_rns_drift = drift_percent(in_finalexp_rns, up_finalexp_rns)
sparse_drift = drift_percent(in_sparse, up_sparse)
miller_double_drift = drift_percent(in_miller_double, up_miller_double)
miller_add_drift = drift_percent(in_miller_add, up_miller_add)
miller_drift = drift_percent(in_miller, up_miller)
pairing_blst_drift = drift_percent(in_pairing_blst, up_pairing_blst)
pairing_rns_drift = drift_percent(in_pairing_rns, up_pairing_rns)

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
    f"upstream {up_inv_blst_name}: {up_inv_blst:.3f} {inv_blst_unit}",
    f"in-tree  {in_inv_blst_name}: {in_inv_blst:.3f} {inv_blst_unit}",
    f"drift BM_FP12_Inverse_BLST: {inv_blst_drift:+.3f}%",
    f"upstream {up_inv_rns_name}: {up_inv_rns:.3f} {inv_rns_unit}",
    f"in-tree  {in_inv_rns_name}: {in_inv_rns:.3f} {inv_rns_unit}",
    f"drift BM_FP12_Inverse_RNS_BLST: {inv_rns_drift:+.3f}%",
    f"upstream {up_conj_name}: {up_conj:.3f} {conj_unit}",
    f"in-tree  {in_conj_name}: {in_conj:.3f} {conj_unit}",
    f"drift BM_FP12_Conjugate: {conj_drift:+.3f}%",
    f"upstream {up_frob1_name}: {up_frob1:.3f} {frob1_unit}",
    f"in-tree  {in_frob1_name}: {in_frob1:.3f} {frob1_unit}",
    f"drift BM_FP12_FrobeniusMap1: {frob1_drift:+.3f}%",
    f"upstream {up_frob2_name}: {up_frob2:.3f} {frob2_unit}",
    f"in-tree  {in_frob2_name}: {in_frob2:.3f} {frob2_unit}",
    f"drift BM_FP12_FrobeniusMap2: {frob2_drift:+.3f}%",
    f"upstream {up_frob3_name}: {up_frob3:.3f} {frob3_unit}",
    f"in-tree  {in_frob3_name}: {in_frob3:.3f} {frob3_unit}",
    f"drift BM_FP12_FrobeniusMap3: {frob3_drift:+.3f}%",
    f"upstream {up_finalexp_name}: {up_finalexp:.3f} {finalexp_unit}",
    f"in-tree  {in_finalexp_name}: {in_finalexp:.3f} {finalexp_unit}",
    f"drift BM_FP12_FinalExpAfterInverse: {finalexp_drift:+.3f}%",
    f"upstream {up_finalexp_blst_name}: {up_finalexp_blst:.3f} {finalexp_blst_unit}",
    f"in-tree  {in_finalexp_blst_name}: {in_finalexp_blst:.3f} {finalexp_blst_unit}",
    f"drift BM_FP12_FinalExp_BLST_Inverter: {finalexp_blst_drift:+.3f}%",
    f"upstream {up_finalexp_rns_name}: {up_finalexp_rns:.3f} {finalexp_rns_unit}",
    f"in-tree  {in_finalexp_rns_name}: {in_finalexp_rns:.3f} {finalexp_rns_unit}",
    f"drift BM_FP12_FinalExp_RNS_BLST_Inverter: {finalexp_rns_drift:+.3f}%",
    f"upstream {up_sparse_name}: {up_sparse:.3f} {sparse_unit}",
    f"in-tree  {in_sparse_name}: {in_sparse:.3f} {sparse_unit}",
    f"drift BM_FP12_MulByXy0z00: {sparse_drift:+.3f}%",
    f"upstream {up_miller_double_name}: {up_miller_double:.3f} {miller_double_unit}",
    f"in-tree  {in_miller_double_name}: {in_miller_double:.3f} {miller_double_unit}",
    f"drift BM_MillerLoop_DoubleStep: {miller_double_drift:+.3f}%",
    f"upstream {up_miller_add_name}: {up_miller_add:.3f} {miller_add_unit}",
    f"in-tree  {in_miller_add_name}: {in_miller_add:.3f} {miller_add_unit}",
    f"drift BM_MillerLoop_AddStep: {miller_add_drift:+.3f}%",
    f"upstream {up_miller_name}: {up_miller:.3f} {miller_unit}",
    f"in-tree  {in_miller_name}: {in_miller:.3f} {miller_unit}",
    f"drift BM_MillerLoop: {miller_drift:+.3f}%",
    f"upstream {up_pairing_blst_name}: {up_pairing_blst:.3f} {pairing_blst_unit}",
    f"in-tree  {in_pairing_blst_name}: {in_pairing_blst:.3f} {pairing_blst_unit}",
    f"drift BM_Pairing_BLST_Inverter: {pairing_blst_drift:+.3f}%",
    f"upstream {up_pairing_rns_name}: {up_pairing_rns:.3f} {pairing_rns_unit}",
    f"in-tree  {in_pairing_rns_name}: {in_pairing_rns:.3f} {pairing_rns_unit}",
    f"drift BM_Pairing_RNS_BLST_Inverter: {pairing_rns_drift:+.3f}%",
    f"scalar   {scalar_name}: {scalar_fp12:.3f} {scalar_unit}",
    f"scalar   {scalar_cyclo_name}: {scalar_cyclo:.3f} {scalar_cyclo_unit}",
    f"scalar   {scalar_sparse_name}: {scalar_sparse:.3f} {scalar_sparse_unit}",
    f"scalar   {scalar_pairing_name}: {scalar_pairing:.3f} {scalar_pairing_unit}",
    f"BM_ModMul parity gate: {'PASS' if abs(mod_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_Mul parity gate: {'PASS' if abs(fp12_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_CyclotomicSqr parity gate: {'PASS' if abs(cyclo_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_Inverse_BLST parity gate: {'PASS' if abs(inv_blst_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_Inverse_RNS_BLST parity gate: {'PASS' if abs(inv_rns_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_Conjugate parity gate: {'PASS' if abs(conj_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_FrobeniusMap1 parity gate: {'PASS' if abs(frob1_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_FrobeniusMap2 parity gate: {'PASS' if abs(frob2_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_FrobeniusMap3 parity gate: {'PASS' if abs(frob3_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_FinalExpAfterInverse parity gate: {'PASS' if abs(finalexp_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_FinalExp_BLST_Inverter parity gate: {'PASS' if abs(finalexp_blst_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_FinalExp_RNS_BLST_Inverter parity gate: {'PASS' if abs(finalexp_rns_drift) <= 5.0 else 'FAIL'}",
    f"BM_FP12_MulByXy0z00 parity gate: {'PASS' if abs(sparse_drift) <= 5.0 else 'FAIL'}",
    f"BM_MillerLoop_DoubleStep parity gate: {'PASS' if abs(miller_double_drift) <= 5.0 else 'FAIL'}",
    f"BM_MillerLoop_AddStep parity gate: {'PASS' if abs(miller_add_drift) <= 5.0 else 'FAIL'}",
    f"BM_MillerLoop parity gate: {'PASS' if abs(miller_drift) <= 5.0 else 'FAIL'}",
    f"BM_Pairing_BLST_Inverter parity gate: {'PASS' if abs(pairing_blst_drift) <= 5.0 else 'FAIL'}",
    f"BM_Pairing_RNS_BLST_Inverter parity gate: {'PASS' if abs(pairing_rns_drift) <= 5.0 else 'FAIL'}",
]

summary_path.write_text("\n".join(lines) + "\n")
print("\n".join(lines))
PY

python3 - "${UPSTREAM_INVERSE_JSON}" "${INTREE_INVERSE_JSON}" "${UPSTREAM_BATCH_REDUCE_1_JSON}" "${INTREE_BATCH_REDUCE_1_JSON}" "${UPSTREAM_BATCH_REDUCE_2_JSON}" "${INTREE_BATCH_REDUCE_2_JSON}" "${UPSTREAM_BATCH_REDUCE_3_JSON}" "${INTREE_BATCH_REDUCE_3_JSON}" "${UPSTREAM_BATCH_REDUCE_4_JSON}" "${INTREE_BATCH_REDUCE_4_JSON}" "${UPSTREAM_BATCH_REDUCE_5_JSON}" "${INTREE_BATCH_REDUCE_5_JSON}" "${UPSTREAM_BATCH_REDUCE_6_JSON}" "${INTREE_BATCH_REDUCE_6_JSON}" "${UPSTREAM_BATCH_REDUCE_12_JSON}" "${INTREE_BATCH_REDUCE_12_JSON}" "${SUMMARY_PATH}" <<'PY'
import json
import pathlib
import sys

def load_median(path_str):
    data = json.loads(pathlib.Path(path_str).read_text())
    for bench in data["benchmarks"]:
        if bench.get("aggregate_name") == "median":
            return bench["name"], float(bench["cpu_time"]), bench["time_unit"]
    raise RuntimeError(f"median aggregate not found in {path_str}")

def drift_percent(current, baseline):
    return ((current - baseline) / baseline) * 100.0

summary_path = pathlib.Path(sys.argv[17])
existing = summary_path.read_text()

entries = [
    ("BM_Inverse",) + load_median(sys.argv[1]) + load_median(sys.argv[2])[1:],
    ("BM_BatchReduceExpand_1",) + load_median(sys.argv[3]) + load_median(sys.argv[4])[1:],
    ("BM_BatchReduceExpand_2",) + load_median(sys.argv[5]) + load_median(sys.argv[6])[1:],
    ("BM_BatchReduceExpand_3",) + load_median(sys.argv[7]) + load_median(sys.argv[8])[1:],
    ("BM_BatchReduceExpand_4",) + load_median(sys.argv[9]) + load_median(sys.argv[10])[1:],
    ("BM_BatchReduceExpand_5",) + load_median(sys.argv[11]) + load_median(sys.argv[12])[1:],
    ("BM_BatchReduceExpand_6",) + load_median(sys.argv[13]) + load_median(sys.argv[14])[1:],
    ("BM_BatchReduceExpand_12",) + load_median(sys.argv[15]) + load_median(sys.argv[16])[1:],
]

lines = []
for benchmark_name, upstream_name, upstream_time, unit, intree_time, _ in entries:
    drift = drift_percent(intree_time, upstream_time)
    lines.extend([
        f"upstream {upstream_name}: {upstream_time:.3f} {unit}",
        f"in-tree  {benchmark_name}_median: {intree_time:.3f} {unit}",
        f"drift {benchmark_name}: {drift:+.3f}%",
        f"{benchmark_name} parity gate: {'PASS' if abs(drift) <= 5.0 else 'FAIL'}",
    ])

summary_path.write_text(existing + "\n".join(lines) + "\n")
print("\n".join(lines))
PY

python3 - "${UPSTREAM_G1_ADD_JSON}" "${INTREE_G1_ADD_JSON}" "${UPSTREAM_G1_MIXED_ADD_JSON}" "${INTREE_G1_MIXED_ADD_JSON}" "${UPSTREAM_G1_DOUBLE_JSON}" "${INTREE_G1_DOUBLE_JSON}" "${UPSTREAM_G2_ADD_JSON}" "${INTREE_G2_ADD_JSON}" "${UPSTREAM_G2_MIXED_ADD_JSON}" "${INTREE_G2_MIXED_ADD_JSON}" "${UPSTREAM_G2_DOUBLE_JSON}" "${INTREE_G2_DOUBLE_JSON}" "${SUMMARY_PATH}" <<'PY'
import json
import pathlib
import sys

def load_median(path_str):
    data = json.loads(pathlib.Path(path_str).read_text())
    for bench in data["benchmarks"]:
        if bench.get("aggregate_name") == "median":
            return bench["name"], float(bench["cpu_time"]), bench["time_unit"]
    raise RuntimeError(f"median aggregate not found in {path_str}")

def drift_percent(current, baseline):
    return ((current - baseline) / baseline) * 100.0

summary_path = pathlib.Path(sys.argv[13])
existing = summary_path.read_text()

entries = [
    ("BM_G1_Add",) + load_median(sys.argv[1]) + load_median(sys.argv[2])[1:],
    ("BM_G1_MixedAdd",) + load_median(sys.argv[3]) + load_median(sys.argv[4])[1:],
    ("BM_G1_Double",) + load_median(sys.argv[5]) + load_median(sys.argv[6])[1:],
    ("BM_G2_Add",) + load_median(sys.argv[7]) + load_median(sys.argv[8])[1:],
    ("BM_G2_MixedAdd",) + load_median(sys.argv[9]) + load_median(sys.argv[10])[1:],
    ("BM_G2_Double",) + load_median(sys.argv[11]) + load_median(sys.argv[12])[1:],
]

lines = []
for benchmark_name, upstream_name, upstream_time, unit, intree_time, _ in entries:
    drift = drift_percent(intree_time, upstream_time)
    lines.extend([
        f"upstream {upstream_name}: {upstream_time:.3f} {unit}",
        f"in-tree  {benchmark_name}_median: {intree_time:.3f} {unit}",
        f"drift {benchmark_name}: {drift:+.3f}%",
        f"{benchmark_name} parity gate: {'PASS' if abs(drift) <= 5.0 else 'FAIL'}",
    ])

summary_path.write_text(existing + "\n".join(lines) + "\n")
print("\n".join(lines))
PY

python3 - "${UPSTREAM_G1_SCALAR_MULT_JSON}" "${INTREE_G1_SCALAR_MULT_JSON}" "${UPSTREAM_G2_SCALAR_MULT_JSON}" "${INTREE_G2_SCALAR_MULT_JSON}" "${SUMMARY_PATH}" <<'PY'
import json
import pathlib
import sys

def load_median(path_str):
    data = json.loads(pathlib.Path(path_str).read_text())
    for bench in data["benchmarks"]:
        if bench.get("aggregate_name") == "median":
            return bench["name"], float(bench["cpu_time"]), bench["time_unit"]
    raise RuntimeError(f"median aggregate not found in {path_str}")

def drift_percent(current, baseline):
    return ((current - baseline) / baseline) * 100.0

summary_path = pathlib.Path(sys.argv[5])
existing = summary_path.read_text()

entries = [
    ("BM_G1_ScalarMult_GLV",) + load_median(sys.argv[1]) + load_median(sys.argv[2])[1:],
    ("BM_G2_ScalarMult_GLS",) + load_median(sys.argv[3]) + load_median(sys.argv[4])[1:],
]

lines = []
for benchmark_name, upstream_name, upstream_time, unit, intree_time, _ in entries:
    drift = drift_percent(intree_time, upstream_time)
    lines.extend([
        f"upstream {upstream_name}: {upstream_time:.3f} {unit}",
        f"in-tree  {benchmark_name}_median: {intree_time:.3f} {unit}",
        f"drift {benchmark_name}: {drift:+.3f}%",
        f"{benchmark_name} parity gate: {'PASS' if abs(drift) <= 5.0 else 'FAIL'}",
    ])

summary_path.write_text(existing + "\n".join(lines) + "\n")
print("\n".join(lines))
PY

python3 - "${UPSTREAM_G1_SCALAR_MULT_50_JSON}" "${INTREE_G1_SCALAR_MULT_50_JSON}" "${UPSTREAM_G2_SCALAR_MULT_50_JSON}" "${INTREE_G2_SCALAR_MULT_50_JSON}" "${SUMMARY_PATH}" <<'PY'
import json
import pathlib
import sys

def load_median(path_str):
    data = json.loads(pathlib.Path(path_str).read_text())
    for bench in data["benchmarks"]:
        if bench.get("aggregate_name") == "median":
            return bench["name"], float(bench["cpu_time"]), bench["time_unit"]
    raise RuntimeError(f"median aggregate not found in {path_str}")

def drift_percent(current, baseline):
    return ((current - baseline) / baseline) * 100.0

summary_path = pathlib.Path(sys.argv[5])
existing = summary_path.read_text()

entries = [
    ("BM_G1_ScalarMult_GLV_50bit",) + load_median(sys.argv[1]) + load_median(sys.argv[2])[1:],
    ("BM_G2_ScalarMult_GLS_50bit",) + load_median(sys.argv[3]) + load_median(sys.argv[4])[1:],
]

lines = []
for benchmark_name, upstream_name, upstream_time, unit, intree_time, _ in entries:
    drift = drift_percent(intree_time, upstream_time)
    lines.extend([
        f"upstream {upstream_name}: {upstream_time:.3f} {unit}",
        f"in-tree  {benchmark_name}_median: {intree_time:.3f} {unit}",
        f"drift {benchmark_name}: {drift:+.3f}%",
        f"{benchmark_name} parity gate: {'PASS' if abs(drift) <= 5.0 else 'FAIL'}",
    ])

summary_path.write_text(existing + "\n".join(lines) + "\n")
print("\n".join(lines))
PY

python3 - "${UPSTREAM_MODMUL_50_JSON}" "${INTREE_MODMUL_50_JSON}" "${UPSTREAM_BATCH_MODMUL_1_50_JSON}" "${INTREE_BATCH_MODMUL_1_50_JSON}" "${UPSTREAM_BATCH_MODMUL_2_50_JSON}" "${INTREE_BATCH_MODMUL_2_50_JSON}" "${UPSTREAM_BATCH_MODMUL_3_50_JSON}" "${INTREE_BATCH_MODMUL_3_50_JSON}" "${UPSTREAM_BATCH_MODMUL_4_50_JSON}" "${INTREE_BATCH_MODMUL_4_50_JSON}" "${UPSTREAM_BATCH_MODMUL_5_50_JSON}" "${INTREE_BATCH_MODMUL_5_50_JSON}" "${UPSTREAM_BATCH_MODMUL_6_50_JSON}" "${INTREE_BATCH_MODMUL_6_50_JSON}" "${UPSTREAM_BATCH_MODMUL_7_50_JSON}" "${INTREE_BATCH_MODMUL_7_50_JSON}" "${UPSTREAM_BATCH_MODMUL_8_50_JSON}" "${INTREE_BATCH_MODMUL_8_50_JSON}" "${UPSTREAM_BATCH_MODMUL_10_50_JSON}" "${INTREE_BATCH_MODMUL_10_50_JSON}" "${UPSTREAM_BATCH_MODMUL_12_50_JSON}" "${INTREE_BATCH_MODMUL_12_50_JSON}" "${UPSTREAM_BATCH_REDUCE_1_50_JSON}" "${INTREE_BATCH_REDUCE_1_50_JSON}" "${UPSTREAM_BATCH_REDUCE_2_50_JSON}" "${INTREE_BATCH_REDUCE_2_50_JSON}" "${UPSTREAM_BATCH_REDUCE_3_50_JSON}" "${INTREE_BATCH_REDUCE_3_50_JSON}" "${UPSTREAM_BATCH_REDUCE_4_50_JSON}" "${INTREE_BATCH_REDUCE_4_50_JSON}" "${UPSTREAM_BATCH_REDUCE_5_50_JSON}" "${INTREE_BATCH_REDUCE_5_50_JSON}" "${UPSTREAM_BATCH_REDUCE_6_50_JSON}" "${INTREE_BATCH_REDUCE_6_50_JSON}" "${UPSTREAM_BATCH_REDUCE_12_50_JSON}" "${INTREE_BATCH_REDUCE_12_50_JSON}" "${UPSTREAM_FP2_MUL_50_JSON}" "${INTREE_FP2_MUL_50_JSON}" "${UPSTREAM_FP2_SQR_50_JSON}" "${INTREE_FP2_SQR_50_JSON}" "${UPSTREAM_INVERSE_50_JSON}" "${INTREE_INVERSE_50_JSON}" "${SUMMARY_PATH}" <<'PY'
import json
import pathlib
import sys

def load_median(path_str):
    data = json.loads(pathlib.Path(path_str).read_text())
    for bench in data["benchmarks"]:
        if bench.get("aggregate_name") == "median":
            return bench["name"], float(bench["cpu_time"]), bench["time_unit"]
    raise RuntimeError(f"median aggregate not found in {path_str}")

def drift_percent(current, baseline):
    return ((current - baseline) / baseline) * 100.0

summary_path = pathlib.Path(sys.argv[43])
existing = summary_path.read_text()

entries = [
    ("BM_ModMul_50bit",) + load_median(sys.argv[1]) + load_median(sys.argv[2])[1:],
    ("BM_BatchModMul_1_50bit",) + load_median(sys.argv[3]) + load_median(sys.argv[4])[1:],
    ("BM_BatchModMul_2_50bit",) + load_median(sys.argv[5]) + load_median(sys.argv[6])[1:],
    ("BM_BatchModMul_3_50bit",) + load_median(sys.argv[7]) + load_median(sys.argv[8])[1:],
    ("BM_BatchModMul_4_50bit",) + load_median(sys.argv[9]) + load_median(sys.argv[10])[1:],
    ("BM_BatchModMul_5_50bit",) + load_median(sys.argv[11]) + load_median(sys.argv[12])[1:],
    ("BM_BatchModMul_6_50bit",) + load_median(sys.argv[13]) + load_median(sys.argv[14])[1:],
    ("BM_BatchModMul_7_50bit",) + load_median(sys.argv[15]) + load_median(sys.argv[16])[1:],
    ("BM_BatchModMul_8_50bit",) + load_median(sys.argv[17]) + load_median(sys.argv[18])[1:],
    ("BM_BatchModMul_10_50bit",) + load_median(sys.argv[19]) + load_median(sys.argv[20])[1:],
    ("BM_BatchModMul_12_50bit",) + load_median(sys.argv[21]) + load_median(sys.argv[22])[1:],
    ("BM_BatchReduceExpand_1_50bit",) + load_median(sys.argv[23]) + load_median(sys.argv[24])[1:],
    ("BM_BatchReduceExpand_2_50bit",) + load_median(sys.argv[25]) + load_median(sys.argv[26])[1:],
    ("BM_BatchReduceExpand_3_50bit",) + load_median(sys.argv[27]) + load_median(sys.argv[28])[1:],
    ("BM_BatchReduceExpand_4_50bit",) + load_median(sys.argv[29]) + load_median(sys.argv[30])[1:],
    ("BM_BatchReduceExpand_5_50bit",) + load_median(sys.argv[31]) + load_median(sys.argv[32])[1:],
    ("BM_BatchReduceExpand_6_50bit",) + load_median(sys.argv[33]) + load_median(sys.argv[34])[1:],
    ("BM_BatchReduceExpand_12_50bit",) + load_median(sys.argv[35]) + load_median(sys.argv[36])[1:],
    ("BM_FP2_Mul_50bit",) + load_median(sys.argv[37]) + load_median(sys.argv[38])[1:],
    ("BM_FP2_Sqr_50bit",) + load_median(sys.argv[39]) + load_median(sys.argv[40])[1:],
    ("BM_Inverse_50bit",) + load_median(sys.argv[41]) + load_median(sys.argv[42])[1:],
]

lines = []
for benchmark_name, upstream_name, upstream_time, unit, intree_time, _ in entries:
    drift = drift_percent(intree_time, upstream_time)
    lines.extend([
        f"upstream {upstream_name}: {upstream_time:.3f} {unit}",
        f"in-tree  {benchmark_name}_median: {intree_time:.3f} {unit}",
        f"drift {benchmark_name}: {drift:+.3f}%",
        f"{benchmark_name} parity gate: {'PASS' if abs(drift) <= 5.0 else 'FAIL'}",
    ])

summary_path.write_text(existing + "\n".join(lines) + "\n")
print("\n".join(lines))
PY

python3 - "${UPSTREAM_G1_ADD_50_JSON}" "${INTREE_G1_ADD_50_JSON}" "${UPSTREAM_G1_MIXED_ADD_50_JSON}" "${INTREE_G1_MIXED_ADD_50_JSON}" "${UPSTREAM_G1_DOUBLE_50_JSON}" "${INTREE_G1_DOUBLE_50_JSON}" "${UPSTREAM_G2_ADD_50_JSON}" "${INTREE_G2_ADD_50_JSON}" "${UPSTREAM_G2_MIXED_ADD_50_JSON}" "${INTREE_G2_MIXED_ADD_50_JSON}" "${UPSTREAM_G2_DOUBLE_50_JSON}" "${INTREE_G2_DOUBLE_50_JSON}" "${SUMMARY_PATH}" <<'PY'
import json
import pathlib
import sys

def load_median(path_str):
    data = json.loads(pathlib.Path(path_str).read_text())
    for bench in data["benchmarks"]:
        if bench.get("aggregate_name") == "median":
            return bench["name"], float(bench["cpu_time"]), bench["time_unit"]
    raise RuntimeError(f"median aggregate not found in {path_str}")

def drift_percent(current, baseline):
    return ((current - baseline) / baseline) * 100.0

summary_path = pathlib.Path(sys.argv[13])
existing = summary_path.read_text()

entries = [
    ("BM_G1_Add_50bit",) + load_median(sys.argv[1]) + load_median(sys.argv[2])[1:],
    ("BM_G1_MixedAdd_50bit",) + load_median(sys.argv[3]) + load_median(sys.argv[4])[1:],
    ("BM_G1_Double_50bit",) + load_median(sys.argv[5]) + load_median(sys.argv[6])[1:],
    ("BM_G2_Add_50bit",) + load_median(sys.argv[7]) + load_median(sys.argv[8])[1:],
    ("BM_G2_MixedAdd_50bit",) + load_median(sys.argv[9]) + load_median(sys.argv[10])[1:],
    ("BM_G2_Double_50bit",) + load_median(sys.argv[11]) + load_median(sys.argv[12])[1:],
]

lines = []
for benchmark_name, upstream_name, upstream_time, unit, intree_time, _ in entries:
    drift = drift_percent(intree_time, upstream_time)
    lines.extend([
        f"upstream {upstream_name}: {upstream_time:.3f} {unit}",
        f"in-tree  {benchmark_name}_median: {intree_time:.3f} {unit}",
        f"drift {benchmark_name}: {drift:+.3f}%",
        f"{benchmark_name} parity gate: {'PASS' if abs(drift) <= 5.0 else 'FAIL'}",
    ])

summary_path.write_text(existing + "\n".join(lines) + "\n")
print("\n".join(lines))
PY

python3 - "${UPSTREAM_FP12_50_JSON}" "${INTREE_FP12_50_JSON}" "${UPSTREAM_FP12_SQR_50_JSON}" "${INTREE_FP12_SQR_50_JSON}" "${UPSTREAM_FP12_CYCLO_50_JSON}" "${INTREE_FP12_CYCLO_50_JSON}" "${UPSTREAM_FP12_CONJ_50_JSON}" "${INTREE_FP12_CONJ_50_JSON}" "${UPSTREAM_FP12_FROB1_50_JSON}" "${INTREE_FP12_FROB1_50_JSON}" "${UPSTREAM_FP12_FROB2_50_JSON}" "${INTREE_FP12_FROB2_50_JSON}" "${UPSTREAM_FP12_FROB3_50_JSON}" "${INTREE_FP12_FROB3_50_JSON}" "${UPSTREAM_FP12_FINALEXP_50_JSON}" "${INTREE_FP12_FINALEXP_50_JSON}" "${UPSTREAM_FP12_SPARSE_50_JSON}" "${INTREE_FP12_SPARSE_50_JSON}" "${UPSTREAM_FP12_INV_BLST_50_JSON}" "${INTREE_FP12_INV_BLST_50_JSON}" "${UPSTREAM_FP12_INV_RNS_50_JSON}" "${INTREE_FP12_INV_RNS_50_JSON}" "${UPSTREAM_FP12_FINALEXP_BLST_50_JSON}" "${INTREE_FP12_FINALEXP_BLST_50_JSON}" "${UPSTREAM_FP12_FINALEXP_RNS_50_JSON}" "${INTREE_FP12_FINALEXP_RNS_50_JSON}" "${UPSTREAM_PAIRING_BLST_50_JSON}" "${INTREE_PAIRING_BLST_50_JSON}" "${UPSTREAM_PAIRING_RNS_50_JSON}" "${INTREE_PAIRING_RNS_50_JSON}" "${SUMMARY_PATH}" <<'PY'
import json
import pathlib
import sys

def load_median(path_str):
    data = json.loads(pathlib.Path(path_str).read_text())
    for bench in data["benchmarks"]:
        if bench.get("aggregate_name") == "median":
            return bench["name"], float(bench["cpu_time"]), bench["time_unit"]
    raise RuntimeError(f"median aggregate not found in {path_str}")

def drift_percent(current, baseline):
    return ((current - baseline) / baseline) * 100.0

summary_path = pathlib.Path(sys.argv[31])
existing = summary_path.read_text()

entries = [
    ("BM_FP12_Mul_50bit",) + load_median(sys.argv[1]) + load_median(sys.argv[2])[1:],
    ("BM_FP12_Sqr_50bit",) + load_median(sys.argv[3]) + load_median(sys.argv[4])[1:],
    ("BM_FP12_CyclotomicSqr_50bit",) + load_median(sys.argv[5]) + load_median(sys.argv[6])[1:],
    ("BM_FP12_Conjugate_50bit",) + load_median(sys.argv[7]) + load_median(sys.argv[8])[1:],
    ("BM_FP12_FrobeniusMap1_50bit",) + load_median(sys.argv[9]) + load_median(sys.argv[10])[1:],
    ("BM_FP12_FrobeniusMap2_50bit",) + load_median(sys.argv[11]) + load_median(sys.argv[12])[1:],
    ("BM_FP12_FrobeniusMap3_50bit",) + load_median(sys.argv[13]) + load_median(sys.argv[14])[1:],
    ("BM_FP12_FinalExpAfterInverse_50bit",) + load_median(sys.argv[15]) + load_median(sys.argv[16])[1:],
    ("BM_FP12_MulByXy0z00_50bit",) + load_median(sys.argv[17]) + load_median(sys.argv[18])[1:],
    ("BM_FP12_Inverse_BLST_50bit",) + load_median(sys.argv[19]) + load_median(sys.argv[20])[1:],
    ("BM_FP12_Inverse_RNS_BLST_50bit",) + load_median(sys.argv[21]) + load_median(sys.argv[22])[1:],
    ("BM_FP12_FinalExp_BLST_Inverter_50bit",) + load_median(sys.argv[23]) + load_median(sys.argv[24])[1:],
    ("BM_FP12_FinalExp_RNS_BLST_Inverter_50bit",) + load_median(sys.argv[25]) + load_median(sys.argv[26])[1:],
    ("BM_Pairing_BLST_Inverter_50bit",) + load_median(sys.argv[27]) + load_median(sys.argv[28])[1:],
    ("BM_Pairing_RNS_BLST_Inverter_50bit",) + load_median(sys.argv[29]) + load_median(sys.argv[30])[1:],
]

lines = []
for benchmark_name, upstream_name, upstream_time, unit, intree_time, _ in entries:
    drift = drift_percent(intree_time, upstream_time)
    lines.extend([
        f"upstream {upstream_name}: {upstream_time:.3f} {unit}",
        f"in-tree  {benchmark_name}_median: {intree_time:.3f} {unit}",
        f"drift {benchmark_name}: {drift:+.3f}%",
        f"{benchmark_name} parity gate: {'PASS' if abs(drift) <= 5.0 else 'FAIL'}",
    ])

summary_path.write_text(existing + "\n".join(lines) + "\n")
print("\n".join(lines))
PY

python3 - "${UPSTREAM_MILLER_DOUBLE_50_JSON}" "${INTREE_MILLER_DOUBLE_50_JSON}" "${UPSTREAM_MILLER_ADD_50_JSON}" "${INTREE_MILLER_ADD_50_JSON}" "${UPSTREAM_MILLER_50_JSON}" "${INTREE_MILLER_50_JSON}" "${SUMMARY_PATH}" <<'PY'
import json
import pathlib
import sys

def load_median(path_str):
    data = json.loads(pathlib.Path(path_str).read_text())
    for bench in data["benchmarks"]:
        if bench.get("aggregate_name") == "median":
            return bench["name"], float(bench["cpu_time"]), bench["time_unit"]
    raise RuntimeError(f"median aggregate not found in {path_str}")

def drift_percent(current, baseline):
    return ((current - baseline) / baseline) * 100.0

summary_path = pathlib.Path(sys.argv[7])
existing = summary_path.read_text()

entries = [
    ("BM_MillerLoop_DoubleStep_50bit",) + load_median(sys.argv[1]) + load_median(sys.argv[2])[1:],
    ("BM_MillerLoop_AddStep_50bit",) + load_median(sys.argv[3]) + load_median(sys.argv[4])[1:],
    ("BM_MillerLoop_50bit",) + load_median(sys.argv[5]) + load_median(sys.argv[6])[1:],
]

lines = []
for benchmark_name, upstream_name, upstream_time, unit, intree_time, _ in entries:
    drift = drift_percent(intree_time, upstream_time)
    lines.extend([
        f"upstream {upstream_name}: {upstream_time:.3f} {unit}",
        f"in-tree  {benchmark_name}_median: {intree_time:.3f} {unit}",
        f"drift {benchmark_name}: {drift:+.3f}%",
        f"{benchmark_name} parity gate: {'PASS' if abs(drift) <= 5.0 else 'FAIL'}",
    ])

summary_path.write_text(existing + "\n".join(lines) + "\n")
print("\n".join(lines))
PY

echo "json reports written under ${LOCAL_REPORT_DIR}"
echo "summary written to ${SUMMARY_PATH}"

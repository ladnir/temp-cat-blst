#!/usr/bin/env bash
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

blst_fuzz_target_limits() {
    case "$1" in
        blst_fuzz_decode) BLST_FUZZ_MAX_LEN=192; BLST_FUZZ_TIMEOUT=2 ;;
        blst_fuzz_bytes) BLST_FUZZ_MAX_LEN=1024; BLST_FUZZ_TIMEOUT=2 ;;
        blst_fuzz_fields) BLST_FUZZ_MAX_LEN=256; BLST_FUZZ_TIMEOUT=2 ;;
        blst_fuzz_groups) BLST_FUZZ_MAX_LEN=64; BLST_FUZZ_TIMEOUT=2 ;;
        blst_fuzz_protocol) BLST_FUZZ_MAX_LEN=2048; BLST_FUZZ_TIMEOUT=10 ;;
        blst_fuzz_pairing_state) BLST_FUZZ_MAX_LEN=1024; BLST_FUZZ_TIMEOUT=10 ;;
        blst_fuzz_uniq) BLST_FUZZ_MAX_LEN=4097; BLST_FUZZ_TIMEOUT=10 ;;
        blst_fuzz_fft_middle) BLST_FUZZ_MAX_LEN=128; BLST_FUZZ_TIMEOUT=10 ;;
        blst_fuzz_msm) BLST_FUZZ_MAX_LEN=64; BLST_FUZZ_TIMEOUT=10 ;;
        blst_fuzz_avx) BLST_FUZZ_MAX_LEN=192; BLST_FUZZ_TIMEOUT=10 ;;
        *)
            echo "unknown fuzz target: $1" >&2
            return 2
            ;;
    esac
}

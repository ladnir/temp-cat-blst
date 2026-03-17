# VROOM Backend Reset Plan Status

This document records the concrete execution baseline and current status for
`docs/vroom-backend-reset-plan.md`.

## Phase 0 Baseline

- upstream VROOM source commit: `31e98bbf7f32e1345299d4276dd6a5bc44dc0b69`
- local bring-up build: `x64-Debug` preset with `BLST_AVX512_PORTABLE=ON`
- parity compiler baseline: Clang 21 with `-O3 -march=native`
- correctness reference for parity tests: scalar `blst`
- intended backend split:
  - packed AVX owns `Fp`
  - packed AVX owns `Fp2`
  - VROOM owns `Fp12` and higher extension-field work
- integration rule during reset:
  - scalar code is not modified as a backend target
  - inactive packed-AVX tower code remains out of the build
  - no API integration work is allowed before the replacement AVX/VROOM backend
    surface is fully present

## Phase 1 Status

- remove from active build:
  - packed AVX tower files above `Fp2`
  - current hand-ported VROOM integration
  - tests tied only to the removed backends
  - benchmarks tied only to the removed backends
- keep in-tree but inactive:
  - packed AVX tower source files above `Fp2`, with explicit warning banners
- local validation target:
  - configure and build the `x64-Debug` preset
  - run the remaining relevant unit tests
- local validation result:
  - `cmake --preset x64-Debug` succeeded
  - `cmake --build out/build/x64-Debug --config Debug` succeeded
  - `ctest --test-dir out/build/x64-Debug --output-on-failure -C Debug` passed
  - passing test count: `240/240`

## Remote Parity Baseline

- benchmark server access: via the configured PuTTY/`plink.exe` workflow
- target compiler task:
  - prefer an official or prebuilt Clang 21 install
  - fall back to a source build only if a suitable packaged binary is not available
- installed compiler result:
  - remote host used: `18.188.231.170`
  - installed compiler version: Clang `21.1.7`
  - install root: `/opt/llvm-21`
  - command path: `/usr/local/bin/clang` and `/usr/local/bin/clang++`
  - smoke compile result: pass
- benchmark policy:
  - parity is measured on the remote AVX server
  - local Windows runs are for build and correctness only

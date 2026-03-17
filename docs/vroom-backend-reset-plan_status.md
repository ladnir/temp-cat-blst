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

## Phase 2 Status

- created landing-zone root:
  - `blst/vroom/`
- created copied-core directories:
  - `blst/vroom/cpu/vector/`
  - `blst/vroom/cpu/reduction/`
  - `blst/vroom/cpu/precompute/`
  - `blst/vroom/src/`
- created `blst`-owned directories:
  - `blst/vroom/adapters/`
  - `blst/vroom/tests/`
  - `blst/vroom/bench/`
- added a stub-free CMake hook:
  - `blst/CMakeLists.txt` now includes `blst/vroom/CMakeLists.txt`
  - `blst/vroom/CMakeLists.txt` defines no targets yet
- recorded the copied-core versus adapter boundary in:
  - `blst/vroom/README.md`

## Phase 3 Status

- copied literal upstream core from commit:
  - `31e98bbf7f32e1345299d4276dd6a5bc44dc0b69`
- import mechanism:
  - `git archive` from the pinned commit
  - copied from commit contents, not the sibling working tree
- imported directories:
  - `blst/vroom/cpu/precompute/`
  - `blst/vroom/cpu/reduction/`
  - `blst/vroom/cpu/vector/`
  - `blst/vroom/src/`
  - `blst/vroom/test_data/`
- imported file count:
  - upstream `cpu/` + `src/`: `39`
  - in-tree copied `blst/vroom/cpu/` + `blst/vroom/src/`: `39`
- current build state after import:
  - no VROOM targets are enabled yet
  - `cmake --build out/build/x64-Debug --config Debug` still succeeds unchanged
- current gaps before parity work:
  - no include-path rebinding has been done yet
  - no symbol collision review has been done yet
  - no `Fp12` parity targets are wired yet

## Phase 4 Status

- imported the upstream vector fixture directory needed by copied base-field
  headers:
  - `blst/vroom/test_data/`
- enabled a guarded copied-VROOM dependency target:
  - `blst/vroom/CMakeLists.txt` now defines `blst_vroom_headers`
  - portable local builds define `USE_AVX512_FALLBACK`
  - Windows local builds bind against the installed `mpir` headers/libs under
    `../vcpkg/installed/x64-windows/`
  - non-Windows builds bind against discovered `gmpxx` and `gmp`
- wired the first copied-VROOM base-field targets:
  - `test/CMakeLists.txt` adds `blst/vroom/tests/Vroom_basefield_tests.cpp`
  - `benchmark/CMakeLists.txt` adds `blst/vroom/bench/Vroom_basefield_bench.cpp`
- added Windows-local portability glue needed for the copied base-field layer:
  - copied VROOM preprocess macros now map to MSVC inline/noinline forms
  - copied VROOM fallback/vector code now uses an in-tree 128-bit compatibility
    shim on MSVC
  - copied GMP wrapper now exposes explicit 64-bit constructors and `to_u64()`
  - Windows test and benchmark targets stage `mpir.dll` into their output
    directories for runtime discovery
- phase scope of the new coverage:
  - deterministic `from_bigint` RNS vector checks for copied VROOM `Fp`
  - deterministic `modmul` known-answer checks for copied VROOM `Fp`
  - randomized `batch_modmul<8>` parity checks against scalar `blst::Fp`
  - randomized `Fp2` multiply parity checks against scalar `blst::Fp2`
  - copied-shape base-field benchmarks for `modmul` and `batch_modmul<8>`
- later `Fp2` gap-closure added on the same copied layer:
  - randomized `Fp2` square parity checks against scalar `blst::Fp2`
  - deterministic `Fp2` edge-value and identity checks for:
    - zero
    - one
    - negative-one
    - maximal reduced inputs
  - copied-shape direct `Fp2` benchmarks for:
    - `BM_FP2_Mul`
    - `BM_FP2_Sqr`
- copied-core build-note for Linux parity:
  - the extra `BigInt(long long)` and `BigInt(unsigned long long)` shim
    constructors in `blst/vroom/cpu/precompute/gmp_wrapper.hpp` are now gated
    to MSVC only so Linux `gmpxx` does not see ambiguous overloads
- local validation result after Phase 4 wiring:
  - `cmake --preset x64-Debug` succeeded
  - `cmake --build out/build/x64-Debug --config Debug` succeeded
  - `ctest --test-dir out/build/x64-Debug --output-on-failure -C Debug` passed
  - passing test count: `244/244`
  - benchmark registration smoke passed for:
    - `BM_ModMul`
    - `BM_BatchModMul_8`

## Phase 5 Status

- added direct copied-VROOM `Fp12` arithmetic coverage:
  - `blst/vroom/tests/Vroom_fp12_tests.cpp`
- functionality now covered in-tree:
  - randomized parity for `FP12Ring::fp12_flat_mul`
  - randomized parity for `FP12Ring::fp12_flat_sqr`
  - deterministic fixed-vector parity for `FP12Ring::fp12_flat_mul`
  - deterministic boundary round-trips through the current `Fp12` helper layer
  - scalar-reference coefficient mapping between copied VROOM flat `Fp12` and
    scalar `blst::Fp12`
- added copied-VROOM `Fp12` benchmark hooks:
  - `blst/vroom/bench/Vroom_fp12_bench.cpp`
  - registered `BM_FP12_Mul`
  - registered `BM_FP12_Sqr`
- aligned local `Fp12` benchmark input seeding to upstream benchmark shape:
  - multiply seeds: `42` and `123`
  - square seed: `42`
- local validation result after Phase 5 wiring:
  - `cmake --build out/build/x64-Debug --config Debug` succeeded
  - `ctest --test-dir out/build/x64-Debug --output-on-failure -C Debug` passed
  - passing test count: `250/250`
  - benchmark smoke passed for:
    - `BM_FP2_Mul`
    - `BM_FP2_Sqr`
    - `BM_FP12_Mul`
    - `BM_FP12_Sqr`

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

## Remote Parity Checkpoint 1

- machine:
  - remote AVX host `18.188.231.170`
  - single-core pinned runs via `taskset -c 0`
- compiler:
  - `/opt/llvm-21/bin/clang++`
  - version `21.1.7`
- in-tree build:
  - source root: `/home/ec2-user/parity_blst_current`
  - configure shape:
    - `Release`
    - `-O3 -march=native`
    - `BLST_PORTABLE=OFF`
    - `BLST_AVX512=ON`
    - `BLST_DEFAULT_AVX=ON`
    - `BLST_AVX512_PORTABLE=OFF`
- upstream build:
  - source root: `/home/ec2-user/parity_VROOM_upstream`
  - benchmark binary: `/home/ec2-user/parity_VROOM_upstream/src/bench_pairing`
  - benchmark headers/libs reused from the in-tree Clang 21 benchmark build
- benchmark command shape:
  - `taskset -c 0 <binary> --benchmark_filter='^BM_ModMul$' --benchmark_repetitions=5 --benchmark_report_aggregates_only=true --benchmark_min_time=0.05s`
  - `taskset -c 0 <binary> --benchmark_filter='^BM_FP12_Mul$' --benchmark_repetitions=5 --benchmark_report_aggregates_only=true --benchmark_min_time=0.05s`
- results:
  - upstream `BM_ModMul` median: `36.0 ns`
  - in-tree `BM_ModMul` median: `35.9 ns`
  - median drift: about `-0.3%`
  - parity result: pass
  - upstream `BM_FP12_Mul` median: `432 ns`
  - in-tree `BM_FP12_Mul` median: `432 ns`
  - median drift: `0.0%`
  - parity result: pass
- gate interpretation:
  - both required checkpoint medians are within the `5%` default pass rule
  - no persistent directional drift is visible in this checkpoint

## Phase 6 Status

- added a scalar `Fp12` multiplication baseline benchmark:
  - `benchmark/Vroom_fp12_scalar_baseline_bench.cpp`
  - benchmark name: `BM_FP12_ScalarMulBaseline`
- aligned direct VROOM and scalar `Fp12` benchmark inputs through one shared
  seeded coefficient helper:
  - `blst/vroom/tests/vroom_basefield_test_utils.hpp`
  - `blst/vroom/bench/Vroom_fp12_bench.cpp`
- added a repeatable remote parity runner:
  - `scripts/run_vroom_remote_parity.sh`
- parity runner behavior:
  - packages the current `blst` worktree snapshot over `HEAD`
  - packages upstream VROOM from pinned commit
    `31e98bbf7f32e1345299d4276dd6a5bc44dc0b69`
  - uploads both snapshots to the AVX-512 server
  - rebuilds both sides with Clang 21 and `-O3 -march=native`
  - runs:
    - upstream `BM_ModMul`
    - upstream `BM_FP12_Mul`
    - in-tree `BM_ModMul`
    - in-tree `BM_FP12_Mul`
    - in-tree scalar `BM_FP12_ScalarMulBaseline`
  - writes JSON reports plus a text summary under `out/vroom-parity/<run-id>/`
- latest scripted remote parity result:
  - run id: `20260317T193705Z`
  - report directory: `out/vroom-parity/20260317T193705Z/`
  - upstream `BM_ModMul` median: `35.938 ns`
  - in-tree `BM_ModMul` median: `35.969 ns`
  - median drift: `+0.085%`
  - upstream `BM_FP12_Mul` median: `431.805 ns`
  - in-tree `BM_FP12_Mul` median: `432.173 ns`
  - median drift: `+0.085%`
  - scalar `BM_FP12_ScalarMulBaseline` median: `1836.087 ns`
  - parity result: pass for both required gates
- local validation after Phase 6 wiring:
  - `cmake --build out/build/x64-Debug --config Debug --target blst_benchmarks`
    succeeded
  - `ctest --test-dir out/build/x64-Debug --output-on-failure -C Debug` passed
  - passing test count: `250/250`

## Phase 7 Status

- added the next isolated copied-VROOM `Fp12` operation:
  - cyclotomic square via `FP12Ring::fp12_flat_cyclotomic_sqr`
- added deterministic correctness coverage:
  - `blst/vroom/tests/Vroom_fp12_tests.cpp`
  - test name: `VroomFp12.CyclotomicSqrMatchesScalarBlstFp12`
  - input shape:
    - seeded scalar `Fp12`
    - scalar `final_exp` to enter the cyclotomic subgroup
    - convert back to copied VROOM coefficients
- added direct benchmark parity coverage:
  - `blst/vroom/bench/Vroom_fp12_cyclotomic_bench.cpp`
  - benchmark name: `BM_FP12_CyclotomicSqr`
- added scalar reference benchmark:
  - `benchmark/Vroom_fp12_scalar_cyclotomic_baseline_bench.cpp`
  - benchmark name: `BM_FP12_ScalarCyclotomicSqrBaseline`
- benchmark-layout note:
  - adding `BM_FP12_CyclotomicSqr` into the existing `Vroom_fp12_bench.cpp`
    translation unit caused a persistent `BM_FP12_Mul` parity miss just above the
    `5%` gate
  - the cyclotomic benchmark and scalar cyclotomic baseline were moved into their
    own translation units to restore the previously established `BM_FP12_Mul`
    benchmark shape
  - after that split, `BM_FP12_Mul` parity returned to pass
- latest scripted remote parity result:
  - run id: `20260317T200410Z`
  - report directory: `out/vroom-parity/20260317T200410Z/`
  - upstream `BM_ModMul` median: `35.939 ns`
  - in-tree `BM_ModMul` median: `35.986 ns`
  - median drift: `+0.130%`
  - upstream `BM_FP12_Mul` median: `431.938 ns`
  - in-tree `BM_FP12_Mul` median: `431.931 ns`
  - median drift: `-0.002%`
  - upstream `BM_FP12_CyclotomicSqr` median: `349.855 ns`
  - in-tree `BM_FP12_CyclotomicSqr` median: `350.132 ns`
  - median drift: `+0.079%`
  - scalar `BM_FP12_ScalarMulBaseline` median: `1834.359 ns`
  - scalar `BM_FP12_ScalarCyclotomicSqrBaseline` median: `825.861 ns`
  - parity result: pass for all tracked gates
- local validation after Phase 7 wiring:
  - `cmake --build out/build/x64-Debug --config Debug` succeeded
  - `ctest --test-dir out/build/x64-Debug --output-on-failure -C Debug` passed
  - passing test count: `251/251`

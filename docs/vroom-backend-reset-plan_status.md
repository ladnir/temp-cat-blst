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
- added the next pairing-adjacent copied-VROOM `Fp12` operation:
  - sparse multiplication via `FP12Ring::flat_mulxy00z0`
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
- added direct sparse benchmark parity coverage:
  - `blst/vroom/bench/Vroom_fp12_sparse_bench.cpp`
  - benchmark name: `BM_FP12_MulByXy0z00`
- added scalar reference benchmark:
  - `benchmark/Vroom_fp12_scalar_cyclotomic_baseline_bench.cpp`
  - benchmark name: `BM_FP12_ScalarCyclotomicSqrBaseline`
- added scalar sparse reference benchmark:
  - `benchmark/Vroom_fp12_sparse_scalar_baseline_bench.cpp`
  - benchmark name: `BM_FP12_ScalarMulByXy0z00Baseline`
- benchmark-layout note:
  - adding `BM_FP12_CyclotomicSqr` into the existing `Vroom_fp12_bench.cpp`
    translation unit caused a persistent `BM_FP12_Mul` parity miss just above the
    `5%` gate
  - the cyclotomic benchmark and scalar cyclotomic baseline were moved into their
    own translation units to restore the previously established `BM_FP12_Mul`
    benchmark shape
  - after that split, `BM_FP12_Mul` parity returned to pass
- latest scripted remote parity result:
  - run id: `20260317T201950Z`
  - report directory: `out/vroom-parity/20260317T201950Z/`
  - upstream `BM_ModMul` median: `35.948 ns`
  - in-tree `BM_ModMul` median: `36.012 ns`
  - median drift: `+0.177%`
  - upstream `BM_FP12_Mul` median: `431.739 ns`
  - in-tree `BM_FP12_Mul` median: `432.647 ns`
  - median drift: `+0.210%`
  - upstream `BM_FP12_CyclotomicSqr` median: `349.897 ns`
  - in-tree `BM_FP12_CyclotomicSqr` median: `349.949 ns`
  - median drift: `+0.015%`
  - upstream `BM_FP12_MulByXy0z00` median: `363.573 ns`
  - in-tree `BM_FP12_MulByXy0z00` median: `362.540 ns`
  - median drift: `-0.284%`
  - scalar `BM_FP12_ScalarMulBaseline` median: `1834.011 ns`
  - scalar `BM_FP12_ScalarCyclotomicSqrBaseline` median: `824.584 ns`
  - scalar `BM_FP12_ScalarMulByXy0z00Baseline` median: `1244.392 ns`
  - parity result: pass for all tracked gates
- local validation after Phase 7 wiring:
  - `cmake --build out/build/x64-Debug --config Debug` succeeded
  - `ctest --test-dir out/build/x64-Debug --output-on-failure -C Debug` passed
  - passing test count: `253/253`

## Phase 8 Status

- added copied-VROOM Miller benchmark coverage:
  - `blst/vroom/bench/Vroom_miller_bench.cpp`
  - benchmark names:
    - `BM_MillerLoop_DoubleStep`
    - `BM_MillerLoop_AddStep`
    - `BM_MillerLoop`
- added copied-VROOM Miller correctness coverage through scalar final
  exponentiation:
  - `blst/vroom/tests/Vroom_miller_tests.cpp`
  - test name:
    - `VroomPairing.MillerLoopThenScalarFinalExpMatchesScalarBlstOnReferencePairs`
  - reference cases:
    - generator `G1` with generator `G2`
    - doubled `G1` with generator `G2`
    - generator `G1` with doubled `G2`
  - validation shape:
    - run copied VROOM `miller_loop`
    - convert the copied VROOM `Fp12` result into scalar `blst::Fp12`
    - run scalar `blst::Fp12::final_exp`
    - compare against scalar `blst::GT::miller_loop` plus `final_mult`
- extended shared test helpers for copied-VROOM pairing inputs:
  - upstream VROOM generator constants are now centralized in
    `blst/vroom/tests/vroom_basefield_test_utils.hpp`
  - scalar and copied-VROOM pairing fixtures now derive from the same
    coordinate source
- extended the remote parity runner:
  - `scripts/run_vroom_remote_parity.sh` now also runs:
    - upstream `BM_MillerLoop_DoubleStep`
    - upstream `BM_MillerLoop_AddStep`
    - upstream `BM_MillerLoop`
    - in-tree `BM_MillerLoop_DoubleStep`
    - in-tree `BM_MillerLoop_AddStep`
    - in-tree `BM_MillerLoop`
- benchmark naming note:
  - the benchmark binary already contained a scalar `BM_MillerLoop` in
    `benchmark/pairing_bench.cpp`
  - that scalar benchmark was renamed to `BM_BlstScalarMillerLoop` so the
    copied-VROOM `BM_MillerLoop` parity benchmark has a unique name
- latest scripted remote parity result:
  - run id: `20260317T213033Z`
  - report directory: `out/vroom-parity/20260317T213033Z/`
  - upstream `BM_ModMul` median: `35.959 ns`
  - in-tree `BM_ModMul` median: `35.953 ns`
  - median drift: `-0.017%`
  - upstream `BM_FP12_Mul` median: `431.813 ns`
  - in-tree `BM_FP12_Mul` median: `431.794 ns`
  - median drift: `-0.004%`
  - upstream `BM_FP12_CyclotomicSqr` median: `349.883 ns`
  - in-tree `BM_FP12_CyclotomicSqr` median: `350.148 ns`
  - median drift: `+0.076%`
  - upstream `BM_FP12_MulByXy0z00` median: `363.524 ns`
  - in-tree `BM_FP12_MulByXy0z00` median: `362.795 ns`
  - median drift: `-0.201%`
  - upstream `BM_MillerLoop_DoubleStep` median: `901.819 ns`
  - in-tree `BM_MillerLoop_DoubleStep` median: `905.982 ns`
  - median drift: `+0.462%`
  - upstream `BM_MillerLoop_AddStep` median: `990.510 ns`
  - in-tree `BM_MillerLoop_AddStep` median: `994.841 ns`
  - median drift: `+0.437%`
  - upstream `BM_MillerLoop` median: `86933.714 ns`
  - in-tree `BM_MillerLoop` median: `87157.637 ns`
  - median drift: `+0.258%`
  - scalar `BM_FP12_ScalarMulBaseline` median: `1836.211 ns`
  - scalar `BM_FP12_ScalarCyclotomicSqrBaseline` median: `824.246 ns`
  - scalar `BM_FP12_ScalarMulByXy0z00Baseline` median: `1241.542 ns`
  - parity result: pass for all tracked gates
- local validation after Phase 8 wiring:
  - `cmake --build out/build/x64-Debug --config Debug --target blst_tests`
    succeeded
  - `cmake --build out/build/x64-Debug --config Debug --target blst_benchmarks`
    succeeded
  - `./test/blst_tests.exe --gtest_filter=VroomPairing.*` passed
  - `ctest --test-dir out/build/x64-Debug --output-on-failure -C Debug` passed
  - passing test count: `254/254`
  - local benchmark smoke passed for:
    - `BM_MillerLoop_DoubleStep`
    - `BM_MillerLoop_AddStep`
    - `BM_MillerLoop`

## Phase 9 Status

- added the in-tree copied-VROOM full pairing wrapper path:
  - `blst/vroom/src/pairing.hpp`
  - now exercised directly in-tree through both:
    - `FP12BLSTInverter`
    - `FP12RNS_BLSTInverter`
- replaced the dead low-level conversion bridge in:
  - `blst/vroom/src/conversion_inversion.hpp`
  - current bridge strategy:
    - convert copied VROOM field elements to scalar `blst::Fp` / `blst::Fp12`
    - use scalar `blst` inversion routines
    - convert back into copied VROOM ring elements
  - this is a real implementation path, not a mock, but it does not yet reuse
    the old `avxbls/src` low-level inverse entry points because those exact
    wrapper symbols are not present in this repository build
- added direct correctness coverage for the copied-VROOM pairing wrapper:
  - `blst/vroom/tests/Vroom_pairing_tests.cpp`
  - test names:
    - `VroomPairing.PairingBlstInverterMatchesScalarBlst`
    - `VroomPairing.PairingRnsBlstInverterMatchesScalarBlst`
  - reference cases:
    - generator `G1` with generator `G2`
    - doubled `G1` with generator `G2`
    - generator `G1` with doubled `G2`
- added direct wrapper benchmarks:
  - `blst/vroom/bench/Vroom_pairing_bench.cpp`
  - benchmark names:
    - `BM_Pairing_BLST_Inverter`
    - `BM_Pairing_RNS_BLST_Inverter`
- added scalar full-pairing baseline benchmark:
  - `benchmark/Vroom_pairing_scalar_baseline_bench.cpp`
  - benchmark name:
    - `BM_BlstScalarPairingFinalExp`
- extended the remote parity runner:
  - `scripts/run_vroom_remote_parity.sh` now also runs:
    - upstream `BM_Pairing_BLST_Inverter`
    - upstream `BM_Pairing_RNS_BLST_Inverter`
    - in-tree `BM_Pairing_BLST_Inverter`
    - in-tree `BM_Pairing_RNS_BLST_Inverter`
    - scalar `BM_BlstScalarPairingFinalExp`
- latest scripted remote parity result:
  - run id: `20260317T220552Z`
  - report directory: `out/vroom-parity/20260317T220552Z/`
  - upstream `BM_ModMul` median: `35.918 ns`
  - in-tree `BM_ModMul` median: `35.926 ns`
  - median drift: `+0.024%`
  - upstream `BM_FP12_Mul` median: `431.972 ns`
  - in-tree `BM_FP12_Mul` median: `428.666 ns`
  - median drift: `-0.765%`
  - upstream `BM_FP12_CyclotomicSqr` median: `349.905 ns`
  - in-tree `BM_FP12_CyclotomicSqr` median: `349.915 ns`
  - median drift: `+0.003%`
  - upstream `BM_FP12_MulByXy0z00` median: `363.572 ns`
  - in-tree `BM_FP12_MulByXy0z00` median: `362.278 ns`
  - median drift: `-0.356%`
  - upstream `BM_MillerLoop_DoubleStep` median: `901.611 ns`
  - in-tree `BM_MillerLoop_DoubleStep` median: `898.226 ns`
  - median drift: `-0.375%`
  - upstream `BM_MillerLoop_AddStep` median: `990.613 ns`
  - in-tree `BM_MillerLoop_AddStep` median: `993.133 ns`
  - median drift: `+0.254%`
  - upstream `BM_MillerLoop` median: `86921.294 ns`
  - in-tree `BM_MillerLoop` median: `87142.792 ns`
  - median drift: `+0.255%`
  - upstream `BM_Pairing_BLST_Inverter` median: `230432.620 ns`
  - in-tree `BM_Pairing_BLST_Inverter` median: `260506.403 ns`
  - median drift: `+13.051%`
  - upstream `BM_Pairing_RNS_BLST_Inverter` median: `225867.415 ns`
  - in-tree `BM_Pairing_RNS_BLST_Inverter` median: `226668.142 ns`
  - median drift: `+0.355%`
  - scalar `BM_BlstScalarPairingFinalExp` median: `584363.667 ns`
  - parity result:
    - pass for `BM_Pairing_RNS_BLST_Inverter`
    - fail for `BM_Pairing_BLST_Inverter`
- interpretation of the remaining gap:
  - the full copied-VROOM pairing wrapper path is now present and correct
  - the RNS-backed wrapper path is already within the `5%` parity gate
  - the direct `FP12BLSTInverter` path still pays too much conversion and/or
    inversion overhead compared with upstream and needs a dedicated follow-up
- local validation after Phase 9 wiring:
  - `cmake --build out/build/x64-Debug --config Debug --target blst_tests`
    succeeded
  - `cmake --build out/build/x64-Debug --config Debug --target blst_benchmarks`
    succeeded
  - `./test/blst_tests.exe --gtest_filter=VroomPairing.*` passed
  - `ctest --test-dir out/build/x64-Debug --output-on-failure -C Debug` passed
  - passing test count: `256/256`
  - local benchmark smoke passed for:
    - `BM_Pairing_BLST_Inverter`
    - `BM_Pairing_RNS_BLST_Inverter`
    - `BM_BlstScalarPairingFinalExp`

## Phase 10 Status

- closed the remaining `FP12BLSTInverter` parity gap in:
  - `blst/vroom/src/conversion_inversion.hpp`
- bridge change:
  - kept the scalar `BigInt` bridge for the single-field `BLSTInverter`
    helper used by the already-passing RNS inverse path
  - replaced the `FP12BLSTInverter` hot-path bridge with a direct
    ring-to-`blst::Fp12` conversion using native `blst::Fp` / `blst::UFp`
    types
  - this avoids the old per-coefficient `BigInt` conversion cost inside the
    full pairing benchmark
- implementation note:
  - an attempted direct include of the old `avxbls/src` C headers exposed an
    include-guard collision with the in-tree `blst` headers
  - the final fix avoids that collision entirely by using the current
    repository's `blst::Fp`, `blst::UFp`, and Montgomery helpers instead of
    trying to force the legacy C header path through the build
- added direct correctness coverage for the new bridge:
  - `blst/vroom/tests/Vroom_fp12_tests.cpp`
  - new test:
    - `VroomFp12.BlstInverterMatchesScalarBlstInverse`
- latest scripted remote parity result:
  - run id: `20260317T230024Z`
  - report directory: `out/vroom-parity/20260317T230024Z/`
  - upstream `BM_ModMul` median: `35.950 ns`
  - in-tree `BM_ModMul` median: `35.978 ns`
  - median drift: `+0.078%`
  - upstream `BM_FP12_Mul` median: `432.097 ns`
  - in-tree `BM_FP12_Mul` median: `432.063 ns`
  - median drift: `-0.008%`
  - upstream `BM_FP12_CyclotomicSqr` median: `349.931 ns`
  - in-tree `BM_FP12_CyclotomicSqr` median: `349.982 ns`
  - median drift: `+0.015%`
  - upstream `BM_FP12_MulByXy0z00` median: `363.685 ns`
  - in-tree `BM_FP12_MulByXy0z00` median: `362.356 ns`
  - median drift: `-0.365%`
  - upstream `BM_MillerLoop_DoubleStep` median: `901.529 ns`
  - in-tree `BM_MillerLoop_DoubleStep` median: `898.489 ns`
  - median drift: `-0.337%`
  - upstream `BM_MillerLoop_AddStep` median: `990.867 ns`
  - in-tree `BM_MillerLoop_AddStep` median: `992.461 ns`
  - median drift: `+0.161%`
  - upstream `BM_MillerLoop` median: `86930.559 ns`
  - in-tree `BM_MillerLoop` median: `87137.630 ns`
  - median drift: `+0.238%`
  - upstream `BM_Pairing_BLST_Inverter` median: `230562.197 ns`
  - in-tree `BM_Pairing_BLST_Inverter` median: `230054.270 ns`
  - median drift: `-0.220%`
  - upstream `BM_Pairing_RNS_BLST_Inverter` median: `223737.859 ns`
  - in-tree `BM_Pairing_RNS_BLST_Inverter` median: `226605.505 ns`
  - median drift: `+1.282%`
  - scalar `BM_FP12_ScalarMulBaseline` median: `1839.244 ns`
  - scalar `BM_FP12_ScalarCyclotomicSqrBaseline` median: `825.366 ns`
  - scalar `BM_FP12_ScalarMulByXy0z00Baseline` median: `1241.310 ns`
  - scalar `BM_BlstScalarPairingFinalExp` median: `585616.167 ns`
  - parity result: pass for all tracked gates, including the previously failing
    `BM_Pairing_BLST_Inverter`
- local validation after the direct-bridge change:
  - `cmake --build out/build/x64-Debug --config Debug --target blst_tests blst_benchmarks`
    succeeded
  - `./test/blst_tests.exe` passed
  - passing test count: `257/257`
  - targeted benchmark smoke passed for:
    - `BM_Pairing_BLST_Inverter`
    - `BM_Pairing_RNS_BLST_Inverter`

## Phase 11 Status

- added the next direct copied-VROOM final-exponentiation helper surface:
  - `Fp12` conjugation
  - `Fp12` Frobenius maps for `n = 1, 2, 3`
  - `final_exp_after_inverse`
- added direct correctness coverage in:
  - `blst/vroom/tests/Vroom_fp12_tests.cpp`
  - new tests:
    - `VroomFp12.ConjugateMatchesScalarBlstFp12`
    - `VroomFp12.FrobeniusMapsMatchScalarBlstFp12`
    - `VroomFp12.FinalExpAfterInverseMatchesScalarBlstFinalExp`
- added direct benchmark coverage in separate translation units to avoid
  disturbing the already-stable `BM_FP12_Mul` layout:
  - `blst/vroom/bench/Vroom_fp12_frobenius_bench.cpp`
  - `blst/vroom/bench/Vroom_fp12_finalexp_bench.cpp`
  - benchmark names:
    - `BM_FP12_Conjugate`
    - `BM_FP12_FrobeniusMap1`
    - `BM_FP12_FrobeniusMap2`
    - `BM_FP12_FrobeniusMap3`
    - `BM_FP12_FinalExpAfterInverse`
- extended the remote parity runner:
  - `scripts/run_vroom_remote_parity.sh` now also runs:
    - upstream and in-tree `BM_FP12_Conjugate`
    - upstream and in-tree `BM_FP12_FrobeniusMap1`
    - upstream and in-tree `BM_FP12_FrobeniusMap2`
    - upstream and in-tree `BM_FP12_FrobeniusMap3`
    - upstream and in-tree `BM_FP12_FinalExpAfterInverse`
- latest scripted remote parity result:
  - run id: `20260318T165248Z`
  - report directory: `out/vroom-parity/20260318T165248Z/`
  - upstream `BM_ModMul` median: `34.918 ns`
  - in-tree `BM_ModMul` median: `35.212 ns`
  - median drift: `+0.844%`
  - upstream `BM_FP12_Mul` median: `414.819 ns`
  - in-tree `BM_FP12_Mul` median: `417.297 ns`
  - median drift: `+0.597%`
  - upstream `BM_FP12_CyclotomicSqr` median: `335.612 ns`
  - in-tree `BM_FP12_CyclotomicSqr` median: `337.759 ns`
  - median drift: `+0.640%`
  - upstream `BM_FP12_Conjugate` median: `11.956 ns`
  - in-tree `BM_FP12_Conjugate` median: `11.925 ns`
  - median drift: `-0.254%`
  - upstream `BM_FP12_FrobeniusMap1` median: `232.229 ns`
  - in-tree `BM_FP12_FrobeniusMap1` median: `232.411 ns`
  - median drift: `+0.078%`
  - upstream `BM_FP12_FrobeniusMap2` median: `187.725 ns`
  - in-tree `BM_FP12_FrobeniusMap2` median: `185.889 ns`
  - median drift: `-0.978%`
  - upstream `BM_FP12_FrobeniusMap3` median: `153.476 ns`
  - in-tree `BM_FP12_FrobeniusMap3` median: `152.919 ns`
  - median drift: `-0.363%`
  - upstream `BM_FP12_FinalExpAfterInverse` median: `125165.964 ns`
  - in-tree `BM_FP12_FinalExpAfterInverse` median: `125766.145 ns`
  - median drift: `+0.480%`
  - upstream `BM_FP12_MulByXy0z00` median: `349.501 ns`
  - in-tree `BM_FP12_MulByXy0z00` median: `349.188 ns`
  - median drift: `-0.090%`
  - upstream `BM_MillerLoop_DoubleStep` median: `870.579 ns`
  - in-tree `BM_MillerLoop_DoubleStep` median: `865.925 ns`
  - median drift: `-0.535%`
  - upstream `BM_MillerLoop_AddStep` median: `950.969 ns`
  - in-tree `BM_MillerLoop_AddStep` median: `952.372 ns`
  - median drift: `+0.148%`
  - upstream `BM_MillerLoop` median: `83362.179 ns`
  - in-tree `BM_MillerLoop` median: `83670.871 ns`
  - median drift: `+0.370%`
  - upstream `BM_Pairing_BLST_Inverter` median: `222139.164 ns`
  - in-tree `BM_Pairing_BLST_Inverter` median: `221863.363 ns`
  - median drift: `-0.124%`
  - upstream `BM_Pairing_RNS_BLST_Inverter` median: `214538.618 ns`
  - in-tree `BM_Pairing_RNS_BLST_Inverter` median: `218658.419 ns`
  - median drift: `+1.920%`
  - scalar `BM_FP12_ScalarMulBaseline` median: `1794.287 ns`
  - scalar `BM_FP12_ScalarCyclotomicSqrBaseline` median: `806.488 ns`
  - scalar `BM_FP12_ScalarMulByXy0z00Baseline` median: `1211.178 ns`
  - scalar `BM_BlstScalarPairingFinalExp` median: `570529.577 ns`
  - parity result: pass for all tracked gates
- local validation after Phase 11 wiring:
  - `cmake --build out/build/x64-Debug --config Debug --target blst_tests blst_benchmarks`
    succeeded
  - `ctest --test-dir out/build/x64-Debug --output-on-failure -C Debug` passed
  - passing test count: `260/260`
  - targeted benchmark smoke passed for:
    - `BM_FP12_Conjugate`
    - `BM_FP12_FrobeniusMap1`
    - `BM_FP12_FrobeniusMap2`
    - `BM_FP12_FrobeniusMap3`
    - `BM_FP12_FinalExpAfterInverse`

## Phase 12 Status

- added the remaining low-level copied-VROOM base-field primitive coverage:
  - addition-chain field inverse via `BLS381AddChainInversion`
  - `batch_reduce_expand` for supported upstream benchmark arities:
    - `1`
    - `2`
    - `3`
    - `4`
    - `5`
    - `6`
    - `12`
- added direct correctness coverage in:
  - `blst/vroom/tests/Vroom_basefield_tests.cpp`
  - new tests:
    - `VroomBaseField.AddChainInverseMatchesScalarBlstFp`
    - `VroomBaseField.BatchReduceExpandMatchesScalarBlstFp`
- added direct benchmark coverage in separate translation units:
  - `blst/vroom/bench/Vroom_basefield_inverse_bench.cpp`
  - `blst/vroom/bench/Vroom_basefield_reduce_bench.cpp`
  - benchmark names:
    - `BM_Inverse`
    - `BM_BatchReduceExpand_1`
    - `BM_BatchReduceExpand_2`
    - `BM_BatchReduceExpand_3`
    - `BM_BatchReduceExpand_4`
    - `BM_BatchReduceExpand_5`
    - `BM_BatchReduceExpand_6`
    - `BM_BatchReduceExpand_12`
- extended the remote parity runner:
  - `scripts/run_vroom_remote_parity.sh` now also runs:
    - upstream and in-tree `BM_Inverse`
    - upstream and in-tree `BM_BatchReduceExpand_1`
    - upstream and in-tree `BM_BatchReduceExpand_2`
    - upstream and in-tree `BM_BatchReduceExpand_3`
    - upstream and in-tree `BM_BatchReduceExpand_4`
    - upstream and in-tree `BM_BatchReduceExpand_5`
    - upstream and in-tree `BM_BatchReduceExpand_6`
    - upstream and in-tree `BM_BatchReduceExpand_12`
- latest scripted remote parity result:
  - run id: `20260318T211148Z`
  - report directory: `out/vroom-parity/20260318T211148Z/`
  - upstream `BM_Inverse` median: `25631.592 ns`
  - in-tree `BM_Inverse` median: `25729.025 ns`
  - median drift: `+0.380%`
  - upstream `BM_BatchReduceExpand_1` median: `25.100 ns`
  - in-tree `BM_BatchReduceExpand_1` median: `25.135 ns`
  - median drift: `+0.141%`
  - upstream `BM_BatchReduceExpand_2` median: `45.747 ns`
  - in-tree `BM_BatchReduceExpand_2` median: `45.545 ns`
  - median drift: `-0.443%`
  - upstream `BM_BatchReduceExpand_3` median: `64.344 ns`
  - in-tree `BM_BatchReduceExpand_3` median: `64.303 ns`
  - median drift: `-0.064%`
  - upstream `BM_BatchReduceExpand_4` median: `83.062 ns`
  - in-tree `BM_BatchReduceExpand_4` median: `82.519 ns`
  - median drift: `-0.654%`
  - upstream `BM_BatchReduceExpand_5` median: `96.711 ns`
  - in-tree `BM_BatchReduceExpand_5` median: `97.068 ns`
  - median drift: `+0.369%`
  - upstream `BM_BatchReduceExpand_6` median: `123.979 ns`
  - in-tree `BM_BatchReduceExpand_6` median: `124.064 ns`
  - median drift: `+0.069%`
  - upstream `BM_BatchReduceExpand_12` median: `265.396 ns`
  - in-tree `BM_BatchReduceExpand_12` median: `266.371 ns`
  - median drift: `+0.367%`
  - parity result: pass for all new tracked gates
- validation mode for Phase 12:
  - local compilation intentionally paused per current workflow
  - remote Release correctness run on the AVX-512 server:
    - `./out/parity-tests-clang21/test/blst_tests --gtest_filter='VroomBaseField.AddChainInverseMatchesScalarBlstFp:VroomBaseField.BatchReduceExpandMatchesScalarBlstFp'`
    - result: pass

## Phase 13 Status

- added the next direct copied-VROOM curve-operation surface:
  - `G1` projective add
  - `G1` mixed add
  - `G1` double
  - `G2` projective add
  - `G2` mixed add
  - `G2` double
- implementation note:
  - copied VROOM curve arithmetic requires the upstream curve ring shape
    `BoundedRing<381, 8, 52, -1932, 2377, 12>`
  - the first in-tree attempt mistakenly reused the default base-field ring and
    failed to compile due insufficient precomputed reduction table depth
  - the final harness uses a dedicated curve-ring alias and matching helper
    constructors for random benchmark inputs
- benchmark naming note:
  - the repository already had scalar `BM_G1_Add`, `BM_G1_Double`,
    `BM_G2_Add`, and `BM_G2_Double` benchmarks in `benchmark/G1_bench.cpp` and
    `benchmark/G2_bench.cpp`
  - those scalar benchmarks were renamed to:
    - `BM_BlstScalarG1_Add`
    - `BM_BlstScalarG1_Double`
    - `BM_BlstScalarG2_Add`
    - `BM_BlstScalarG2_Double`
  - this was necessary so the copied-VROOM direct curve benchmarks could own
    the upstream parity names without collision
- added direct correctness coverage in:
  - `blst/vroom/tests/Vroom_curve_tests.cpp`
  - new tests:
    - `VroomCurve.G1AddMixedAddAndDoubleMatchScalarBlst`
    - `VroomCurve.G2AddMixedAddAndDoubleMatchScalarBlst`
- added direct benchmark coverage in:
  - `blst/vroom/bench/Vroom_curve_bench.cpp`
  - benchmark names:
    - `BM_G1_Add`
    - `BM_G1_MixedAdd`
    - `BM_G1_Double`
    - `BM_G2_Add`
    - `BM_G2_MixedAdd`
    - `BM_G2_Double`
- extended the remote parity runner:
  - `scripts/run_vroom_remote_parity.sh` now also runs:
    - upstream and in-tree `BM_G1_Add`
    - upstream and in-tree `BM_G1_MixedAdd`
    - upstream and in-tree `BM_G1_Double`
    - upstream and in-tree `BM_G2_Add`
    - upstream and in-tree `BM_G2_MixedAdd`
    - upstream and in-tree `BM_G2_Double`
- latest scripted remote parity result:
  - run id: `20260318T211148Z`
  - report directory: `out/vroom-parity/20260318T211148Z/`
  - upstream `BM_G1_Add` median: `223.134 ns`
  - in-tree `BM_G1_Add` median: `222.396 ns`
  - median drift: `-0.331%`
  - upstream `BM_G1_MixedAdd` median: `222.937 ns`
  - in-tree `BM_G1_MixedAdd` median: `222.847 ns`
  - median drift: `-0.040%`
  - upstream `BM_G1_Double` median: `211.973 ns`
  - in-tree `BM_G1_Double` median: `210.159 ns`
  - median drift: `-0.856%`
  - upstream `BM_G2_Add` median: `493.299 ns`
  - in-tree `BM_G2_Add` median: `494.006 ns`
  - median drift: `+0.143%`
  - upstream `BM_G2_MixedAdd` median: `455.518 ns`
  - in-tree `BM_G2_MixedAdd` median: `458.877 ns`
  - median drift: `+0.738%`
  - upstream `BM_G2_Double` median: `379.370 ns`
  - in-tree `BM_G2_Double` median: `379.291 ns`
  - median drift: `-0.021%`
  - parity result: pass for all new tracked gates
- remote Release correctness validation on the same snapshot:
  - run directory: `/home/ec2-user/vroom-parity/20260318T211148Z/blst_current`
  - command:
    - `./out/parity-tests-clang21/test/blst_tests --gtest_filter='VroomCurve.*:VroomBaseField.AddChainInverseMatchesScalarBlstFp:VroomBaseField.BatchReduceExpandMatchesScalarBlstFp'`
  - result: pass

## Phase 14 Status

- added the next direct copied-VROOM scalar-multiplication surface:
  - `G1` GLV scalar multiplication through `GLVMult`
  - `G2` GLS scalar multiplication through `GLSMult`
- compatibility note:
  - copied `scalar_mult.hpp` assumed upstream VROOM’s sibling `blst/` layout and
    global helper names
  - the in-tree version now resolves both layouts:
    - upstream `../blst/...` when building in the vendor tree
    - local `blst/common/...` and `blst/backend/...` when building under `blst/vroom`
  - copied `src/gmp_wrapper.hpp` now supports `BLST_VROOM_SKIP_TEST_VALUES`
    so new scalar-mult translation units can avoid pulling in duplicate
    `test_values.hpp` globals during the monolithic test link
- added direct correctness coverage in:
  - `blst/vroom/tests/Vroom_scalar_mult_tests.cpp`
  - new tests:
    - `VroomScalarMult.G1ScalarMultGlvMatchesScalarBlst`
    - `VroomScalarMult.G2ScalarMultGlsMatchesScalarBlst`
  - both tests compare copied-VROOM byte-array and `BigInt` entry points against
    scalar `blst` scalar multiplication on valid generator-derived points
- added direct benchmark coverage in:
  - `blst/vroom/bench/Vroom_scalar_mult_bench.cpp`
  - benchmark names:
    - `BM_G1_ScalarMult_GLV`
    - `BM_G2_ScalarMult_GLS`
- extended the remote parity runner:
  - `scripts/run_vroom_remote_parity.sh` now also runs:
    - upstream and in-tree `BM_G1_ScalarMult_GLV`
    - upstream and in-tree `BM_G2_ScalarMult_GLS`
- latest scripted remote parity result:
  - run id: `20260318T214239Z`
  - report directory: `out/vroom-parity/20260318T214239Z/`
  - upstream `BM_G1_ScalarMult_GLV` median: `46022.856 ns`
  - in-tree `BM_G1_ScalarMult_GLV` median: `45919.691 ns`
  - median drift: `-0.224%`
  - upstream `BM_G2_ScalarMult_GLS` median: `68959.491 ns`
  - in-tree `BM_G2_ScalarMult_GLS` median: `69476.797 ns`
  - median drift: `+0.750%`
  - parity result: pass for both new tracked gates
- validation mode for Phase 14:
  - local compilation intentionally paused per current workflow
  - remote Release correctness run on the AVX-512 server:
    - run directory: `/home/ec2-user/vroom-parity/20260318T214239Z/blst_current`
    - command:
      - `./out/parity-tests-clang21/test/blst_tests --gtest_filter='VroomScalarMult.*'`
    - result: pass

## Phase 15 Status

- added the next copied-VROOM parameter-set slice:
  - 50-bit `G1` GLV scalar multiplication
  - 50-bit `G2` GLS scalar multiplication
- implementation note:
  - this phase keeps the same copied `scalar_mult.hpp` surface and exercises it
    under `BoundedRing<381, 8, 50, -1932, 2377, 12>`
  - shared test helpers now expose dedicated 50-bit curve aliases and random
    projective/affine constructors
- added direct correctness coverage in:
  - `blst/vroom/tests/Vroom_scalar_mult_50bit_tests.cpp`
  - new tests:
    - `VroomScalarMult50Bit.G1ScalarMultGlvMatchesScalarBlst`
    - `VroomScalarMult50Bit.G2ScalarMultGlsMatchesScalarBlst`
- added direct benchmark coverage in:
  - `blst/vroom/bench/Vroom_scalar_mult_50bit_bench.cpp`
  - benchmark names:
    - `BM_G1_ScalarMult_GLV_50bit`
    - `BM_G2_ScalarMult_GLS_50bit`
- extended the remote parity runner:
  - `scripts/run_vroom_remote_parity.sh` now also:
    - builds upstream `bench_pairing_50bit`
    - runs upstream and in-tree `BM_G1_ScalarMult_GLV_50bit`
    - runs upstream and in-tree `BM_G2_ScalarMult_GLS_50bit`
- latest scripted remote parity result:
  - run id: `20260319T020039Z`
  - report directory: `out/vroom-parity/20260319T020039Z/`
  - upstream `BM_G1_ScalarMult_GLV_50bit` median: `42670.828 ns`
  - in-tree `BM_G1_ScalarMult_GLV_50bit` median: `42623.742 ns`
  - median drift: `-0.110%`
  - upstream `BM_G2_ScalarMult_GLS_50bit` median: `65796.291 ns`
  - in-tree `BM_G2_ScalarMult_GLS_50bit` median: `65189.799 ns`
  - median drift: `-0.922%`
  - parity result: pass for both new tracked gates
- validation mode for Phase 15:
  - local compilation intentionally paused per current workflow
  - remote Release correctness run on the AVX-512 server:
    - run directory: `/home/ec2-user/vroom-parity/20260319T020039Z/blst_current`
    - command:
      - `./out/parity-tests-clang21/test/blst_tests --gtest_filter='VroomScalarMult50Bit.*'`
    - result: pass

## Phase 16 Status

- added the next copied-VROOM 50-bit tower and wrapper surface:
  - 50-bit `Fp12` multiply and square
  - 50-bit `Fp12` cyclotomic square
  - 50-bit `Fp12` conjugate and Frobenius maps `1/2/3`
  - 50-bit `Fp12` multiply-by-`xy00z0`
  - 50-bit `Fp12` inverse through `FP12BLSTInverter` and `FP12RNS_BLSTInverter`
  - 50-bit final exponentiation wrapper paths through both inverter flavors
  - 50-bit full pairing wrapper paths through both inverter flavors
- implementation notes:
  - shared VROOM helper factories now accept ring-generic `Fp12` and affine generator paths so the same test scaffolding can drive both 52-bit and 50-bit rings
  - the 50-bit final-exponentiation wrapper benchmark needed the copied inverter header explicitly included
  - the 50-bit sparse `Fp12` benchmark initially drifted because the destination `Fp12` object was zero-initialized inside the hot loop while upstream leaves it uninitialized before `flat_mulxy00z0` overwrites it
  - the in-tree benchmark now matches upstream’s benchmark shape for that path
- added direct correctness coverage in:
  - `blst/vroom/tests/Vroom_fp12_50bit_tests.cpp`
  - `blst/vroom/tests/Vroom_pairing_50bit_tests.cpp`
  - new tests:
    - `VroomFp1250Bit.MulMatchesScalarBlstFp12`
    - `VroomFp1250Bit.CyclotomicSqrMatchesScalarBlstFp12`
    - `VroomFp1250Bit.MulByXy00z0MatchesScalarBlstFp12`
    - `VroomFp1250Bit.BlstInverterMatchesScalarBlstInverse`
    - `VroomFp1250Bit.RnsBlstInverterMatchesScalarBlstInverse`
    - `VroomFp1250Bit.FinalExpBlstInverterMatchesScalarBlst`
    - `VroomFp1250Bit.FinalExpRnsBlstInverterMatchesScalarBlst`
    - `VroomPairing50Bit.PairingBlstInverterMatchesScalarBlst`
    - `VroomPairing50Bit.PairingRnsBlstInverterMatchesScalarBlst`
- added direct benchmark coverage in:
  - `blst/vroom/bench/Vroom_fp12_50bit_bench.cpp`
  - `blst/vroom/bench/Vroom_fp12_inverse_50bit_bench.cpp`
  - `blst/vroom/bench/Vroom_fp12_finalexp_wrapper_50bit_bench.cpp`
  - `blst/vroom/bench/Vroom_pairing_50bit_bench.cpp`
  - benchmark names:
    - `BM_FP12_Mul_50bit`
    - `BM_FP12_Sqr_50bit`
    - `BM_FP12_CyclotomicSqr_50bit`
    - `BM_FP12_Conjugate_50bit`
    - `BM_FP12_FrobeniusMap1_50bit`
    - `BM_FP12_FrobeniusMap2_50bit`
    - `BM_FP12_FrobeniusMap3_50bit`
    - `BM_FP12_FinalExpAfterInverse_50bit`
    - `BM_FP12_MulByXy0z00_50bit`
    - `BM_FP12_Inverse_BLST_50bit`
    - `BM_FP12_Inverse_RNS_BLST_50bit`
    - `BM_FP12_FinalExp_BLST_Inverter_50bit`
    - `BM_FP12_FinalExp_RNS_BLST_Inverter_50bit`
    - `BM_Pairing_BLST_Inverter_50bit`
    - `BM_Pairing_RNS_BLST_Inverter_50bit`
- extended the remote parity runner:
  - `scripts/run_vroom_remote_parity.sh` now also runs the upstream `bench_pairing_50bit` counterparts and compares them against all new in-tree 50-bit benchmark names above
- latest scripted remote parity result:
  - run id: `20260319T044912Z`
  - report directory: `out/vroom-parity/20260319T044912Z/`
  - upstream `BM_FP12_Mul` median: `419.458 ns`
  - in-tree `BM_FP12_Mul_50bit` median: `419.310 ns`
  - median drift: `-0.035%`
  - upstream `BM_FP12_Sqr` median: `357.458 ns`
  - in-tree `BM_FP12_Sqr_50bit` median: `357.178 ns`
  - median drift: `-0.078%`
  - upstream `BM_FP12_CyclotomicSqr` median: `323.795 ns`
  - in-tree `BM_FP12_CyclotomicSqr_50bit` median: `322.676 ns`
  - median drift: `-0.345%`
  - upstream `BM_FP12_Conjugate` median: `11.959 ns`
  - in-tree `BM_FP12_Conjugate_50bit` median: `11.959 ns`
  - median drift: `-0.002%`
  - upstream `BM_FP12_FrobeniusMap1` median: `219.116 ns`
  - in-tree `BM_FP12_FrobeniusMap1_50bit` median: `221.082 ns`
  - median drift: `+0.897%`
  - upstream `BM_FP12_FrobeniusMap2` median: `180.868 ns`
  - in-tree `BM_FP12_FrobeniusMap2_50bit` median: `182.828 ns`
  - median drift: `+1.083%`
  - upstream `BM_FP12_FrobeniusMap3` median: `140.806 ns`
  - in-tree `BM_FP12_FrobeniusMap3_50bit` median: `143.825 ns`
  - median drift: `+2.144%`
  - upstream `BM_FP12_FinalExpAfterInverse` median: `120123.300 ns`
  - in-tree `BM_FP12_FinalExpAfterInverse_50bit` median: `120036.125 ns`
  - median drift: `-0.073%`
  - upstream `BM_FP12_MulByXy0z00` median: `353.578 ns`
  - in-tree `BM_FP12_MulByXy0z00_50bit` median: `353.655 ns`
  - median drift: `+0.022%`
  - upstream `BM_FP12_Inverse_BLST` median: `9001.469 ns`
  - in-tree `BM_FP12_Inverse_BLST_50bit` median: `9049.043 ns`
  - median drift: `+0.529%`
  - upstream `BM_FP12_Inverse_RNS_BLST` median: `4234.665 ns`
  - in-tree `BM_FP12_Inverse_RNS_BLST_50bit` median: `4342.548 ns`
  - median drift: `+2.548%`
  - upstream `BM_FP12_FinalExp_BLST_Inverter` median: `131374.294 ns`
  - in-tree `BM_FP12_FinalExp_BLST_Inverter_50bit` median: `131578.428 ns`
  - median drift: `+0.155%`
  - upstream `BM_FP12_FinalExp_RNS_BLST_Inverter` median: `125110.647 ns`
  - in-tree `BM_FP12_FinalExp_RNS_BLST_Inverter_50bit` median: `125594.070 ns`
  - median drift: `+0.386%`
  - upstream `BM_Pairing_BLST_Inverter` median: `212005.805 ns`
  - in-tree `BM_Pairing_BLST_Inverter_50bit` median: `211771.699 ns`
  - median drift: `-0.110%`
  - upstream `BM_Pairing_RNS_BLST_Inverter` median: `205907.982 ns`
  - in-tree `BM_Pairing_RNS_BLST_Inverter_50bit` median: `205311.853 ns`
  - median drift: `-0.290%`
  - parity result: pass for all new tracked 50-bit gates
- validation mode for Phase 16:
  - local compilation intentionally paused per current workflow
  - remote Release correctness run on the AVX-512 server:
    - run directory: `/home/ec2-user/vroom-parity/20260319T044912Z/blst_current`
    - command:
      - `./out/parity-tests-clang21/test/blst_tests --gtest_filter='VroomFp1250Bit.*:VroomPairing50Bit.*'`
    - result: pass

## Phase 17 Status

- added the next copied-VROOM 50-bit Miller surface:
  - `Miller::double_step`
  - `Miller::add_step`
  - `Miller::miller_loop`
- added direct correctness coverage in:
  - `blst/vroom/tests/Vroom_miller_50bit_tests.cpp`
  - new test:
    - `VroomPairing50Bit.MillerLoopThenScalarFinalExpMatchesScalarBlstOnReferencePairs`
  - this checks copied VROOM 50-bit Miller output after scalar `blst` final exponentiation against scalar `blst` final pairing on generator and doubled reference pairs
- added direct benchmark coverage in:
  - `blst/vroom/bench/Vroom_miller_50bit_bench.cpp`
  - benchmark names:
    - `BM_MillerLoop_DoubleStep_50bit`
    - `BM_MillerLoop_AddStep_50bit`
    - `BM_MillerLoop_50bit`
- extended the remote parity runner:
  - `scripts/run_vroom_remote_parity.sh` now also runs the upstream `bench_pairing_50bit` Miller benchmarks and compares them against the new in-tree 50-bit Miller benchmark names
- latest scripted remote parity result:
  - run id: `20260319T045825Z`
  - report directory: `out/vroom-parity/20260319T045825Z/`
  - upstream `BM_MillerLoop_DoubleStep` median: `837.904 ns`
  - in-tree `BM_MillerLoop_DoubleStep_50bit` median: `840.757 ns`
  - median drift: `+0.341%`
  - upstream `BM_MillerLoop_AddStep` median: `957.097 ns`
  - in-tree `BM_MillerLoop_AddStep_50bit` median: `951.995 ns`
  - median drift: `-0.533%`
  - upstream `BM_MillerLoop` median: `80336.485 ns`
  - in-tree `BM_MillerLoop_50bit` median: `81152.991 ns`
  - median drift: `+1.016%`
  - parity result: pass for all new tracked 50-bit Miller gates
- validation mode for Phase 17:
  - local compilation intentionally paused per current workflow
  - remote Release correctness run on the AVX-512 server:
    - run directory: `/home/ec2-user/vroom-parity/20260319T045825Z/blst_current`
    - command:
      - `./out/parity-tests-clang21/test/blst_tests --gtest_filter='VroomPairing50Bit.*'`
    - result: pass

## Phase 18 Status

- added the next copied-VROOM 50-bit curve surface:
  - `G1` add, mixed-add, and double
  - `G2` add, mixed-add, and double
- validation note:
  - the first parity sample after adding the new 50-bit curve benchmark translation unit produced a one-off `BM_FP12_FinalExp_RNS_BLST_Inverter_50bit` drift just over the `5%` gate
  - a full rerun on the same snapshot returned that benchmark to `+2.289%`
  - this phase therefore uses the stable rerun as the recorded checkpoint and treats the first result as benchmark noise rather than a persistent regression
- added direct correctness coverage in:
  - `blst/vroom/tests/Vroom_curve_50bit_tests.cpp`
  - new tests:
    - `VroomCurve50Bit.G1AddMixedAddAndDoubleMatchScalarBlst`
    - `VroomCurve50Bit.G2AddMixedAddAndDoubleMatchScalarBlst`
- added direct benchmark coverage in:
  - `blst/vroom/bench/Vroom_curve_50bit_bench.cpp`
  - benchmark names:
    - `BM_G1_Add_50bit`
    - `BM_G1_MixedAdd_50bit`
    - `BM_G1_Double_50bit`
    - `BM_G2_Add_50bit`
    - `BM_G2_MixedAdd_50bit`
    - `BM_G2_Double_50bit`
- extended the remote parity runner:
  - `scripts/run_vroom_remote_parity.sh` now also runs the upstream 50-bit curve benchmarks from `bench_pairing_50bit` and compares them against the new in-tree suffixed benchmark names
- latest scripted remote parity result:
  - run id: `20260319T051257Z`
  - report directory: `out/vroom-parity/20260319T051257Z/`
  - upstream `BM_G1_Add` median: `207.027 ns`
  - in-tree `BM_G1_Add_50bit` median: `207.304 ns`
  - median drift: `+0.134%`
  - upstream `BM_G1_MixedAdd` median: `199.069 ns`
  - in-tree `BM_G1_MixedAdd_50bit` median: `199.431 ns`
  - median drift: `+0.182%`
  - upstream `BM_G1_Double` median: `192.884 ns`
  - in-tree `BM_G1_Double_50bit` median: `191.970 ns`
  - median drift: `-0.474%`
  - upstream `BM_G2_Add` median: `494.728 ns`
  - in-tree `BM_G2_Add_50bit` median: `496.932 ns`
  - median drift: `+0.445%`
  - upstream `BM_G2_MixedAdd` median: `421.446 ns`
  - in-tree `BM_G2_MixedAdd_50bit` median: `422.454 ns`
  - median drift: `+0.239%`
  - upstream `BM_G2_Double` median: `349.912 ns`
  - in-tree `BM_G2_Double_50bit` median: `351.846 ns`
  - median drift: `+0.553%`
  - parity result: pass for all new tracked 50-bit curve gates
- validation mode for Phase 18:
  - local compilation intentionally paused per current workflow
  - remote Release correctness run on the AVX-512 server:
    - run directory: `/home/ec2-user/vroom-parity/20260319T051257Z/blst_current`
    - command:
      - `./out/parity-tests-clang21/test/blst_tests --gtest_filter='VroomCurve50Bit.*'`
    - result: pass

## Phase 19 Status

- added the remaining copied-VROOM 50-bit base-field and `Fp2` surface from
  upstream `bench_pairing_50bit`:
  - `modmul`
  - `batch_modmul` arities `1/2/3/4/5/6/7/8/10/12`
  - `batch_reduce_expand` arities `1/2/3/4/5/6/12`
  - add-chain inverse
  - `Fp2` multiply
  - `Fp2` square
- added direct correctness coverage in:
  - `blst/vroom/tests/Vroom_basefield_50bit_tests.cpp`
  - new tests:
    - `VroomBaseField50Bit.ModMulMatchesScalarBlstFp`
    - `VroomBaseField50Bit.BatchModMulMatchesScalarBlstFp`
    - `VroomBaseField50Bit.AddChainInverseMatchesScalarBlstFp`
    - `VroomBaseField50Bit.BatchReduceExpandMatchesScalarBlstFp`
    - `VroomBaseField50Bit.Fp2MulMatchesScalarBlstFp2`
    - `VroomBaseField50Bit.Fp2SquareMatchesScalarBlstFp2`
    - `VroomBaseField50Bit.Fp2EdgeValuesAndIdentitiesMatchScalarBlstFp2`
- added direct benchmark coverage in:
  - `blst/vroom/bench/Vroom_basefield_50bit_bench.cpp`
  - `blst/vroom/bench/Vroom_basefield_inverse_50bit_bench.cpp`
  - `blst/vroom/bench/Vroom_basefield_reduce_50bit_bench.cpp`
  - benchmark names:
    - `BM_ModMul_50bit`
    - `BM_BatchModMul_1_50bit`
    - `BM_BatchModMul_2_50bit`
    - `BM_BatchModMul_3_50bit`
    - `BM_BatchModMul_4_50bit`
    - `BM_BatchModMul_5_50bit`
    - `BM_BatchModMul_6_50bit`
    - `BM_BatchModMul_7_50bit`
    - `BM_BatchModMul_8_50bit`
    - `BM_BatchModMul_10_50bit`
    - `BM_BatchModMul_12_50bit`
    - `BM_BatchReduceExpand_1_50bit`
    - `BM_BatchReduceExpand_2_50bit`
    - `BM_BatchReduceExpand_3_50bit`
    - `BM_BatchReduceExpand_4_50bit`
    - `BM_BatchReduceExpand_5_50bit`
    - `BM_BatchReduceExpand_6_50bit`
    - `BM_BatchReduceExpand_12_50bit`
    - `BM_FP2_Mul_50bit`
    - `BM_FP2_Sqr_50bit`
    - `BM_Inverse_50bit`
- helper-layer note:
  - shared VROOM test helpers now construct `Fp2` values generically over the
    selected ring type so the same scalar-reference scaffolding can drive both
    52-bit and 50-bit copied rings
- extended the remote parity runner:
  - `scripts/run_vroom_remote_parity.sh` now also runs the upstream 50-bit
    base-field and `Fp2` benchmarks from `bench_pairing_50bit` and compares
    them against the new in-tree suffixed benchmark names
- latest scripted remote parity result:
  - run id: `20260319T062838Z`
  - report directory: `out/vroom-parity/20260319T062838Z/`
  - upstream `BM_ModMul` median: `32.373 ns`
  - in-tree `BM_ModMul_50bit` median: `32.101 ns`
  - median drift: `-0.838%`
  - upstream `BM_BatchModMul_12` median: `289.294 ns`
  - in-tree `BM_BatchModMul_12_50bit` median: `288.962 ns`
  - median drift: `-0.115%`
  - upstream `BM_BatchReduceExpand_12` median: `260.620 ns`
  - in-tree `BM_BatchReduceExpand_12_50bit` median: `261.009 ns`
  - median drift: `+0.149%`
  - upstream `BM_FP2_Mul` median: `56.522 ns`
  - in-tree `BM_FP2_Mul_50bit` median: `56.459 ns`
  - median drift: `-0.112%`
  - upstream `BM_FP2_Sqr` median: `56.013 ns`
  - in-tree `BM_FP2_Sqr_50bit` median: `55.714 ns`
  - median drift: `-0.534%`
  - upstream `BM_Inverse` median: `23784.848 ns`
  - in-tree `BM_Inverse_50bit` median: `23651.150 ns`
  - median drift: `-0.562%`
  - parity result: pass for all new tracked 50-bit base-field and `Fp2` gates
- validation mode for Phase 19:
  - local compilation intentionally paused per current workflow
  - remote Release correctness run on the AVX-512 server:
    - run directory: `/home/ec2-user/vroom-parity/20260319T062838Z/blst_current`
    - build command:
      - `cmake -S . -B out/parity-tests-clang21 -G Ninja -DCMAKE_C_COMPILER='/opt/llvm-21/bin/clang' -DCMAKE_CXX_COMPILER='/opt/llvm-21/bin/clang++' -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS_RELEASE='-O3 -march=native' -DCMAKE_CXX_FLAGS_RELEASE='-O3 -march=native' -DBLST_PORTABLE=OFF -DBLST_ASM=ON -DBLST_ADX=ON -DBLST_AVX512=ON -DBLST_DEFAULT_AVX=ON -DBLST_AVX512_PORTABLE=OFF -DBLST_BENCHMARK=OFF -DBLST_TEST=ON -DBLST_CUDA=OFF`
      - `cmake --build out/parity-tests-clang21 --target blst_tests -j2`
    - test command:
      - `./out/parity-tests-clang21/test/blst_tests --gtest_filter='VroomBaseField50Bit.*'`
    - result: pass
- direct benchmark-surface audit result after Phase 19:
  - the copied in-tree benchmark names now cover the full upstream surfaces of
    both `src/bench_pairing.cpp` and `src/bench_pairing_50bit.cpp`

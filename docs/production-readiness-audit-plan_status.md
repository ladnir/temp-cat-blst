# Production Readiness Audit Status

## Current State

This audit pass is no longer read-only. The supported C API and the scalar/core
library paths behind it were reviewed, the concrete hazards found in that pass
were remediated in code, and regression coverage was added.

Scope for this pass:

- `bindings/blst_scalar.h`
- `blst/c_api.cpp`
- `blst/exports.cpp`
- `blst/scalar/curve/GT.cpp`
- `blst/scalar/curve/multi_scalar.cpp`
- `blst/scalar/curve/multi_scalar.h`
- `blst/scalar/fields/hash_to_field.h`
- `test/C_api_tests.cpp`

Legacy surfaces intended for removal, such as `bindings/blst.hpp`, were not
treated as primary audit targets except where they affect cleanup scope.

## Findings Closed In This Pass

### 1. Exported integer/Fp/Fr conversion helpers no longer rely on UB type punning

Risk category:

- undefined behavior risk
- public C API correctness risk
- layout and aliasing risk

Root cause:

- low-level conversion entry points in `blst/exports.cpp` reinterpreted
  `uint32_t[]` and `uint64_t[]` buffers as `Fp*` or `Fr*`
- the `uint32_t[12]` path could be validly only 4-byte aligned at the ABI while
  the implementation read it as an `Fp`

Fix landed:

- replaced the typed dereference punning with explicit limb assembly and
  extraction in [blst/exports.cpp](/C:/Users/peter/repo/blst-codex-middle-product/blst/exports.cpp)

Regression coverage:

- added a C API test that round-trips through a buffer that is only
  `uint32_t`-aligned, not `Fp`-aligned, in
  [test/C_api_tests.cpp](/C:/Users/peter/repo/blst-codex-middle-product/test/C_api_tests.cpp)

### 2. Public Pippenger MSM wrappers now define and guard their edge cases

Risk category:

- public API contract gap
- missing defensive checks
- correctness risk

Root cause:

- the wrappers forwarded unsupported inputs into raw bucket code
- `npoints == 1`, `scratch == nullptr`, and `window == 0` were not defended at
  the public boundary even though the header did not forbid them

Fix landed:

- `blst_p1s_mult_pippenger` and `blst_p2s_mult_pippenger` now route
  `npoints == 1` through single-point multiplication
- multi-point calls now reject null scratch explicitly
- tiled wrappers reject zero-window and out-of-range window invocations
- public header docs were tightened in
  [bindings/blst_scalar.h](/C:/Users/peter/repo/blst-codex-middle-product/bindings/blst_scalar.h)

Regression coverage:

- added G1 and G2 tests for single-point fallback, null-scratch rejection, and
  zero-window tile rejection in
  [test/C_api_tests.cpp](/C:/Users/peter/repo/blst-codex-middle-product/test/C_api_tests.cpp)

### 3. Line-precompute null-as-generator behavior is now an explicit supported contract

Risk category:

- stale surface
- public API contract gap

Root cause:

- `blst_precompute_lines` and `blst_miller_loop_lines` already treated null as
  “use the generator,” but the public header did not document that behavior

Fix landed:

- documented null-as-generator semantics in
  [bindings/blst_scalar.h](/C:/Users/peter/repo/blst-codex-middle-product/bindings/blst_scalar.h)

Regression coverage:

- added a C API test that compares null-input and explicit-generator behavior in
  [test/C_api_tests.cpp](/C:/Users/peter/repo/blst-codex-middle-product/test/C_api_tests.cpp)

Decision:

- treat this as accepted compatibility behavior, not as a stale undocumented
  quirk going forward

### 4. Scalar contiguous `miller_loop_n` no longer silently fails above 16 pairs

Risk category:

- correctness bug
- public C API contract gap
- stale internal limit

Root cause:

- the scalar contiguous overload in `GT.cpp` returned early once `n` exceeded
  `MILLER_LOOP_N_MAX`
- the supported C API path routed arbitrary-length input through that overload

Fix landed:

- batched the contiguous scalar overload across repeated chunks instead of
  returning early in
  [blst/scalar/curve/GT.cpp](/C:/Users/peter/repo/blst-codex-middle-product/blst/scalar/curve/GT.cpp)

Regression coverage:

- added a `blst_miller_loop_n` test for `n == 17` against the naive product in
  [test/C_api_tests.cpp](/C:/Users/peter/repo/blst-codex-middle-product/test/C_api_tests.cpp)

### 5. Public XMD helper now has bounded failure handling and fully clears output

Risk category:

- correctness bug
- missing defensive checks
- memory-contract risk

Root cause:

- the helper did not harden the rounded-length arithmetic
- allocation failure was unchecked
- oversize outputs had no defined failure result
- failure cleanup used `vec_zero`, which only clears whole limbs and left
  trailing bytes dirty for non-limb-sized outputs

Fix landed:

- hardened length checks and allocation handling in
  [blst/scalar/fields/hash_to_field.h](/C:/Users/peter/repo/blst-codex-middle-product/blst/scalar/fields/hash_to_field.h)
- documented the failure contract in
  [bindings/blst_scalar.h](/C:/Users/peter/repo/blst-codex-middle-product/bindings/blst_scalar.h)
- replaced failure-path clearing with byte-accurate zeroing for arbitrary output
  lengths

Regression coverage:

- added a non-multiple-of-32 success-path test and an oversize-output
  zeroization test in
  [test/C_api_tests.cpp](/C:/Users/peter/repo/blst-codex-middle-product/test/C_api_tests.cpp)

### 6. Wbits precompute tail handling and size arithmetic are now guarded

Risk category:

- undefined behavior risk
- correctness bug
- missing defensive checks

Root cause:

- the leftover-tail path could continue after chunking had reduced `npoints` to
  zero
- the tail allocation was unchecked
- the public/internal size arithmetic could overflow

Fix landed:

- added zero-tail early returns and allocation checks in
  [blst/scalar/curve/multi_scalar.cpp](/C:/Users/peter/repo/blst-codex-middle-product/blst/scalar/curve/multi_scalar.cpp)
- added overflow-checked size helpers in
  [blst/scalar/curve/multi_scalar.h](/C:/Users/peter/repo/blst-codex-middle-product/blst/scalar/curve/multi_scalar.h)
- propagated the size-helper behavior into the public docs in
  [bindings/blst_scalar.h](/C:/Users/peter/repo/blst-codex-middle-product/bindings/blst_scalar.h)

Regression coverage:

- added a high-window `wbits` MSM regression and an overflow-size negative test
  in [test/C_api_tests.cpp](/C:/Users/peter/repo/blst-codex-middle-product/test/C_api_tests.cpp)

## Validation Added

- C API regression coverage for low-alignment integer/Fp conversion buffers
- C API regression coverage for G1/G2 Pippenger single-point and invalid-input
  behavior
- C API regression coverage for high-window wbits precompute paths
- C API regression coverage for `blst_miller_loop_n` beyond 16 pairs
- C API regression coverage for line-precompute null-as-generator behavior
- C API regression coverage for XMD short-output and oversize-failure handling

## Validation Performed

This remediation pass was validated locally on Windows:

- `cmake --build out/build/x64-Debug --target blst_tests`
- `ctest --test-dir out/build/x64-Debug --output-on-failure`

Outcome:

- Windows debug test run: `47` passed, `0` failed

Previously recorded validation from the earlier audit pass still applies as
separate evidence:

- non-sanitized Linux test run: `240` passed, `3` skipped, `0` failed
- sanitized Linux test run: `240` passed, `3` skipped, `0` failed

The three skipped tests are the AVX expert C API tests in configurations where
that surface is not enabled.

## Accepted Risks / Intentional Behavior Confirmed

- opaque `blst_pairing` storage remains caller-managed and alignment-sensitive
  by design; callers are expected to use `blst_pairing_sizeof()` together with
  `blst_pairing_alignof()`
- `bindings/blst.h` and `bindings/blst_aux.h` remain intentional compatibility
  umbrella headers
- the assert-only packed Fp12 AVX stubs remain an intentionally disabled
  surface during the VROOM reset

## Remaining Open Risks

- decode and deserialize surfaces still need dedicated fuzzing
- Go and C# bindings still need comparable smoke coverage for tightened public
  behavior
- sanitizer coverage is still strongest on Linux scalar/instrumentable paths
  and does not replace native validation of asm- or AVX-heavy configurations
- signedness and narrowing assumptions in arithmetic helpers still merit a
  deeper follow-on audit

## Recommended Next Directions

1. Add fuzz targets for public decode and deserialize surfaces.
2. Add Go and C# smoke coverage for the tightened public contracts.
3. Run a portable scalar Linux Valgrind pass as a secondary dynamic check.
4. Continue a signedness-focused audit in scalar and AVX helper code.

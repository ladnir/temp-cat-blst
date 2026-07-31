# Phase 1 Hardening Evidence

## Record

- Date: 2026-07-30
- Baseline commit: `67a2dfcadf75f5452d503acc41ab3ebbb4568a5d`
- State: validated working-tree change set; not yet committed
- Scope: `AUD-004`, `AUD-007`, `AUD-008`, `AUD-012`, and review-discovered
  `AUD-013`
- Benchmarks: disabled and not run
- CUDA: disabled

This record closes findings, not whole-file reviews. Production and test files
remain pending in `FILE_MANIFEST.csv` until their one-file-at-a-time checklist
records are complete.

## Implemented changes

### AUD-004: empty v5 salt

`blst/sign/keygen.cpp` no longer calls the assembly byte-copy primitive when
the salt length is zero. The public `NULL, 0` rejection contract remains
unchanged. C and both Rust signature variants test a non-null empty salt
against an independently derived exact secret-key vector.

### AUD-007: 256-bit modular right shift

Both right-shift operators in `blst/blst_t.hpp` now use modular right shift.
Assembly-default builds select `rshift_mod_256`; portable builds select the
fixed-size `cc_rshift_mod_n` template at compile time. A zero-bit shift returns
the original value without invoking primitives that require a nonzero count.

The direct tests compare against an independent limb-by-limb modular-halving
oracle over odd, even, cross-limb, repeated, in-place, and non-mutating cases.

### AUD-008: Fr Montgomery operand contract

`Fr::to_mont` now passes the fixed, reduced R-squared constant as the first
operand and the potentially unreduced input as the second operand. This adds
no branches, allocation, or arithmetic and matches the refined assembly
operand contract.

Tests cover modulus boundaries through `2r + 1`, full-width and maximum
256-bit inputs, aliasing, round trips, and 4,096 deterministic wide inputs.

### AUD-012: compressed infinity

New regressions use the exact compressed infinity encodings: 48 bytes for G1
and 96 bytes for G2, with byte zero set to `0xc0`. Coverage reaches scalar C,
native C++, Rust, and portable AVX-512 C wrappers.

### AUD-013: standalone native C++ header and layout

Independent review found that the AUD-007 test initially supplied compatibility
aliases that hid a production defect: `blst_t.hpp` did not compile standalone
after the fork namespaced and renamed its vector types. The header now
explicitly sources its historical global surface from the fork's namespaced
types and helpers.

The `n`, `nbits`, and `degree` constants are now static compile-time metadata.
They no longer inflate each field element or cause methods using
`sizeof(element)` to process non-limb metadata. Compile-time assertions require
the 256- and 384-bit wrapper sizes to equal their limb-array sizes. Independent
re-review found the remediation adequate with no remaining concrete blocker.

## Validation matrix

| Profile | Configuration | Result |
| --- | --- | --- |
| Scalar portable | MSVC Debug; `BLST_ASM=OFF`, `BLST_ADX=OFF`, AVX/CUDA/benchmarks off | Focused Phase 1 executable: 10 passed, AVX-only case skipped |
| Scalar assembly | MSVC Debug; `BLST_ASM=ON`, `BLST_DEFAULT_ASM=ON`, `BLST_ADX=ON`, AVX/CUDA/benchmarks off | Focused Phase 1 executable: 10 passed, AVX-only case skipped; full CTest: 319/319 passed |
| Portable AVX-512 | MSVC Debug; assembly/ADX on, `BLST_AVX512_PORTABLE=ON`, default AVX off, CUDA/benchmarks off | Focused Phase 1 executable: 11/11 passed, including AVX C decoding |
| Rust retained API | `cargo test --features full-rust-bindings` focused filters | Compressed infinity 1/1; empty-salt key generation 2/2 |

Independent lane validation also built and ran the focused CTest registrations
in separate portable and assembly/ADX directories before final integration.
The final header/layout adjustment was rebuilt and rerun in both profiles.

## Residual validation

- Native AVX-512 hardware execution remains for the AVX audit matrix; portable
  AVX-512 executed here.
- Linux GCC/Clang and AArch64 assembly remain for cross-platform CI.
- `keygen_v4_5` shares the guarded implementation but lacks its own explicit
  empty-salt regression.
- The rest of upstream `bf008ec` remains coupled to the assembly generator
  chain and is not closed by the validated Fr change.
- A full portable Debug CTest run exceeded the interactive validation window;
  all focused portable Phase 1 tests passed.
- The changes still need separable commits that exclude the pre-existing BTx
  working-tree edits before the Phase 1 process gate is closed.

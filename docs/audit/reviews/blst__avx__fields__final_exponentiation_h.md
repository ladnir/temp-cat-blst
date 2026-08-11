# `blst/avx/fields/final_exponentiation.h` Review

## Scope

This record covers the complete selected AVX final-exponentiation owner at
commit `835ffbb`. It traces the easy part, fixed hard-exponent chain, inversion
boundary, production C wrapper, scalar BLST reference schedule, and accepted
Lean correspondence. BTx/BTE and CUDA are excluded.

## Contracts

- `final_exp` accepts a valid nonzero Fp12 value. It is allocation-free and
  fixed-shape after inversion, but inherits the selected inverter's
  variable-time behavior. It performs no nonzero or field-membership check in
  Release.
- `final_exp_after_inverse` additionally requires `f_inverse` to be the exact
  multiplicative inverse of `f`. The two inputs are read-only complete values.
- The result is BLST's established cubed pairing final exponent,
  `f^(3 * ((p^12 - 1) / r))`, and therefore lies in GT for a valid nonzero
  input.
- `raise_to_z` and `raise_to_z_div_by_2` accept cyclotomic/GT values and execute
  fixed public addition chains. They perform no membership validation.
- The current returning helpers permit exact caller-side storage reuse only at
  the assignment boundary. Their internal input and output owners must remain
  distinct while a cyclotomic square is evaluated because both split halves
  read the complete input.

## Correctness result

No arithmetic defect is confirmed. The source schedule matches scalar BLST's
`Fp12::final_exp`, including the inherited cubed exponent. Accepted
`AVX-P3-005` proves the literal easy/hard schedule and exact operation ledger;
accepted `AVX-P3-008` proves the concrete inversion lift assuming one correct
nonzero base reciprocal.

## Findings

- `AUD-274`: every repeated cyclotomic square materializes and copies a full
  1,536-byte Fp12 owner in optimized native builds.
- `AUD-275`: runtime-count helpers and SSA-like stage ownership obscure a
  compile-time-fixed schedule and retain 35--38 KiB `final_exp_after_inverse`
  stack frames.
- `AUD-276`: provenance and the nonzero, inverse, cyclotomic, timing, and public
  C contracts are incomplete.

## Code-generation evidence

Compile-only inspection used native AVX-512/IFMA configurations; no executable
or benchmark ran.

- MSVC `/O2 /DNDEBUG /arch:AVX512` gives `mul_n_sqr` a 1,664-byte frame and a
  1,536-byte `memcpy` per square iteration. The complete five-chain schedule
  executes at least 305 such copies, totaling 468,480 bytes before stage-return
  copies. `final_exp_after_inverse` has a 38,528-byte frame and eight direct
  `memcpy` calls.
- Clang 18 `-O3 -DNDEBUG -mavx512f -mavx512ifma` independently retains the
  defect. `final_exp_after_inverse` has a 35,456-byte frame and five direct
  `memcpy` calls; specialized raise-chain bodies retain repeated full-element
  copies.

## Proposed remediation

Introduce a distinct-output cyclotomic-square boundary and use two fixed
ping-pong owners for repeated squares. Then express the hard part with a small,
reviewed mutable accumulator set and fixed/compile-time square counts. Preserve
the exact formulas, multiplication/square ledger, allocation-free behavior,
and lack of runtime dispatch. Reinspect MSVC and Clang code generation before
requesting a single serial Peach benchmark.

## Residual risk

The file remains `in_review` until the ownership repair, focused portable-Debug
differentials, exact Lean/source rebind, and optimized native code-generation
inspection pass. P0 compiler/IFMA semantics, sanitizer/fuzz, Apple, and broad
release evidence remain external milestone work.

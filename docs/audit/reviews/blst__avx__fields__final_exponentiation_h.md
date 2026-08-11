# `blst/avx/fields/final_exponentiation.h` Review

## Scope

This record covers the complete selected AVX final-exponentiation owner through
commit `35516ff`. It traces the easy part, fixed hard-exponent chain, inversion
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

- `AUD-274`: repeated cyclotomic squares materialized full 1,536-byte copies.
  The fixed chains now use one local distinct-output boundary and two ping-pong
  owners; exact in-place split evaluation remains explicitly forbidden.
- `AUD-275`: runtime counts and SSA-like ownership retained 35--38 KiB frames.
  Square counts are now compile-time constants and the hard part has a bounded
  mutable accumulator set behind deliberate noinline long-chain kernels.
- `AUD-276`: provenance and the nonzero, inverse, cyclotomic, timing, and public
  C contracts are now explicit.

## Code-generation evidence

Compile-only inspection used native AVX-512/IFMA configurations; no executable
or benchmark ran.

- Before repair, MSVC emitted at least 305 loop-level copies totaling 468,480
  bytes and used a 38,528-byte hard-part frame. Clang used 35,456 bytes and
  independently retained the repeated copies.
- At `35516ff`, MSVC gives `final_exp_after_inverse` a 12,416-byte frame. Its
  largest nested multiply/square kernel uses 3,192 bytes, for about 15.6 KiB
  peak local stack. Each raise kernel has zero `memcpy`; its five product
  boundaries retain one copy each, so the five-chain final exponent has 25
  nested product copies rather than at least 305 square-loop copies.
- Clang 18 gives the repaired hard part a 14,720-byte frame. Its raise kernel
  likewise has zero `memcpy`; the ownership result is cross-toolchain.

## Remediation

The distinct-output square is private to `final_exp_detail`, so the already
proved general Fp12 owner and API remain unchanged. Two ping-pong owners remove
loop copies. Fixed template counts replace runtime counts, and noinline long
chains prevent compiler inlining from inflating the hard-part frame. The hard
schedule uses four wide and four standard owners while preserving the exact
formula and call ledger.

## Residual risk

Five focused portable-Debug Fp12/final-exponentiation/Miller differentials pass,
including direct full/half fixed-z comparisons. The three focused Lean modules
build, and `AVX-P3-005` is rebound to the exact source hash while retaining its
35 multiplications, 315 cyclotomic squares, nine conjugations, Frobenius counts
1/2/1, and one inversion. The audit file is complete. A serial native benchmark
would quantify speedup but is not required to establish removal of the copies
and oversized stack. P0 compiler/IFMA semantics, sanitizer/fuzz, Apple, and
broad release evidence remain external milestone work.
